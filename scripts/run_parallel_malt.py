#!/usr/bin/env python3
"""
Parallel MALT on a single node: warm the OS page cache once, then run malt-run
processes concurrently with --memoryMode map so they share the cached index.

Designed for HOPS pipelines where per-sample ``hops -m malt`` reloads the full
index into each JVM (memoryMode=load). After all RMA files are written, HOPS
``-m me_po`` runs separately (Snakemake rule hops_maltex_post).
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

CHUNK_BYTES = 64 * 1024 * 1024


def _log(log_path: Path, message: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    line = message.rstrip()
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def _parse_hops_config(path: Path) -> dict[str, str]:
    cfg: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        cfg[key.strip()] = value.strip()
    return cfg


def _index_target(index_prefix: str) -> Path:
    """Return the filesystem path to warm (directory or file prefix path)."""
    path = Path(index_prefix)
    if path.is_dir():
        return path
    if path.is_file():
        return path
    if path.parent.is_dir() and any(path.parent.glob(f"{path.name}*")):
        return path
    raise FileNotFoundError(f"MALT index not found: {index_prefix}")


def _index_files(index_prefix: str) -> list[Path]:
    target = _index_target(index_prefix)
    if target.is_file():
        return [target]
    if target.is_dir():
        files = sorted(p for p in target.rglob("*") if p.is_file())
        if files:
            return files
        raise FileNotFoundError(f"No files in MALT index directory: {target}")
    parent = target.parent
    stem = target.name
    files = sorted(p for p in parent.glob(f"{stem}*") if p.is_file())
    if not files:
        raise FileNotFoundError(
            f"No MALT index files found for prefix: {index_prefix} (looked in {parent})"
        )
    return files


def _warm_index(index_prefix: str, log_path: Path) -> float:
    start = time.monotonic()
    warm_path = _index_target(index_prefix)
    files = _index_files(index_prefix)
    if shutil.which("vmtouch"):
        cmd = ["vmtouch", "-t", str(warm_path)]
        _log(log_path, f"[warmup] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    else:
        _log(
            log_path,
            f"[warmup] vmtouch not found; reading {len(files)} index file(s) into page cache",
        )
        for path in files:
            _log(log_path, f"[warmup] read {path}")
            with path.open("rb") as fh:
                while fh.read(CHUNK_BYTES):
                    pass
    elapsed = time.monotonic() - start
    _log(log_path, f"[warmup] finished in {elapsed:.1f}s")
    return elapsed


def _expected_rma(bio: str, malt_root: Path) -> Path:
    return malt_root / bio / "malt" / f"{bio}_unaligned.rma6"


def _clear_sample_output(malt_root: Path, bio: str) -> None:
    target = malt_root / bio
    if target.is_symlink() or target.is_file():
        target.unlink()
    elif target.is_dir():
        shutil.rmtree(target)


def _normalize_rma(bio: str, malt_root: Path, log_path: Path | None = None) -> Path | None:
    """Find malt-run output and move it to the HOPS-compatible path if needed."""
    expected = _expected_rma(bio, malt_root)
    if expected.is_file() and expected.stat().st_size > 0:
        return expected

    out_path = malt_root / bio
    found: Path | None = None

    if out_path.is_file() and out_path.stat().st_size > 0:
        # malt-run -o <path> wrote the RMA as a plain file (no .rma6 suffix).
        found = out_path
    elif out_path.is_dir():
        preferred = out_path / "malt" / f"{bio}_unaligned.rma6"
        if preferred.is_file() and preferred.stat().st_size > 0:
            return preferred
        rmas = [p for p in out_path.rglob("*.rma6") if p.is_file() and p.stat().st_size > 0]
        if rmas:
            found = max(rmas, key=lambda p: p.stat().st_size)

    if found is None:
        return None

    if found.resolve() == expected.resolve():
        return expected

    expected.parent.mkdir(parents=True, exist_ok=True)
    if expected.exists():
        expected.unlink()
    shutil.move(str(found), str(expected))
    if log_path is not None:
        _log(log_path, f"[malt] {bio}: normalized RMA {found} -> {expected}")
    return expected


def _malt_run_cmd(
    *,
    malt_bin: str,
    index: str,
    fastq: Path,
    output_dir: Path,
    threads: int,
    heap_gb: int,
    memory_mode: str,
    alignment: str,
    align_type: str,
    percent_id: str,
) -> list[str]:
    return [
        malt_bin,
        "-d",
        index,
        "-i",
        str(fastq),
        "-o",
        str(output_dir),
        "-t",
        str(threads),
        f"-J-Xmx{heap_gb}G",
        "--memoryMode",
        memory_mode,
        "-m",
        alignment,
        "-at",
        align_type,
        "-id",
        percent_id,
    ]


def _run_one(
    *,
    bio: str,
    fastq: Path,
    malt_root: Path,
    malt_bin: str,
    index: str,
    threads: int,
    heap_gb: int,
    memory_mode: str,
    alignment: str,
    align_type: str,
    percent_id: str,
    resume: bool,
    log_path: Path,
) -> tuple[str, int, float, str]:
    out_dir = malt_root / bio
    start = time.monotonic()

    if resume:
        existing = _normalize_rma(bio, malt_root, log_path)
        if existing is not None:
            elapsed = time.monotonic() - start
            _log(log_path, f"[malt] {bio}: reuse existing RMA ({existing}) in {elapsed:.1f}s")
            return bio, 0, elapsed, "reused"

    _clear_sample_output(malt_root, bio)
    # HOPS passes an existing output directory; malt-run then writes malt/<sample>_unaligned.rma6 inside.
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = _malt_run_cmd(
        malt_bin=malt_bin,
        index=index,
        fastq=fastq,
        output_dir=out_dir,
        threads=threads,
        heap_gb=heap_gb,
        memory_mode=memory_mode,
        alignment=alignment,
        align_type=align_type,
        percent_id=percent_id,
    )
    hops_log_dir = log_path.parent / "hops_malt"
    hops_log_dir.mkdir(parents=True, exist_ok=True)
    sample_log = hops_log_dir / f"{bio}.log"
    _log(log_path, f"[malt] {bio}: start {' '.join(cmd)}")
    with sample_log.open("w", encoding="utf-8") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    elapsed = time.monotonic() - start

    if proc.returncode != 0:
        _log(log_path, f"[malt] {bio}: FAILED exit={proc.returncode} in {elapsed:.1f}s")
        return bio, proc.returncode, elapsed, "failed"

    rma = _normalize_rma(bio, malt_root, log_path)
    if rma is None:
        _log(log_path, f"[malt] {bio}: missing RMA under {out_dir} after exit 0")
        return bio, 1, elapsed, "missing_rma"

    _log(log_path, f"[malt] {bio}: OK in {elapsed:.1f}s -> {rma}")
    return bio, 0, elapsed, "ok"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bio-samples", nargs="+", required=True)
    parser.add_argument("--fastq", nargs="+", required=True, help="One FASTQ per bio sample")
    parser.add_argument("--malt-root", required=True)
    parser.add_argument("--config", required=True, help="config_hops_custom.txt for MALT params")
    parser.add_argument("--parallel", type=int, default=2, help="Concurrent malt-run jobs")
    parser.add_argument("--threads", type=int, default=10, help="Threads per malt-run")
    parser.add_argument("--heap-gb", type=int, default=400, help="Java -Xmx per malt-run (GB)")
    parser.add_argument("--memory-mode", default="map", choices=["map", "page", "load"])
    parser.add_argument("--malt-bin", default="malt-run")
    parser.add_argument("--skip-warmup", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Skip samples with valid RMA")
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    if len(args.bio_samples) != len(args.fastq):
        print("ERROR: --bio-samples and --fastq must have the same length", file=sys.stderr)
        return 1

    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")

    malt_root = Path(args.malt_root)
    malt_root.mkdir(parents=True, exist_ok=True)

    hops_cfg = _parse_hops_config(Path(args.config))
    index = hops_cfg.get("index")
    if not index:
        print("ERROR: index= missing from HOPS config", file=sys.stderr)
        return 1

    alignment = hops_cfg.get("m", "BlastN")
    align_type = hops_cfg.get("at", "SemiGlobal")
    percent_id = hops_cfg.get("id", "90.0")

    total_threads = args.parallel * args.threads
    nproc = os.cpu_count() or 1
    if total_threads > nproc:
        _log(
            log_path,
            f"[warn] parallel×threads = {args.parallel}×{args.threads} = {total_threads} "
            f"> CPU count ({nproc}); jobs may contend",
        )

    malt_bin = args.malt_bin
    if shutil.which(malt_bin) is None:
        print(f"ERROR: {malt_bin} not found in PATH", file=sys.stderr)
        return 1

    pipeline_start = time.monotonic()
    warmup_s = 0.0
    if args.skip_warmup:
        _log(log_path, "[warmup] skipped (--skip-warmup)")
    else:
        warmup_s = _warm_index(index, log_path)

    align_start = time.monotonic()
    jobs = list(zip(args.bio_samples, [Path(p) for p in args.fastq]))
    results: list[tuple[str, int, float, str]] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as pool:
        futures = {
            pool.submit(
                _run_one,
                bio=bio,
                fastq=fq,
                malt_root=malt_root,
                malt_bin=malt_bin,
                index=index,
                threads=args.threads,
                heap_gb=args.heap_gb,
                memory_mode=args.memory_mode,
                alignment=alignment,
                align_type=align_type,
                percent_id=percent_id,
                resume=args.resume,
                log_path=log_path,
            ): bio
            for bio, fq in jobs
        }
        for fut in concurrent.futures.as_completed(futures):
            results.append(fut.result())

    align_s = time.monotonic() - align_start
    total_s = time.monotonic() - pipeline_start

    failed = [bio for bio, rc, _, _ in results if rc != 0]
    ok = [bio for bio, rc, _, _ in results if rc == 0]

    for bio, _, elapsed, status in sorted(results, key=lambda x: x[0]):
        _log(log_path, f"[summary] {bio}: {status} ({elapsed:.1f}s)")

    _log(
        log_path,
        f"[summary] warmup={warmup_s:.1f}s align={align_s:.1f}s total={total_s:.1f}s "
        f"ok={len(ok)} failed={len(failed)}",
    )

    for bio in ok:
        done = malt_root / bio / ".malt_done"
        done.parent.mkdir(parents=True, exist_ok=True)
        done.touch()

    if failed:
        print(
            "ERROR: parallel MALT failed for: " + ", ".join(sorted(failed)),
            file=sys.stderr,
        )
        print(f"See {log_path} and logs/hops_malt/<sample>.log", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
