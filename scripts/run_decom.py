#!/usr/bin/env python3
"""Run decOM with Dask tuning, automatic low-resource retry, and clear failure modes."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _append_log(log_path: Path, message: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(message.rstrip() + "\n")


def _dask_env() -> dict[str, str]:
    env = os.environ.copy()
    # Workers killed by OOM often leave Dask clients hanging on shutdown.
    for key, value in {
        "DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT": "300s",
        "DASK_DISTRIBUTED__COMM__TIMEOUTS__TCP": "300s",
        "DASK_DISTRIBUTED__COMM__TIMEOUTS__COMM": "300s",
        "DASK_DISTRIBUTED__COMM__TIMEOUTS__DEATH": "300s",
        "DASK_DISTRIBUTED__ADMIN__LOG_LENGTH": "0",
    }.items():
        env.setdefault(key, value)
    return env


def _run_once(
    *,
    p_sink: Path,
    p_sources: Path,
    p_keys_dir: Path,
    memory: str,
    threads: int,
    output_dir: Path,
    log_path: Path,
    label: str,
) -> int:
    keys_arg = str(p_keys_dir)
    if not keys_arg.endswith(os.sep):
        keys_arg += os.sep

    cmd = [
        "decOM",
        "-p_sinks",
        str(p_sink),
        "-p_sources",
        str(p_sources),
        "-p_keys",
        keys_arg,
        "-mem",
        memory,
        "-t",
        str(threads),
        "-o",
        str(output_dir),
    ]
    _append_log(
        log_path,
        f"[run_decom] {label}: mem={memory} threads={threads} "
        f"(approx RAM budget = mem × threads)\n"
        f"[run_decom] command: {' '.join(cmd)}",
    )
    with log_path.open("a", encoding="utf-8") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=_dask_env())
    return int(proc.returncode)


def _attempts(memory: str, threads: int) -> list[tuple[str, int, str]]:
    """Retry plan: keep configured memory; only reduce threads on retry."""
    plans: list[tuple[str, int, str]] = [(memory, threads, "primary")]
    if threads > 1:
        plans.append((memory, 1, "retry_single_thread"))
    seen: set[tuple[str, int]] = set()
    unique: list[tuple[str, int, str]] = []
    for mem, thr, label in plans:
        key = (mem, thr)
        if key in seen:
            continue
        seen.add(key)
        unique.append((mem, thr, label))
    return unique


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-sink", required=True)
    parser.add_argument("--p-sources", required=True)
    parser.add_argument("--p-keys-dir", required=True)
    parser.add_argument("--memory", default="64GB")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument(
        "--fail-on-error",
        action="store_true",
        help="Exit non-zero when decOM cannot run or all attempts fail (stops Snakemake).",
    )
    args = parser.parse_args()

    p_sink = Path(args.p_sink)
    p_sources = Path(args.p_sources)
    p_keys_dir = Path(args.p_keys_dir)
    output_dir = Path(args.output)
    log_path = Path(args.log)

    def _soft_fail(reason: str) -> int:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "FAILED.txt").write_text(
            f"{reason}\nSee {log_path}\n", encoding="utf-8"
        )
        _append_log(log_path, f"[run_decom] {reason}")
        if log_path.is_file():
            tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-30:]
            for line in tail:
                print(line, file=sys.stderr)
        if args.fail_on_error:
            print(
                "ERROR: decOM failed after all retry attempts. "
                "Set decom_fail_on_error: false to continue the pipeline, "
                "or enable_decom: false to skip decOM.",
                file=sys.stderr,
            )
            return 1
        (output_dir / ".decom_skipped").touch()
        _append_log(log_path, "[run_decom] continuing (decom_fail_on_error=false)")
        return 0

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")

    # Validate inputs (same checks as validate_decom_inputs.py).
    validate_script = Path(__file__).resolve().parent / "validate_decom_inputs.py"
    if validate_script.is_file():
        proc = subprocess.run(
            [
                sys.executable,
                str(validate_script),
                "--p-sink",
                str(p_sink),
                "--p-keys-dir",
                str(p_keys_dir),
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "").strip()
            return _soft_fail(f"decOM input validation failed: {detail or 'see log'}")

    if not p_sources.is_dir():
        return _soft_fail(f"decOM_sources is missing or not a directory: {p_sources}")

    last_rc = 1
    for mem, threads, label in _attempts(args.memory, args.threads):
        if output_dir.exists():
            shutil.rmtree(output_dir)
        # decOM refuses to start if -o already exists; do not pre-create the directory.

        last_rc = _run_once(
            p_sink=p_sink,
            p_sources=p_sources,
            p_keys_dir=p_keys_dir,
            memory=mem,
            threads=threads,
            output_dir=output_dir,
            log_path=log_path,
            label=label,
        )
        if last_rc == 0:
            if not output_dir.is_dir():
                return _soft_fail(f"decOM exited 0 but output directory missing: {output_dir}")
            (output_dir / ".decom_done").touch()
            _append_log(log_path, "[run_decom] decOM completed successfully")
            return 0

        _append_log(log_path, f"[run_decom] {label} failed with exit code {last_rc}")

    return _soft_fail("decOM failed after all retry attempts")


if __name__ == "__main__":
    raise SystemExit(main())
