#!/usr/bin/env python3
"""Stage per-sample MALT .rma6 files into a shared directory for HOPS MaltExtract."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bio-samples", nargs="+", required=True)
    parser.add_argument("--malt-root", default="results/metagenomics/hops/malt")
    parser.add_argument("--rma-root", default="results/metagenomics/hops/rma")
    parser.add_argument("--output", nargs="+", required=True, help="Expected staged RMA paths")
    args = parser.parse_args()

    malt_root = Path(args.malt_root)
    rma_root = Path(args.rma_root)
    rma_root.mkdir(parents=True, exist_ok=True)

    for bio in args.bio_samples:
        candidates = [
            malt_root / bio / "malt" / f"{bio}_unaligned.rma6",
            malt_root / bio / f"{bio}_unaligned.rma6",
        ]
        src = next((p for p in candidates if p.is_file()), None)
        dst = rma_root / f"{bio}_unaligned.rma6"
        if src is None:
            print(f"ERROR: missing MALT output for {bio}; checked:", file=sys.stderr)
            for p in candidates:
                print(f"  - {p}", file=sys.stderr)
            return 1
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src.resolve())

    missing = [p for p in args.output if not Path(p).exists()]
    if missing:
        for p in missing:
            print(f"ERROR: staged RMA not created: {p}", file=sys.stderr)
        return 1

    print(f"Staged {len(args.bio_samples)} RMA files under {rma_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
