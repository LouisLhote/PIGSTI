#!/usr/bin/env bash
# Render Snakemake DAG to SVG (Graphviz dot required).
# Usage: from PIGSTI_publication/ :  bash scripts/render_snakemake_dag.sh [output.svg]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

OUT="${1:-docs/images/snakemake_dag.svg}"
mkdir -p "$(dirname "${OUT}")"

if ! command -v dot >/dev/null 2>&1; then
  echo "[render_snakemake_dag] ERROR: dot (Graphviz) not found. Install graphviz package." >&2
  exit 1
fi

# --dag prints DOT to stdout. Optionally bound the graph, e.g.
#   snakemake --snakefile Snakefile --configfile config/config.yaml --cores 1 --dag results/workflow/ \
#   | dot -Tsvg -o docs/images/snakemake_dag.svg
snakemake \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --cores 1 \
  --use-conda \
  --conda-prefix "${SNAKEMAKE_CONDA_PREFIX:-.snakemake/conda}" \
  --dag 2>/dev/null \
  | dot -Tsvg -o "${OUT}"

echo "[render_snakemake_dag] Wrote ${OUT}"
