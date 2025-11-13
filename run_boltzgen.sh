#!/usr/bin/env bash

#
# Run BoltzGen on a design specification YAML.

# # Protein binders
# ./run_boltzgen.sh protein_vs_2vsm.yaml protein-anything

# # Nanobody binders
# ./run_boltzgen.sh nanobody_vs_2vsm.yaml nanobody-anything
#


set -euo pipefail

run_in() {
  local env="$1"; shift
  if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n "$env" "$@"
  elif command -v conda >/dev/null 2>&1; then
    conda run -n "$env" "$@"
  else
    echo "ERROR: neither micromamba nor conda found in PATH." >&2
    exit 1
  fi
}

YAML="${1:-}"
if [[ -z "$YAML" || ! -f "$YAML" ]]; then
  echo "Usage: $0 <design_spec.yaml>  [protocol]  [out_dir]"
  echo "  protocol: protein-anything | nanobody-anything (default: protein-anything)"
  exit 2
fi

PROTOCOL="${2:-nanobody-anything}"          # or: nanobody-anything/protein-anything
OUT="${3:-./boltzgen_out_10000_prot/$(basename "$YAML" .yaml)_${PROTOCOL}}"
NUM_DESIGNS="${NUM_DESIGNS:-10000}"           # sanity check run; scale to 10k–60k later
BUDGET="${BUDGET:-10}"
DEVICES="${DEVICES:-1}"                    # number of GPUs
BG_ENV="${RFD_ENV:-bg}"

# --- Hugging Face offline mode ---
# export HF_HUB_OFFLINE=1
# export TRANSFORMERS_OFFLINE=1
# export HF_DATASETS_OFFLINE=1
# export HF_HUB_ENABLE_HF_TRANSFER=0
# # Optional: ensure a consistent local cache (change if you keep models elsewhere)
# : "${HF_HOME:=$HOME/.cache/huggingface}"
# export HF_HOME
# ---------------------------------

# Optional: control cache location (first run downloads ~6 GB of weights)
# export HF_HOME=/path/to/cache

echo "[check] $YAML"
run_in "$BG_ENV" boltzgen check "$YAML"                     # quick mmCIF + binding-site sanity check. :contentReference[oaicite:2]{index=2}

echo "[run] $PROTOCOL  ->  $OUT"
run_in "$BG_ENV" boltzgen run "$YAML" \
  --output "$OUT" \
  --protocol "$PROTOCOL" \
  --num_designs "$NUM_DESIGNS" \
  --budget "$BUDGET" \
  --devices "$DEVICES" \
  --use_kernels false \
  --reuse


echo
echo "[done] Designs in: $OUT"
echo "Tip: fast re-filtering (tune thresholds without re-designing):"
echo "boltzgen run \"$YAML\" --output \"$OUT\" --protocol \"$PROTOCOL\" --steps filtering --alpha 0.2"
