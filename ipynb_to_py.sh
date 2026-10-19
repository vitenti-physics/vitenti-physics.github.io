#!/usr/bin/env bash
# Post-render: export every rendered notebook to a plain Python script.
set -e

out_dir="${QUARTO_PROJECT_OUTPUT_DIR:-_site}"

find "$out_dir" -name "*.ipynb" -print0 |
  while IFS= read -r -d '' nb; do
    jupyter nbconvert --to python "$nb" --log-level WARN
  done
