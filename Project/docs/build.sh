#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python "$ROOT/scripts/generate_all.py"
cd "$ROOT/docs"
if command -v bibtex >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
elif [[ -x /usr/bin/bibtex.original ]]; then
  latexmk -e '$bibtex="/usr/bin/bibtex.original %O %B"' \
    -pdf -interaction=nonstopmode -halt-on-error main.tex
else
  echo "A working BibTeX executable is required." >&2
  exit 1
fi
