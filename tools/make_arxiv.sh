#!/usr/bin/env bash
# Build a self-contained arXiv source bundle (paper/arxiv.tar.gz) and check that it compiles on its own.
# Usage: tools/make_arxiv.sh   (run after `make all`)
set -euo pipefail
cd "$(dirname "$0")/.."

out=$(mktemp -d)
cp paper/main.tex paper/numbers.tex paper/tab_hypotheses.tex paper/tab_robustness.tex "$out"/
cp figures/fig1_heat_timeline.png figures/fig2_event_study_dev2.png "$out"/
# figures sit next to main.tex in the bundle
sed -i 's#\\graphicspath{{../figures/}}#\\graphicspath{{./}}#' "$out/main.tex"

# compile check in a clean copy
check=$(mktemp -d)
cp "$out"/* "$check"/
(cd "$check" && pdflatex -interaction=nonstopmode main.tex >/dev/null && pdflatex -interaction=nonstopmode main.tex >/dev/null)
if grep -E "^!|undefined" "$check/main.log"; then
  echo "bundle does not compile cleanly" >&2
  exit 1
fi

tar -czf paper/arxiv.tar.gz -C "$out" .
echo "wrote paper/arxiv.tar.gz:"
tar -tzf paper/arxiv.tar.gz
rm -rf "$out" "$check"
