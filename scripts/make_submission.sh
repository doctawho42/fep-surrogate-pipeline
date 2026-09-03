#!/usr/bin/env bash
# Assemble the ACS Paragon Plus submission set.
#
# ACS accepts a LaTeX submission as two items: an author-generated PDF, and a ZIP archive
# holding every source the manuscript references (the other .tex files, the .bib, the graphics).
# This builds both, plus the Supporting Information PDF and the cover letter.
#
# The archive is verified by compiling it in a scratch directory before it is packed. That check
# is not ceremony: the manuscript and the supplement cross-reference each other through xr, and a
# file left out of the archive shows up as a stale number rather than as an error.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/submission"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

mkdir -p "$OUT" "$STAGE/figures"
cd "$ROOT/docs"
cp paper_jctc.tex paper_jctc_si.tex paper_body.tex paper_si.tex refs.bib \
   tab_auditability.tex tab_design.tex tab_ground.tex tab_metrics.tex "$STAGE/"
for g in $(grep -oh "includegraphics\[[^]]*\]{[^}]*}" paper_body.tex paper_si.tex paper_jctc.tex \
           | sed 's/.*{//;s/}//' | sort -u); do
  cp "$ROOT/figs/$g" "$STAGE/figures/$g"
done
sed -i.bak 's|\\graphicspath{{\.\./figs/}}|\\graphicspath{{figures/}}|' \
  "$STAGE/paper_jctc.tex" "$STAGE/paper_jctc_si.tex"
rm -f "$STAGE"/*.bak
cp "$ROOT/docs/submission_readme.txt" "$STAGE/README.txt"

cd "$STAGE"
for _ in 1 2; do
  pdflatex -interaction=nonstopmode paper_jctc    >/dev/null 2>&1 || true
  bibtex   paper_jctc                             >/dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode paper_jctc_si >/dev/null 2>&1 || true
  bibtex   paper_jctc_si                          >/dev/null 2>&1 || true
done
for _ in 1 2 3; do
  pdflatex -interaction=nonstopmode paper_jctc    >/dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode paper_jctc_si >/dev/null 2>&1 || true
done
for f in paper_jctc paper_jctc_si; do
  test -f "$f.pdf" || { echo "FAIL: $f did not compile from the archive"; exit 1; }
  n=$(pdftotext "$f.pdf" - | grep -c '??' || true)
  test "$n" -eq 0 || { echo "FAIL: $f carries $n unresolved reference(s) when built from the archive"; exit 1; }
done
echo "archive compiles clean in isolation"

# keep sources only; every build product goes
rm -f ./*.aux ./*.bbl ./*.blg ./*.out ./*.pdf ./*.fls ./*.fdb_latexmk ./*.toc ./*.lo? ./acs-*.bib
rm -f "$OUT/manuscript_latex.zip"
zip -rq "$OUT/manuscript_latex.zip" . -x '.*'

cp "$ROOT/docs/paper_jctc.pdf"    "$OUT/manuscript.pdf"
cp "$ROOT/docs/paper_jctc_si.pdf" "$OUT/supporting_information.pdf"
cp "$ROOT/docs/cover_letter.pdf"  "$OUT/cover_letter.pdf"
echo "wrote $OUT: manuscript.pdf, supporting_information.pdf, manuscript_latex.zip, cover_letter.pdf"
