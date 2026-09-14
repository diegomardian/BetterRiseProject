#!/usr/bin/env bash
# Bundle exactly what Overleaf needs, and nothing else.
#
# Overleaf compiles LaTeX. It has no result tables, no Python and no repo, so the
# figure scripts and the shell guards stay here — the figures go in as the built
# PDFs they already are. Regenerate those locally (make_fig1.py, make_fig3.py)
# and re-run this whenever a number changes.
#
# neurips_2026.sty is NOT included: it is not in this repository, deliberately.
# Download the official workshop style and add it in Overleaf.
set -euo pipefail
cd "$(dirname "$0")"
OUT="${1:-overleaf}"
rm -rf "$OUT" "$OUT.zip"
mkdir -p "$OUT/sections" "$OUT/figures"

cp main.tex main_short.tex refs.bib "$OUT/"
cp -R sections/. "$OUT/sections/"
cp -R figures/. "$OUT/figures/"

cat > "$OUT/READ_ME_FIRST.txt" <<'EOF'
WMHS @ NeurIPS 2026 — Overleaf project
======================================

BEFORE IT WILL COMPILE
  Add neurips_2026.sty to this project's root. It is the official workshop
  style and is not bundled here. Do not modify it — the style file warns that
  tweaking it risks desk rejection.

TWO DOCUMENTS
  main.tex        revised full paper (main-text limit 9 pages)
  main_short.tex  original extended abstract (main-text limit 4 pages)

  The revised full paper reads sections/full/, including its condensed appendix.
  The extended abstract retains sections/. Figures and bibliography are shared.
  Select the document under Menu > Settings > Main document.

COMPILER
  pdfLaTeX. Overleaf runs bibtex automatically; if references show as [?],
  recompile once more.

BEFORE YOU SUBMIT
  Use the official workshop style, rebuild, and check the page limit and
  anonymity. The responsible-use statement is in section 5. No artifact URL
  is included in this submission.

EOF

if command -v zip >/dev/null 2>&1; then
  ( cd "$OUT" && zip -qr "../$OUT.zip" . )
  echo "wrote $OUT.zip  ($(find "$OUT" -type f | wc -l | tr -d ' ') files, $(du -h "$OUT.zip" | cut -f1))"
else
  echo "wrote $OUT/ (zip not installed — compress it yourself)"
fi
