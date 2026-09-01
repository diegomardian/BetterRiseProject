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
cp sections/*.tex "$OUT/sections/"
cp figures/*.pdf "$OUT/figures/"

cat > "$OUT/READ_ME_FIRST.txt" <<'EOF'
WMHS @ NeurIPS 2026 — Overleaf project
======================================

BEFORE IT WILL COMPILE
  Add neurips_2026.sty to this project's root. It is the official workshop
  style and is not bundled here. Do not modify it — the style file warns that
  tweaking it risks desk rejection.

TWO DOCUMENTS, ONE SOURCE
  main.tex        full paper        6 pages of main text (limit 9)
  main_short.tex  extended abstract 4 pages of main text (limit 4)

  Both \input the same sections/ files and differ only in the \iffull flag on
  line 33. The short version is a strict subset of the same prose, so a number
  cannot differ between them. Pick which one to compile in Overleaf under
  Menu > Settings > Main document.

COMPILER
  pdfLaTeX. Overleaf runs bibtex automatically; if references show as [?],
  recompile once more.

BEFORE YOU SUBMIT
  §5 contains a red [ANONYMISED-ARTIFACT-URL] placeholder. Replace
  \artifacturl in main.tex (and main_short.tex) with the anonymised mirror
  URL. It is red so it cannot ship unnoticed.
EOF

if command -v zip >/dev/null 2>&1; then
  ( cd "$OUT" && zip -qr "../$OUT.zip" . )
  echo "wrote $OUT.zip  ($(find "$OUT" -type f | wc -l | tr -d ' ') files, $(du -h "$OUT.zip" | cut -f1))"
else
  echo "wrote $OUT/ (zip not installed — compress it yourself)"
fi
