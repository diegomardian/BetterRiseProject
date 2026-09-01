#!/usr/bin/env bash
# Bundle exactly what Overleaf needs for the ICBINB submission, and nothing else.
#
# Overleaf compiles LaTeX. It has no result tables and no Python, so make_figs.py
# and the shell guards stay here — the figures go in as the built PDFs they
# already are. Regenerate them locally and re-run this whenever a number changes.
#
# neurips_2026.sty is NOT included: it is not in this repository, deliberately.
# Download the official workshop style and add it in Overleaf.
set -euo pipefail
cd "$(dirname "$0")"
OUT="${1:-overleaf}"
rm -rf "$OUT" "$OUT.zip"
mkdir -p "$OUT/sections" "$OUT/figures"

cp main.tex refs.bib "$OUT/"
cp sections/*.tex "$OUT/sections/"
cp figures/*.pdf "$OUT/figures/"

cat > "$OUT/READ_ME_FIRST.txt" <<'EOF'
ICBINB-BIO @ NeurIPS 2026 — Overleaf project
============================================

DEADLINE: 2 September 2026, 11:59pm AoE.

BEFORE IT WILL COMPILE
  Add neurips_2026.sty to this project's root. It is the official workshop
  style and is not bundled here. Do not modify it — the style file warns that
  tweaking it risks desk rejection.

PAGE LIMIT
  Full papers: at most EIGHT pages of main text. References and appendices do
  not count, and neither do the ethics and reproducibility statements — which
  is why those sit after the bibliography in main.tex. Do not move them back.

  The paper is currently at exactly 8 of 8. There is no slack: anything you add
  needs something cut. Check by finding which page the bibliography starts on.

COMPILER
  pdfLaTeX. Overleaf runs bibtex automatically; if references show as [?],
  recompile once more.

FIGURES
  All four are pre-built PDFs. Two are generated from committed result tables by
  paper/icbinb/make_figs.py and two are shared with the WMHS submission. If a
  number changes, regenerate them in the repository and re-bundle — editing a
  number here alone will silently desync it from the table it came from, which
  is the failure this paper is about.
EOF

if command -v zip >/dev/null 2>&1; then
  ( cd "$OUT" && zip -qr "../$OUT.zip" . )
  echo "wrote $OUT.zip  ($(find "$OUT" -type f | wc -l | tr -d ' ') files, $(du -h "$OUT.zip" | cut -f1))"
else
  echo "wrote $OUT/ (zip not installed — compress it yourself)"
fi
