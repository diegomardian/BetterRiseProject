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

  IT FITS, AND ONLY JUST. Compiled against the NeurIPS geometry, the main text
  fills page 8 to the bottom and the References heading sits at the TOP of
  page 9. There is no margin: anything you add needs something cut first.

  Check where the MAIN TEXT ENDS, not where the bibliography starts. Those are
  different questions and it is easy to measure the wrong one -- a draft of
  this paper spilled ~81 words onto page 9 while a checker that looked for the
  References heading still reported "8 pages", because the heading had simply
  moved down the page rather than onto a new one. build.sh measures the spill;
  trust it over eyeballing.

  figures/fig_tiers.pdf is bundled but NOT included by main.tex. It was cut to
  make the limit; section 4 states every number it showed. If the real style
  turns out to leave room, re-adding it is one \includegraphics line.

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
