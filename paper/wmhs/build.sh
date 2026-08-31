#!/usr/bin/env bash
# Build both submissions and enforce their page limits.
#
# main.tex       full paper, 9 pages of main text
# main_short.tex extended abstract, 4 pages of main text
#
# Both \input the same sections/, differing only in \iffull. The short build is
# a strict subset of the same prose, so a number cannot say one thing in one and
# something else in the other — which, given what this paper is about, is the
# one way it must not be wrong. What the short build drops, its appendix carries.
#
# References and appendices do not count toward either limit, so the check is
# "where does the bibliography start", not "how many pages is the PDF".
set -uo pipefail
cd "$(dirname "$0")"

if ! kpsewhich neurips_2026.sty >/dev/null 2>&1 && [ ! -f neurips_2026.sty ]; then
  echo "neurips_2026.sty not found. Download the official workshop style and put"
  echo "it beside main.tex — it is deliberately not vendored here." >&2
  exit 1
fi

build() {
  local doc=$1
  pdflatex -interaction=batchmode "$doc.tex" >/dev/null 2>&1
  bibtex "$doc" >/dev/null 2>&1
  pdflatex -interaction=batchmode "$doc.tex" >/dev/null 2>&1
  pdflatex -interaction=batchmode "$doc.tex" >/dev/null 2>&1
  [ -f "$doc.pdf" ]
}

fail=0
for doc in main main_short; do
  echo "building $doc …"
  if ! build "$doc"; then echo "  FAIL  $doc.pdf was not produced"; fail=1; continue; fi
  errors=$(grep -c '^!' "$doc.log")
  undef=$(grep -c 'undefined' "$doc.log")
  over=$(grep -c 'Overfull' "$doc.log")
  echo "  errors $errors | undefined $undef | overfull $over"
  [ "$errors" -eq 0 ] || fail=1
  [ "$undef"  -eq 0 ] || { echo "  FAIL  undefined references or citations"; fail=1; }
done

python3 - <<'PY' || fail=1
import sys
try:
    from pypdf import PdfReader
except ImportError:
    print("  --    pypdf not installed; page limits NOT checked"); sys.exit(0)
bad = 0
for doc, limit in (("main", 9), ("main_short", 4)):
    pages = [(p.extract_text() or "") for p in PdfReader(f"{doc}.pdf").pages]
    ref = next((i for i, t in enumerate(pages) if "References" in t), len(pages))
    spill = pages[ref].find("References") if ref < len(pages) else 0
    over = ref > limit or (ref == limit and spill > 60)
    print(f"  {doc:11s} main text {ref} page(s), limit {limit}"
          f"{'  FAIL — over' if over else '  ok'}")
    bad |= over
sys.exit(1 if bad else 0)
PY

echo
if [ "$fail" -eq 0 ]; then
  echo "PASS. Next: ./check_anonymity.sh (it must exit 0 before you submit)."
else
  echo "BUILD NOT SUBMITTABLE"
fi
exit "$fail"
