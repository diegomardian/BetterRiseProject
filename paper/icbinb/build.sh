#!/usr/bin/env bash
# Build the ICBINB submission and enforce its page limit.
#
# CFP, verified 2026-08-31: full papers up to EIGHT pages excluding references
# and appendices; tiny papers four. Ethics and reproducibility statements do not
# count toward the limit either. So the check is "where does the bibliography
# start", not "how many pages is the PDF".
#
# Deadline: 2 September 2026, 11:59pm AoE.
set -uo pipefail
cd "$(dirname "$0")"
LIMIT=8

if ! kpsewhich neurips_2026.sty >/dev/null 2>&1 && [ ! -f neurips_2026.sty ]; then
  echo "neurips_2026.sty not found — download the official style and put it here." >&2
  exit 1
fi

pdflatex -interaction=batchmode main.tex >/dev/null 2>&1
bibtex main >/dev/null 2>&1
pdflatex -interaction=batchmode main.tex >/dev/null 2>&1
pdflatex -interaction=batchmode main.tex >/dev/null 2>&1

fail=0
errors=$(grep -c '^!' main.log); undef=$(grep -c 'undefined' main.log)
over=$(grep -c 'Overfull' main.log)
echo "  errors $errors | undefined $undef | overfull $over"
[ "$errors" -eq 0 ] || fail=1
[ "$undef" -eq 0 ] || { echo "  FAIL  undefined references or citations"; fail=1; }

python3 - "$LIMIT" <<'PY' || fail=1
import sys
try:
    from pypdf import PdfReader
except ImportError:
    print("  --    pypdf not installed; page limit NOT checked"); sys.exit(0)
limit = int(sys.argv[1])
pages = [(p.extract_text() or "") for p in PdfReader("main.pdf").pages]
ref = next((i for i, t in enumerate(pages) if "References" in t), len(pages))
spill = pages[ref].find("References") if ref < len(pages) else 0
over = ref > limit or (ref == limit and spill > 60)
print(f"  main text {ref} page(s), limit {limit}{'  FAIL — over' if over else '  ok'}")
sys.exit(1 if over else 0)
PY

echo
[ "$fail" -eq 0 ] && echo "PASS" || echo "NOT SUBMITTABLE"
exit "$fail"
