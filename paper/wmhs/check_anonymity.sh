#!/usr/bin/env bash
# Double-blind guard. Run before every submission build.
#
# The paper's own rule is that a check unable to fail is worse than no check, so
# this one is written to fail: it greps the sources AND the built PDF for the
# things that actually de-anonymise a submission, and checks the PDF's own
# metadata, which is where an author name leaks without appearing on any page.
set -uo pipefail
cd "$(dirname "$0")"
fail=0
note() { echo "  FAIL  $*"; fail=1; }

# Identifying strings. Extend this list, never shorten it.
PATTERNS='diegomardian|BetterRiseProject|bodebosell|[Mm]eowsers|github\.com/[A-Za-z0-9_-]+|/Users/|/home/[a-z]'

echo "sources:"
if hits=$(grep -rniE "$PATTERNS" main.tex sections/ refs.bib 2>/dev/null); then
  note "identifying strings in the LaTeX sources:"; echo "$hits" | sed 's/^/        /'
else
  echo "  ok    no identifying strings in main.tex, sections/, refs.bib"
fi

echo "built PDFs:"
# Both builds ship, so both get checked. A missing PDF is a FAIL, not a note:
# this script's whole job is to be the thing that cannot pass without looking,
# and "not built" was passing.
for pdf in main.pdf main_short.pdf; do
  if [ ! -f "$pdf" ]; then
    note "$pdf not built - run ./build.sh, then re-run this"
    continue
  fi
  if hits=$(grep -aoiE "$PATTERNS" "$pdf" 2>/dev/null | sort -u); then
    note "identifying strings inside $pdf:"; echo "$hits" | sed 's/^/        /'
  else
    echo "  ok    no identifying strings in $pdf"
  fi

  # Metadata is where an author name leaks without appearing on any page, so a
  # missing reader is a FAIL too - silently skipping it is the same defect.
  if command -v pdfinfo >/dev/null 2>&1; then
    if pdfinfo "$pdf" | grep -qiE '^(Author|Keywords) +[^ ]'; then
      note "$pdf carries Author/Keywords metadata"
    else
      echo "  ok    no Author/Keywords metadata in $pdf"
    fi
  elif python3 -c 'import pypdf' >/dev/null 2>&1; then
    meta=$(python3 -c '
import sys
from pypdf import PdfReader
info = PdfReader(sys.argv[1]).metadata or {}
bad = {k: v for k, v in info.items()
       if k in ("/Author", "/Keywords") and str(v).strip()}
print("; ".join(f"{k}={v}" for k, v in bad.items()))
' "$pdf")
    if [ -n "$meta" ]; then
      note "$pdf carries Author/Keywords metadata: $meta"
    else
      echo "  ok    no Author/Keywords metadata in $pdf"
    fi
  else
    note "cannot read $pdf metadata: install poppler (pdfinfo) or pypdf"
  fi
done

echo
if [ "$fail" -eq 0 ]; then echo "PASS — safe to submit double-blind"; else echo "NOT SAFE TO SUBMIT"; fi
exit "$fail"
