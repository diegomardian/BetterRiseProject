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

echo "built PDF:"
if [ -f main.pdf ]; then
  if hits=$(grep -aoiE "$PATTERNS" main.pdf 2>/dev/null | sort -u); then
    note "identifying strings inside main.pdf:"; echo "$hits" | sed 's/^/        /'
  else
    echo "  ok    no identifying strings in main.pdf"
  fi
  if command -v pdfinfo >/dev/null 2>&1; then
    if pdfinfo main.pdf | grep -qiE '^(Author|Keywords) +[^ ]'; then
      note "main.pdf carries Author/Keywords metadata"
    else
      echo "  ok    no Author/Keywords metadata"
    fi
  elif grep -aoE '/(Author|Keywords) ?\(([^)]+)\)' main.pdf >/dev/null 2>&1; then
    # No poppler: read the /Info dictionary out of the raw PDF instead. An
    # empty "/Author ()" is what pdflatex writes when no author is set.
    note "main.pdf carries Author/Keywords metadata:"
    grep -aoE '/(Author|Keywords) ?\(([^)]+)\)' main.pdf | sort -u | sed 's/^/        /'
  else
    echo "  ok    no Author/Keywords metadata (read from the raw PDF;"
    echo "        install poppler for pdfinfo if you want a second opinion)"
  fi
else
  note "main.pdf not built, so neither the PDF text nor its metadata was"
  echo "        checked — and metadata is where an author name leaks without"
  echo "        appearing on any page. Run ./build.sh first."
fi

echo
if [ "$fail" -eq 0 ]; then echo "PASS — safe to submit double-blind"; else echo "NOT SAFE TO SUBMIT"; fi
exit "$fail"
