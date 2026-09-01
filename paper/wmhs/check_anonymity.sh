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
  fi
else
  echo "  --    main.pdf not built; run the build and re-run this"
fi

echo
if [ "$fail" -eq 0 ]; then echo "PASS — safe to submit double-blind"; else echo "NOT SAFE TO SUBMIT"; fi
exit "$fail"
