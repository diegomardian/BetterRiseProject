#!/usr/bin/env bash
# Double-blind guard. Run before every submission build.
#
# The paper's own rule is that a check unable to fail is worse than no check, so
# this one is written to fail: it greps the sources AND the built PDF for the
# things that actually de-anonymise a submission, and it fails on the unfilled
# artifact placeholder too, so §5 cannot promise a code release that points
# nowhere.
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

echo "artifact link:"
if grep -q "ANONYMISED-ARTIFACT-URL" main.tex; then
  note "the artifact URL is still the placeholder — §5 promises a code release."
  echo "        Create an anonymising mirror (e.g. anonymous.4open.science) and"
  echo "        replace \\artifacturl in main.tex. Never the real repository."
else
  echo "  ok    artifact URL filled in"
  if grep -E '\\newcommand\{\\artifacturl\}' main.tex | grep -qiE 'github\.com|gitlab'; then
    note "the artifact URL points at a real forge — that de-anonymises the submission."
  fi
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
