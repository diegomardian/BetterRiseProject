#!/usr/bin/env bash
# Build the BMC manuscript and refuse to pass on a broken build.
#
# WHY THIS EXISTS. WMHS measured a page count on a build with 17 LaTeX errors;
# dropped content made the count a page short, and the wrong number nearly
# shipped. A build is only a measurement if the `errors` line is 0. This is the
# same gate: fail on any LaTeX error and on any undefined reference or citation.
#
# Under `article` this checks that the prose compiles and the bibliography
# resolves. THE BUILT PDF IS NOT SUBMISSION GEOMETRY until the journal's
# `bmcart` class is dropped in beside main.tex; no page count taken here means
# anything, and this script does not take one.
set -uo pipefail
cd "$(dirname "$0")"

doc=main
pdflatex -interaction=batchmode "$doc.tex" >/dev/null 2>&1
bibtex "$doc" >/dev/null 2>&1
pdflatex -interaction=batchmode "$doc.tex" >/dev/null 2>&1
pdflatex -interaction=batchmode "$doc.tex" >/dev/null 2>&1

if [ ! -f "$doc.pdf" ]; then
  echo "FAIL  $doc.pdf was not produced"
  exit 1
fi

errors=$(grep -c '^!' "$doc.log" || true)
undef=$(grep -c 'undefined' "$doc.log" || true)
echo "  $doc: errors $errors | undefined $undef"

status=0
if [ "$errors" -ne 0 ]; then
  echo "  FAIL  LaTeX errors; a page count off this build means nothing" >&2
  status=1
fi
if [ "$undef" -ne 0 ]; then
  echo "  FAIL  undefined references or citations" >&2
  status=1
fi

if [ "$status" -eq 0 ]; then
  echo "PASS. Built under 'article'; not submission geometry until bmcart is present."
fi
exit "$status"
