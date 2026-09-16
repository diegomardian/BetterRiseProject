WMHS @ NeurIPS 2026 — Overleaf project
======================================
Bundled 2026-09-15. Source: papers/whms_bode on branch wmhs/build-fixes.

THE SUBMISSION
  main.tex — the only document here. 9 pages of main text (limit 9),
  15 pages including references and appendices, which do not count.

  "When a recovery curve fails to validate a decision rule"

BEFORE IT WILL COMPILE
  Add neurips_2026.sty to this project's root. It is the official workshop
  style and is deliberately not bundled. Do not modify it — the style file
  warns that tweaking it risks desk rejection.

COMPILER
  pdfLaTeX. Overleaf runs bibtex automatically; if references show as [?],
  recompile once more.

VERIFIED AT BUNDLE TIME
  Compiled from these files alone, in a clean directory:
    0 LaTeX errors
    0 undefined references or citations
    main text 9 of 9 pages
  Swept for mangled control sequences (a \textbf whose backslash-t became a
  tab once reached a reviewer's PDF as visible body text): clean.

PAGE COUNT — THE ONE THING TO RE-CHECK
  Measured against the official neurips_2024.sty, because the 2026 workshop
  style was not published when this was bundled. NeurIPS geometry has been
  stable across those years, so the count should hold — but it is NOT the file
  you will compile with, and the paper is at 9 of 9 with ZERO margin.

  Recompile with the real style and confirm before submitting. If it runs
  over, the cheapest cut is the residual-matrix table, which the appendix can
  carry.

  A page count from a build with errors is not a measurement — LaTeX drops the
  body of a table it cannot align, so a broken build is SHORTER than the real
  one. Check the log is clean first. That mistake has already happened on this
  paper once.

BEFORE YOU SUBMIT
  - Rebuild with the official style; confirm main text <= 9 pages.
  - Confirm the responsible-use statement is present. It is mandatory; a
    submission without one is desk-rejected.
  - Check anonymity: no author names, affiliations or repository URLs in the
    PDF or its metadata.
  - grep the log for "undefined" and confirm zero. Note the log writes
    "Reference `x' on page N undefined", so a naive grep for "' undefined"
    matches nothing — that exact mistake hid two broken references here.
