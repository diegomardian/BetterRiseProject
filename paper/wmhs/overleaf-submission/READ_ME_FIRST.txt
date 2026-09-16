WMHS @ NeurIPS 2026 — Overleaf project
======================================
Bundled 2026-09-15 from branch wmhs/build-fixes, commit b503971.

WHICH DOCUMENT IS THE SUBMISSION
  main.tex        <-- THIS ONE. Full paper, 9 pages of main text (limit 9).
  main_short.tex  Extended abstract, 4 pages of main text (limit 4).

  Both build and both fit. main.tex is the recommendation. Select it under
  Menu > Settings > Main document.

  main.tex reads sections/full/. main_short.tex reads sections/. Figures and
  refs.bib are shared. Both trees are bundled; do not delete either.

BEFORE IT WILL COMPILE
  Add neurips_2026.sty to this project's root. It is the official workshop
  style and is deliberately not bundled. Do not modify it — the style file
  warns that tweaking it risks desk rejection.

COMPILER
  pdfLaTeX. Overleaf runs bibtex automatically; if references show as [?],
  recompile once more.

PAGE COUNTS — READ THIS
  The 9-of-9 and 4-of-4 figures above were measured against the official
  neurips_2024.sty behind a compatibility shim, because the 2026 workshop
  style had not been published when this was bundled. NeurIPS page geometry
  has been stable across those years, so the counts should hold — but they
  are NOT measured against the file you will compile with.

  The full paper is at 9 of 9 with ZERO margin. Recompile with the real style
  and check the page count before submitting. If it runs over, the cheapest
  cut is the residual-matrix table in section 3, which the appendix can carry.

  A page count from a build with errors is not a measurement. Check the log is
  clean first — that mistake has already happened once on this paper.

BEFORE YOU SUBMIT
  - Rebuild with the official style; confirm main text <= 9 pages.
  - Confirm the responsible-use statement is present (section 5). It is
    mandatory; a submission without one is desk-rejected.
  - Check anonymity: no author names, no affiliations, no repo URLs in the
    PDF or its metadata.
  - grep the log for "undefined" and confirm zero.

WHAT WAS VERIFIED AT BUNDLE TIME
  Compiled standalone from these files alone: 0 errors, 0 undefined
  references, 0 missing files, for both documents. 55 tests passing in the
  source repository, which re-derive the quoted numbers from versioned result
  tables.
