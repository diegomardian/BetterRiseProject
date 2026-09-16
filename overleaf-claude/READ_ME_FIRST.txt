WMHS @ NeurIPS 2026 — Overleaf project (alternate revision)
===========================================================
Bundled 2026-09-15 from papers/wmhs_claude on branch wmhs/claude-revision.

THIS IS THE ALTERNATE VERSION, FOR COMPARISON
  There are two revisions of this paper in the repository:

    papers/whms_bode   the main line, on branch wmhs/build-fixes
    papers/wmhs_claude THIS ONE, on branch wmhs/claude-revision

  They share a common ancestor and most of their text. Compare and pick one;
  do not submit both.

WHAT THIS VERSION HAS THAT THE MAIN LINE DOES NOT
  1. The trial-simulator sections restored. Commit ccb6264 deleted
     sections/blind.tex and dropped sections/refdesign.tex from main.tex while
     condensing the appendix, which removed the equality identity, the
     reference-and-design result, the learned-generator experiment and the
     censoring/RMST material from the compiled paper. All are back.
  2. The interval repair study. Five candidate intervals plus a balanced-arm
     design change, none of which reaches the 5% null target at 50 cells, and
     the patient-level candidates which never pass anywhere. This is the direct
     answer to the reviewer's "diagnoses without prescribing".
  3. The assignment-design result (Appendix "Assignment design and the
     equality"): under complete 1:1 permuted blocks even the unadjusted
     contrast reproduces the reference, while Bernoulli assignment at the same
     probability does not.
  4. A correction. An earlier draft reported 9.5% null rejection at 800 cells;
     at 2,000 replicates rather than 200 that cell measures 6.85% +/- 0.56.
     The original figure was Monte Carlo noise. The 50-cell rates, which the
     argument rests on, reproduce.

  The single-cell detail that was in the main text (fixed-pool evaluation, the
  cohort gap, the fixed-reference control) moved to the appendix to make room.
  Appendices do not count toward the page limit.

THE SUBMISSION
  main.tex — 9 pages of main text (limit 9), 18 pages including references
  and appendices.

BEFORE IT WILL COMPILE
  Add neurips_2026.sty to this project's root. It is the official workshop
  style and is deliberately not bundled.

VERIFIED AT BUNDLE TIME
  Compiled from these files alone in a clean directory: 0 LaTeX errors,
  0 undefined references, 0 missing files, main text 9 of 9.
  Swept for mangled control sequences: clean.

PAGE COUNT — RE-CHECK
  Measured against the official neurips_2024.sty, because the 2026 workshop
  style was not published when this was bundled. Geometry has been stable
  across those years, but this is NOT the file you will compile with, and
  there is ZERO margin at 9 of 9.

  A page count from a build with errors is not a measurement: LaTeX drops the
  body of a table it cannot align, so a broken build comes out SHORTER. Check
  the log is clean first. And note the log writes "Reference `x' on page N
  undefined" — a naive grep for "' undefined" matches nothing.
