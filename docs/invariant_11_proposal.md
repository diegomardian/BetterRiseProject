# Proposed invariant 11: no circular programme claim

## Proposed assertion

No transcript-derived label may define a population and then be used to claim
the state of its own defining programme.

This is the testable half of the ML policy in `CONTRIBUTING.md`: ML may provide
measurements, but the resulting biological claim still needs an independent
population definition or an explicit abstention rule.

## Guard design

Before a result writer accepts an analysis specification, it should require two
declared sets:

- `population_definition_genes`: genes or gene programmes used to assign the
  cells, regions, or crypt positions being compared;
- `claim_genes`: genes used as the biological endpoint.

The guard fails when their intersection is non-empty, unless the specification
declares an independently measured endpoint and names the independent assay.
The failure message must name the overlapping genes and say that the analysis
is circular, rather than silently proceeding.

## RESULT — implemented 2026-09-07

**Approved by W2 and implemented as `CLAUDE.md` invariant 11**, with the three
specification refinements the review asked for:

1. **The declaration is mandatory and omission is not a pass.** An analysis with
   gene-free labels declares `genes=()`; leaving the declaration out raises the
   same way a circular one does. This is invariant 1's rule in another place —
   unstated is not none.
2. **The exception is structural, not a name.** An overlap needs four things a
   name cannot fake: the declaration must name the assay that actually produced
   the endpoint, that assay must differ from the labels' assay, the modality
   must differ, and the labels must not be `derived_from` the endpoint's own
   assay. The stated reason is recorded and never evaluated. The fourth
   condition catches the case the first three miss — segment on MxIF CDX2, call
   the regions mature by morphology, then report CDX2 protein in them.
3. **It lives at analysis-specification validation**, in
   `src/common/label_provenance.py`, called by each job's
   `validate_specification()`. Not in `write_versioned_table`: a writer sees a
   table of numbers, by which point the population is chosen and the columns no
   longer say what defined it, so a guard there could only check that two
   fields were filled in.

Both fixtures the proposal asked for exist. The failing one — GUCA2A defines a
mature label and is then claimed — and the passing one, CDX2 protein by MxIF
against a CDX2-transcript label, are in `tests/test_checks_can_fail.py` and
`tests/test_label_provenance.py`. The two live reference jobs (Becker, Crowell
multisection) now carry declarations, and a test asserts the guard can read
both, so the invariant cannot decay into documentation.

**Review round two, 2026-09-07.** Two fixes, both about the gap between a
guard that works and an invariant a reader can see:

- **It fired too late.** `validate_specification()` was called inside
  `verdict()`, after the sections were opened and the population-level
  quantities computed — so it could only object to work already done. Both jobs
  now call it immediately after `parse_args`, before a path is even checked.
  `verdict()` keeps its own call for API callers. The forcing test points each
  job at a nonexistent path with a circular declaration: the circularity must be
  what stops it, which is only true if the guard runs first.
- **The sidecars did not carry it.** `provenance_meta` existed and neither job
  used it, so a reader holding a parquet could not see what defined the
  population. Both now fold it into `extra_meta`, and the test asserts the
  wiring and round-trips the JSON.

**Known limit.** The guard binds where a job calls it. Nothing yet forces a
*new* job to call it, because there is no shared analysis-specification object
in this repository to hang it on — specifications are per-job constants. The
honest statement is that invariant 11 is enforced in the two analyses that have
labels today and is a convention for the next one.

## Deliberate non-implementation — superseded, kept for the record

This proposal changes the project's invariant contract and needs a `shared/...`
PR, a written reason, and two approvals under `CONTRIBUTING.md` §5. It must not
be added ad hoc to `CLAUDE.md`, `src/schema.py`, or a live analysis while a
result is under review. The implementation PR should include a failing fixture
where `GUCA2A` defines a mature-cell label and is then claimed as the endpoint,
and a passing fixture with an independent morphological or protein-based label.
