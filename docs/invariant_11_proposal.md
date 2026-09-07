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

## Deliberate non-implementation

This proposal changes the project's invariant contract and needs a `shared/...`
PR, a written reason, and two approvals under `CONTRIBUTING.md` §5. It must not
be added ad hoc to `CLAUDE.md`, `src/schema.py`, or a live analysis while a
result is under review. The implementation PR should include a failing fixture
where `GUCA2A` defines a mature-cell label and is then claimed as the endpoint,
and a passing fixture with an independent morphological or protein-based label.
