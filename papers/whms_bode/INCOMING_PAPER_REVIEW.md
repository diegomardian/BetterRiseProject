# Review of the incoming paper and adopted changes

Reviewed after pulling `wmhs/build-fixes` from `5d0baef` to `a336594` on
2026-09-15. The incoming Overleaf paper is
`paper/wmhs/overleaf-submission/main.tex`; its README identifies
`papers/whms_bode` as its source. Its substantive additions over our first
critique revision are the external-control demonstration in validation and
appendix sections. I also reviewed the changes in `paper/wmhs/sections/full`.
The later local fixed-pool/false-positive revision was backed up and restored
before integrating the incoming work. That revision remains the foundation.

## Assessment

**The external-control example is worth adding, with substantial qualifications.**
It gives the paper a concrete consequence of reference-population mismatch.
The same two estimators change ranking when the validation reference changes.
It also improves the proposed remedy: potential-outcome information alone is
insufficient if it is averaged over the wrong population. These findings
connect the single-cell target analysis to a trial-like simulation without
claiming that the single-cell experiments validate a clinical world model.

The new main-text table uses bias and RMSE instead of maximum residual and
recovery ratio. These directly answer whether the estimator recovers the
intended effect, and retain the revealing wrong-population control. All
entries come from the pinned results, not transcribed narrative summaries.
The abstract, motivation, conclusion and limitations now acknowledge this
second experiment. The long incoming appendix discussion is replaced with a
compact reproducible specification.

## Findings that required correction or qualification

| Incoming claim | Independent assessment and treatment |
|---|---|
| Six seeds | The module and sidecar use **one seed, 20260915**, and one fixed 50,000-record registry. Corrected. Six refers to cohort sizes, not seed streams. |
| Zero exclusions / 48,000 primary draws | There are **1,200 study draws**, evaluated by 8 estimators against 5 references: 48,000 comparisons. **47,969 are finite and 31 excluded**, all from replicate 53 at n=50. Corrected. The n=1,600 main table has 200 finite observations in every cell. |
| Clinical progression-free-survival benefit in months | The generator uses uncensored Normal outcomes, which can be negative; no clinical registry was used. Reported as abstract outcome units and a synthetic continuous-outcome illustration, not a realistic survival model. |
| Pooled mix is the registry's mix | A 3:1 control ratio makes the pooled file control-dominated, not identical to the registry. It combines enrolled and remaining-registry controls; finite-pool sampling also changes the realized mix. Corrected. |
| Truth is unavailable / theta does not exist | The finite-cohort effect is well defined and available to the simulator through both potential outcomes. It is unavailable to estimators that observe one outcome per record. Population specification is the issue. Corrected. |
| The correct estimator is unbiased and error falls as n^-1/2 | The run shows small measured bias and falling RMSE, not a proof of exact unbiasedness or a convergence law. Reported endpoint RMSEs (0.608 to 0.105) and the measured table values instead. |
| A 1.89-fold overstatement causes a wrong clinical decision | At n=1,600, the median pooled estimate is 1.890 and the median enrolled truth 1.003, on opposite sides of the illustrative 1.5 threshold. The saved shift analysis compares medians; it does **not** estimate a decision-error rate. Explicitly qualified. |
| Every routine safeguard fails / no protocol tabulates target balance | The experiment evaluates specific measured-stratum balance and reference checks. It cannot establish what every protocol does. Retained the measured between-arm versus enrolled-distribution contrast; removed universal practice claims. |
| Breaking equality is necessary and insufficient | Correct target-population specification is essential, but a useful check need not disagree with a correct estimator on every dataset. The zero-shift control itself shows equality with little bias. Removed the general necessity claim. |
| Defective versus correct methods | Pooled standardization is an estimator of a different population effect. It is wrong for the stated enrolled target, not intrinsically an invalid estimator. Made that distinction explicit. |
| Nine falsifiers passing establishes realism | It establishes consistency with those prespecified checks. It does not show clinical realism, prevalence, or a real deployment error. The same population mistake in estimator and reference is an explicit assumption. |

The zero-shift result in the saved table is approximately **-0.0017**, not
+0.0017: the incoming narrative quoted its magnitude as if signed. The paper
now gives the signed rounded bias, -0.002. The incoming 1.82–1.93 estimate
range also mixes estimator summaries; the revision uses the unambiguous
n=1,600 result for the two methods actually tabulated.

## Material not carried over

- Stronger claims about observed clinical harm, widespread practice or a
  generally necessary non-equality rule: not supported by the experiment.
- Regulatory-use assertions and new regulatory citations: they are unnecessary
  for the narrow experimental claim and do not validate this synthetic design.
- The longer algebra-first framing and residual-rate discussion in the other
  full-paper tree: our fixed-pool, null-rejection and candidate-reversal results
  provide stronger evidence. The old count of fourteen tolerance switches is
  still inconsistent with the independently checked count of thirteen.
- The incoming bundle's obsolete style caveat: this revision already uses the
  official unmodified 2026 style, which is included in the portable source ZIP.

## Verification

- Independently reran the full primary experiment: all **240 summary rows**
  agree with the saved table (1,200 study draws; numerical tolerance 1e-11
  absolute and 1e-9 relative).
- Independently reran the 50-draw balance experiment: all **9 summary rows** agree.
- Replayed n=50 to identify the excluded study and reconcile all 31 comparisons.
- **24 paper claim tests and 35 external-control harness tests pass.** New checks
  cover target-dependent rankings, true denominators, zero-shift control,
  balance, provenance checksums and exported aggregate results.
- Four new result tables are version-pinned alongside the previous 38 inputs;
  primary, shift and balance summaries are included in the anonymous source
  bundle. Exact-source provenance remains outside that bundle.
- The final PDF has **8 main-text pages, 2 reference pages and 3 appendix
  pages** (13 total), using the unchanged official style. All 13 pages were
  rendered and visually inspected; anonymity and reference checks pass.
- The source ZIP compiles in a clean directory with no errors or unresolved
  references, and all 13 extracted page texts match the delivered PDF.

The simulator and its committed results were not changed. This is a critical
integration of the useful results into our paper, not an endorsement of all
claims in the incoming version.
