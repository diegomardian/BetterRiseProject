# Patient-t target check, fixed before running

This follow-up tests the critique's proposed resampling-unit explanation. It is
an exploratory re-scoring of the incoming simulated studies, not new patients.
The incoming patient-t procedure is kept unchanged. No method or cutoff is selected.

- Reuse ae8c109's generator, patient-t implementation and original five seed
  streams [20260831, 1, 2, 3, 4], 400 studies each, exactly.
- Four cohort/pool conditions; held-out patient counts 2 and 5; mature tumor
  counts 50 and 800; shifts 0.5 and 1: 32 settings, 64,000 studies.
- Compare three targets for the same returned patient-t interval: the incoming
  sampled-reference benchmark; the fixed cell-weighted source-pool mean contrast;
  and the fixed equally weighted mean contrast over source patients with mature
  cells. Compute both fixed targets before generating either arm.
- The equal-patient target conditions on the selected empirical pool too. It
  is not a population effect over new patients. Patients missing mature cells
  from the source pool are outside that conditional mean; record their count.
- Keep unavailable intervals explicit; score coverage and rejection conditional
  on availability and report denominators, abstention and unconditional detection.
- Reproduce incoming patient-t valid counts, rejection counts and sampled-reference
  inclusion counts on all 32 settings before interpreting the new coverage.
- Report Wilson intervals, source-patient counts, width and fixed-target distances.
  Rescoring cannot change rejection, detection or abstention on the same intervals.
- Code and design are committed before running; results carry that revision,
  input hashes and environments. Export aggregates only, as CSV and parquet.
