"""Target, availability and replay checks for the patient-t critique follow-up."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from patient_target_followup import fixed_targets
ROOT=Path(__file__).resolve().parent


def test_fixed_targets_distinguish_cell_and_patient_weighting():
    values=np.array([1.,1.,1.,9.]);patients=np.array(['a','a','a','b'])
    cell,patient,n=fixed_targets(values,patients,.5)
    assert n==2
    np.testing.assert_allclose([cell,patient],[-.6,-1.])
    # Replicating a patient's cells alters the cell target but not equal-patient weighting.
    more=np.array([1.,1.,1.,1.,1.,9.]);labels=np.array(['a']*5+['b'])
    c,p,_=fixed_targets(more,labels,.5);assert c!=cell and p==patient


def test_null_and_equal_sized_patient_targets_coincide():
    values=np.array([1.,3.,8.,12.]);labels=np.array(['a','a','b','b'])
    c,p,n=fixed_targets(values,labels,.5);assert c==p and n==2
    assert fixed_targets(values,labels,1.)[:2]==(0.,0.)


def test_replay_matches_original_counts_and_seed_aggregation():
    d=pd.read_csv(ROOT/'diagnostics/patient_target_followup.csv');s=pd.read_csv(ROOT/'diagnostics/patient_target_seed_rates.csv')
    keys=['cohort','pool','n_held_out','n_cells_mature','shift'];cols=['n_attempted','n_valid','n_abstained','n_rejected','n_cell_covered','n_patient_covered','n_draw_covered']
    np.testing.assert_array_equal(d.set_index(keys).sort_index()[cols],s.groupby(keys)[cols].sum().sort_index())
    old=pd.read_csv(ROOT/'diagnostics/incoming_interval_repair_pooled.csv');old=old[(old.candidate=='pseudobulk_patient_t')&old.n_cells_mature.isin([50,800])].set_index(keys).sort_index();a=d.set_index(keys).sort_index()
    for new,prior in [('n_valid','n_scored'),('n_abstained','n_abstained'),('n_rejected','n_excludes_zero'),('n_draw_covered','n_covered')]:np.testing.assert_array_equal(a[new],old[prior])
    assert len(d)==32 and d.n_attempted.sum()==64000


def test_targets_do_not_change_detection_or_hide_abstentions():
    d=pd.read_csv(ROOT/'diagnostics/patient_target_followup.csv');null=d[d['shift']==1]
    for flag in ['n_cell_covered','n_patient_covered','n_draw_covered']:np.testing.assert_array_equal(null[flag]+null.n_rejected,null.n_valid)
    np.testing.assert_allclose(d.unconditional_rejection,d.rejection*(1-d.abstention))
    assert (d.n_attempted==d.n_valid+d.n_abstained).all()
    for metric in ['cell_coverage','patient_coverage','draw_inclusion','rejection']:
        assert (d[metric+'_low']<=d[metric]).all() and (d[metric+'_high']>=d[metric]).all()


def test_reported_patient_tradeoff_and_five_patient_counterexample():
    d=pd.read_csv(ROOT/'diagnostics/patient_target_followup.csv');x=d[(d.n_cells_mature==50)&(d.n_held_out==2)&(d['shift']==.5)]
    np.testing.assert_allclose([x.cell_coverage.min(),x.cell_coverage.max()],[1257/1309,1944/1987])
    np.testing.assert_allclose([x.patient_coverage.min(),x.patient_coverage.max()],[1258/1309,1940/1987])
    k=d[(d.cohort=='kul3')&(d.n_cells_mature==50)&(d.n_held_out==5)&(d['shift']==1)]
    np.testing.assert_allclose(sorted(k.rejection),[.087,.181]);assert (k.rejection_low>.05).all()
    # Two-patient pool experiments are not analyses of all source-cohort patients.
    assert set(d.n_held_out)=={2,5}


def test_patient_result_provenance_and_versioned_parquet():
    p=json.loads((ROOT/'diagnostics/patient_target_provenance.json').read_text())
    for name,sha in p['input_sha256'].items():
        path=ROOT/name if name.endswith(('.py','.md')) else ROOT/'diagnostics'/name
        assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
    for name,sha in p['output_sha256'].items():assert hashlib.sha256((ROOT/'diagnostics'/name).read_bytes()).hexdigest()==sha
    for name in ['patient_target_followup','patient_target_seed_rates']:
        csv=pd.read_csv(ROOT/'diagnostics'/f'{name}.csv');parquet=pd.read_parquet(ROOT.parents[1]/'results/2026-09-16_102801f'/f'{name}.parquet')
        pd.testing.assert_frame_equal(csv,parquet,check_dtype=False,check_exact=False)
