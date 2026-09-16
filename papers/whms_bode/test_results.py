"""Claim checks for this revision, using immutable saved result tables.

These tests do not rerun the simulations. They verify denominators, selections,
aggregation, and the new manuscript's quantitative assertions.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from _tables import pinned
from make_tables import render

PAPER=Path(__file__).resolve().parent

def table(name):return pd.read_parquet(pinned(name))
def prose(name):
    from check_submission import source_graph
    return '\n'.join(p.read_text() for p in source_graph(PAPER/'sections'/name))

def test_expanded_learned_design_denominators():
    d=table('trial_learned_generator_primary')
    draw=d[d.reference=='T_draw']
    assert len(draw)==7*10*6==420
    assert set(draw.n_draws)=={178,180}
    assert draw.n_draws.sum()==75580
    assert draw.n_valid_pairs.sum()==74726
    excluded=(draw.n_draws-draw.n_valid_pairs).sum()
    assert excluded==854
    assert draw.n_nonfinite_reference.sum()==410
    assert draw.n_nonfinite_estimate.sum()==731
    assert 410+731-excluded==287
    assert 7*6*6*30==7560
    assert draw.groupby(['generator','n_patients']).n_draws.first().sum()==7558
    assert (draw.n_valid_pairs>0).all()
    a=prose('refdesign.tex')
    for literal in ['75{,}600','75{,}580','74{,}726','854','410','731','287','7{,}558']:
        assert literal in a


def test_learned_table_includes_every_estimator_and_uses_correct_groups():
    assert prose('learned_table.tex')==render()
    d=table('trial_learned_generator_primary')
    d=d[d.reference=='T_draw']
    mixture=d[d.generator.str.startswith('gmm')]
    ridge=mixture[mixture.estimator=='gcomp-ridge-outcome-a1e-5']
    assert f'{ridge.max_residual.max():.1e}'=='7.9e-05'
    mlp=mixture[mixture.estimator=='gcomp-mlp-outcome']
    assert round(mlp.max_residual.max(),1)==3.3
    assert 'G-comp., MLP outcome' in prose('learned_table.tex')


def test_learned_identity_is_nonvacuous():
    d=table('trial_learned_generator_primary')
    d=d[(d.reference=='T_draw')&(d.estimator=='gcomp-saturated')]
    assert len(d)==42 and d.is_learned.sum()==30
    assert (d.max_residual==0).all() and (d.n_valid_pairs>0).all()


def test_tolerance_switch_uses_actual_comparisons_not_a_theta_approximation():
    d=table('trial_learned_generator_primary')
    d=d[d.reference=='T_draw']
    band=d[(~d.allclose_residual_vs_zero)&d.allclose_estimate_vs_reference]
    assert len(band)==41 and band.estimator.str.contains('ridge').all()
    a=d[d.n_patients==100].set_index(['generator','estimator'])
    b=d[d.n_patients==5000].set_index(['generator','estimator'])
    switched=(~a.allclose_residual_vs_zero)&(~a.allclose_estimate_vs_reference)&(~b.allclose_residual_vs_zero)&b.allclose_estimate_vs_reference
    assert switched.sum()==13
    assert 'Thirteen generator--estimator combinations' in prose('refdesign.tex')


def test_model_reference_uncertainty_is_only_monte_carlo_for_mixtures():
    d=table('trial_learned_generator_primary');d=d[d.reference=='T_model']
    mc=d[d.T_model_is_monte_carlo]
    assert mc.generator.str.startswith('gmm').all()
    assert round(mc.T_model_mc_se.min(),3)==.019
    assert round(mc.T_model_mc_se.max(),3)==.027
    assert (d[~d.T_model_is_monte_carlo].T_model_mc_se==0).all()
    assert 'errors are 0.019--0.027' in prose('refdesign.tex')


def test_common_count_figure_deduplicates_identical_grid_views():
    d=table('controlled_grid_rates_r200_b200');d=d[d.binning=='per_count']
    keys=['cohort','pool','seed','n_cells_mature']
    for c in ['coverage','discrimination','n_estimated','n_attempted']:
        assert d.groupby(keys)[c].nunique(dropna=False).max()==1
    d=d.drop_duplicates(keys)
    assert d.groupby(['cohort','pool','n_cells_mature']).seed.nunique().eq(8).all()
    assert (d[d.n_cells_mature==0].n_estimated==0).all()


def test_direct_performance_at_committed_count():
    d=table('controlled_grid_rates_r200_b200')
    d=d[(d.binning=='per_count')&(d.grid=='dense')&(d.n_cells_mature==50)&(d.pool=='reference')]
    s=d[d.cohort=='smc'];k=d[d.cohort=='kul3']
    assert len(s)==len(k)==8
    assert (s.discrimination.min(),s.discrimination.max())==(.755,.81)
    assert (k.discrimination.min(),k.discrimination.max())==(.45,.55)
    assert ((s.coverage>=.9)&(s.discrimination>=.8)).sum()==2
    assert not ((k.coverage>=.9)&(k.discrimination>=.8)).any()
    # Historical seed rates remain checked even after the prose is condensed.



def test_controlled_candidates_and_seed_counts():
    d=table('controlled_grid_crossings_r200_b200')
    d=d[d.criterion=='coverage_and_discrimination']
    assert d[d.pool=='pooled'].candidate.isna().all()
    d=d[d.pool=='reference']
    def candidates(cohort,grid,binning):
        return d[(d.cohort==cohort)&(d.grid==grid)&(d.binning==binning)].candidate.value_counts().to_dict()
    assert candidates('smc','committed','fixed_bins')=={100.:8}
    assert candidates('smc','extended','fixed_bins')=={70.:8}
    assert candidates('smc','dense','fixed_bins')=={70.:8}
    assert candidates('smc','extended','adaptive_bins')=={65.:7,90.:1}
    assert candidates('smc','extended','per_count')=={60.:5,50.:2,80.:1}
    assert set(d[(d.cohort=='smc')&(d.grid=='dense')&(d.binning=='per_count')].candidate)=={45.,50.,60.,80.}
    for grid in d.grid.unique():
        for binning in d.binning.unique():assert candidates('kul3',grid,binning)=={400.:6,800.:2}
    rates=table('controlled_grid_rates_r200_b200')
    assert round(rates[['coverage_mc_se','discrimination_mc_se']].max().max(),3)==.035


def test_patient_influence_and_inner_budget():
    d=table('calibration_inner_budget_crossings_b200-1000-5000')
    d=d[d.criterion=='coverage_and_discrimination']
    assert d[d.cohort=='smc'].sort_values('inner_bootstrap').candidate.tolist()==[70.,70.,70.]
    assert d[d.cohort=='kul3'].sort_values('inner_bootstrap').candidate.tolist()==[400.,300.,400.]
    d=table('calibration_lopo_summary_b200-1000-5000').set_index(['cohort','criterion'])
    for key,value in [(('smc','coverage_and_discrimination'),0),(('kul3','coverage_and_discrimination'),3),(('smc','coverage_only'),5),(('kul3','coverage_only'),0)]:
        assert d.loc[key,'n_changing_conclusion']==value
    assert d.loc[('kul3','coverage_and_discrimination'),'changing_patients']=='KUL01,KUL30,KUL31'
    d=table('calibration_lopo_influence_b200-1000-5000')
    r=d[(d.cohort=='kul3')&(d.criterion=='coverage_and_discrimination')&(d.omitted_patient=='KUL31')].iloc[0]
    assert r.omitted_candidate==200 and r.omitted_status=='lower_bound_unobserved'


def test_clean_control_coverage_examples_match_saved_results():
    d=table('residual_performance_clean_control')
    assert len(d)==8 and (d.n_valid==2000).all()
    text=prose('refdesign.tex')
    for r in d[d.estimator.isin(['empirical-mean-calibrated','same-point-narrow-interval'])].itertuples():
        assert r.max_residual_vs_reference==0
        assert f'{100*r.interval_coverage:.2f}\\%' in text
    assert d.coverage_mc_se.max()*100<=1.08


def test_recovery_counts_and_raw_zero_cell_caveat():
    d=table('calibration_gap_recovery')
    assert d.n_replicates.sum()==119600 and d.n_ratio_undefined.sum()==29900
    assert (d.max_abs_residual_vs_realised==0).all()
    assert d.seed.nunique()==13
    assert 'raw diagnostic values arise from' in prose('appendix.tex')
    assert 'draw-dependent' in prose('calibration.tex')
    assert r'i_{\mathrm{req}}(D_N;s)=f_N(s-1)\bar Y_N' in prose('calibration.tex')


def test_survival_and_cox_claims():
    d=table('trial_survival')
    e=d[d.estimator=='exponential-mle-standardised'];km=d[d.estimator=='km-rmst-standardised']
    assert len(e)==len(km)==28 and (e.max_residual_vs_observed==0).all()
    assert e.groupby('censoring_target').max_residual_vs_latent.max().round(3).tolist()==[0.,.473,1.089,1.586]
    assert f'{km.max_residual_vs_latent.max():.1e}'=='2.8e-15'
    h=table('trial_survival_headline');assert 1200-h.n_replicates.min()==463
    d=table('trial_survival_lifelines_check')
    assert len(d)==8
    # Column names are checked rather than silently accepting an unrelated max.
    diff=[c for c in d if 'diff' in c and ('abs' in c or 'absolute' in c)]
    assert len(diff)==1
    assert f'{d[diff[0]].max():.1e}'=='1.4e-07'


def test_external_audit_scope_and_count():
    d=table('prevalence_audit_verdicts')
    assert len(d)==7 and (d.verdict=='NO').all()
    assert (d.code_executed=='yes').sum()==1
    assert d.max_residual.notna().sum()==1
    assert round(d.max_residual.max()*1e14,1)==4.0
    assert len(table('prevalence_audit_screening'))==19
    assert 'supports no prevalence estimate.' in prose('appendix.tex')
    assert 'already the norm' not in '\n'.join(x.read_text() for x in (PAPER/'sections').glob('*.tex'))


def test_low_count_diagnostic_reproduces_saved_alternative():
    d=pd.read_csv(PAPER/'diagnostics/low_count_diagnostic.csv')
    assert len(d)==32 and d.n_valid.sum()==6400
    assert (d.n_valid==d.n_attempted).all()
    original=table('controlled_grid_rates_r200_b200')
    original=original[(original.binning=='per_count')&(original.seed==20260831)].drop_duplicates(['cohort','pool','n_cells_mature'])
    paired=d[d['shift']==.5].merge(original,on=['cohort','pool','n_cells_mature'],suffixes=('_diagnostic','_original'))
    assert len(paired)==16
    for a,b in [('coverage_diagnostic','coverage_original'),('exclusion_rate','discrimination'),('median_ci_width_diagnostic','median_ci_width_original')]:
        np.testing.assert_allclose(paired[a],paired[b],rtol=1e-12,atol=1e-12)


def test_null_rejection_and_sparse_sample_claims():
    from make_diagnostic_table import render as diagnostic_render
    d=pd.read_csv(PAPER/'diagnostics/low_count_diagnostic.csv')
    null=d[(d['shift']==1)&(d.n_cells_mature==5)]
    assert sorted(null.n_excludes_zero)==[75,123,126,134]
    assert (d.n_zero_width==0).all()
    sparse=d[(d['shift']==.5)&(d.n_cells_mature==5)&(d.pool=='pooled')].set_index('cohort')
    assert sparse.loc['smc','n_tumour_all_zero']==96
    assert sparse.loc['kul3','n_tumour_all_zero']==85
    assert (sparse.n_tumour_all_zero==sparse.n_excludes_zero_and_tumour_all_zero).all()
    assert prose('low_count_table.tex')==diagnostic_render()


def test_later_failures_restrict_first_crossing_interpretation():
    from make_diagnostic_table import reversals
    rows=reversals()
    assert rows==json.loads((PAPER/'diagnostics/candidate_reversals.json').read_text())
    failures=[r for r in rows if r['later_failures']]
    assert len(failures)==3
    assert all(r['cohort']=='smc' and r['pool']=='reference' for r in failures)
    assert sorted((r['first_candidate'],r['all_subsequent_evaluated_candidate']) for r in failures)==[(45.,60.),(50.,70.),(50.,70.)]


def test_documentary_claim_is_inspectable_and_matches_source():
    import hashlib
    for record in json.loads((PAPER/'evidence/case_provenance.json').read_text()):
        source=PAPER.parents[1]/record['source']
        assert hashlib.sha256(source.read_bytes()).hexdigest()==record['sha256']
        assert record['excerpt'] in source.read_text()
    review=json.loads((PAPER/'evidence/review_record.json').read_text())
    assert review['artifact_A']['replicates_per_setting']==6
    original=json.loads((PAPER/'evidence/case_provenance.json').read_text())[0]['excerpt']
    assert review['artifact_A']['original_statement'] in ' '.join(original.replace('**','').split())
    assert review['artifact_B']['minimum_target_coverage']==.9
    assert 'record A' in prose('calibration.tex')


def test_fixed_pool_coverage_and_paired_accounting():
    from make_fixed_pool_tables import fixed_table
    d=pd.read_csv(PAPER/'diagnostics/low_count_diagnostic.csv')
    assert prose('fixed_pool_table.tex')==fixed_table(d)
    assert (d.n_fixed_pool_covered-d.n_draw_benchmark_covered==d.n_fixed_only-d.n_draw_only).all()
    np.testing.assert_allclose(d.fixed_pool_coverage-d.coverage,d.paired_coverage_difference,atol=1e-15)
    n=d[d['shift']==1]
    assert (n.n_fixed_pool_covered==n.n_draw_benchmark_covered).all()
    assert (n.n_fixed_pool_covered+n.n_excludes_zero==n.n_valid).all()
    a=d[(d['shift']==.5)&(d.n_cells_mature==50)].set_index(['cohort','pool'])
    for key,count in [(('smc','pooled'),173),(('smc','reference'),182),(('kul3','pooled'),174),(('kul3','reference'),169)]:
        assert a.loc[key,'n_fixed_pool_covered']==count
    assert not ((a.fixed_pool_coverage>=.9)&(a.exclusion_rate>=.8)).any()


def test_variance_formula_against_an_exact_discrete_distribution():
    from make_fixed_pool_tables import signal_to_noise
    # Source: P(X=0)=P(X=2)=1/2. Half thinning gives P(T=0,1,2)=(5/8,1/4,1/8).
    # Enumerate the independent arm distribution instead of duplicating the formula.
    terms=[(n,t,.5*p) for n in [0.,2.] for t,p in [(0.,.625),(1.,.25),(2.,.125)]]
    fixed_errors=np.array([.4*(t-n)+.2 for n,t,p in terms])
    draw_errors=np.array([.4*(t-.5*n) for n,t,p in terms])
    probabilities=np.array([p for n,t,p in terms])
    v_fixed=np.dot(probabilities,fixed_errors**2)-np.dot(probabilities,fixed_errors)**2
    v_draw=np.dot(probabilities,draw_errors**2)-np.dot(probabilities,draw_errors)**2
    assert np.isclose(v_fixed-v_draw,.4**2*(1-.5**2))
    assert np.isclose(v_fixed/.2**2,6.)
    assert np.isclose(signal_to_noise(1.,1.,1,n_n=1),.2/np.sqrt(v_fixed))


def test_source_pool_comparison_uses_matching_tissue_eligibility():
    from make_fixed_pool_tables import signal_to_noise
    p=pd.read_csv(PAPER/'diagnostics/source_pool_properties.csv')
    for c,n_pairs,total,cv,snr in [('smc',45,662,4.35,3.02),('kul3',15,844,16.99,1.53)]:
        q=p[(p.cohort==c)&(p.pool=='reference')]
        patients=10 if c=='smc' else 6
        assert len(q)==n_pairs and q.n_source_cells.sum()/(patients-1)==total
        assert round(q.squared_cv.median(),2)==cv
        values=[signal_to_noise(r.mean_expression,r.variance_expression,50) for r in q.itertuples()]
        assert round(float(np.median(values)),2)==snr


def test_anonymous_review_record_excludes_local_provenance():
    review=(PAPER/'evidence/review_record.json').read_text()
    from check_submission import source_graph
    text=review+'\n'+'\n'.join(p.read_text() for p in source_graph(PAPER/'main.tex'))
    for token in ['gate_memo_w2','harness_design_spec','BetterRiseProject','bodebosell','/Users/']:
        assert token.lower() not in text.replace(r'\_','_').lower()


def test_external_control_design_denominators_and_target_ranking():
    from make_external_control_table import render as render_external
    d=table('external_control_primary')
    assert len(d)==8*5*6
    assert (d.n_valid+d.n_excluded).eq(200).all()
    assert d.n_valid.sum()==47969 and d.n_excluded.sum()==31
    assert d.loc[d.n_excluded>0,'n_trial'].eq(50).all()
    assert d.loc[d.n_trial==1600,'n_valid'].eq(200).all()
    assert prose('external_control_table.tex')==render_external(d)
    selected=d[d.n_trial==1600].set_index(['reference','estimator'])
    for ref,winner in [('obs-pooled','ate-standardisation'),('po-pooled','ate-standardisation'),('po-trial','att-standardisation')]:
        pair=selected.loc[ref].loc[['ate-standardisation','att-standardisation']]
        assert pair.rmse.idxmin()==winner
    pooled=selected.loc[('po-trial','ate-standardisation')]
    assert round(pooled.bias,3)==.885
    assert round(pooled.estimate_median,3)==1.890
    assert round(pooled.reference_median,3)==1.003
    assert pooled.estimate_median>1.5>pooled.reference_median
    own=d[(d.estimator=='ate-standardisation')&(d.reference=='obs-pooled')]
    assert own.max_residual.eq(0).all()


def test_external_control_null_and_balance_qualify_claim():
    shifts=table('external_control_shift')
    null=shifts[(shifts.estimator=='ate-standardisation')&(shifts['shift']==0)].iloc[0]
    assert null.max_residual_vs_obs_pooled==0
    assert round(null.bias_vs_po_trial,3)==-.002
    assert not null.decision_flips
    b=table('external_control_balance')
    pooled=b[b.weighting=='ate']
    assert pooled.n_replicates.eq(50).all()
    assert pooled.max_abs_smd_between_arms.max()<1e-14
    assert round(pooled.max_abs_smd_vs_enrolled.max(),3)==.964
    assert not table('external_control_falsifiers').fired.any()


def test_external_control_provenance_is_pinned_and_does_not_inflate_replication():
    import hashlib
    from _tables import REPO_ROOT
    provenance=json.loads((PAPER/'diagnostics/external_control_provenance.json').read_text())
    assert provenance['seed_streams']==1 and provenance['primary_study_draws']==1200
    for name,digest in provenance['input_sha256'].items():
        assert hashlib.sha256((REPO_ROOT/name).read_bytes()).hexdigest()==digest
    assert 'One seed (20260915)' in prose('appendix.tex')
    assert '48{,}000 comparisons' in prose('appendix.tex')
    for name in ['primary','shift','balance','falsifiers']:
        exported=pd.read_csv(PAPER/f'diagnostics/external_control_{name}.csv')
        pd.testing.assert_frame_equal(exported,table(f'external_control_{name}'),check_exact=False,rtol=1e-12,atol=1e-15)


def test_merged_control_panel_claims_use_correct_strata():
    d=table('retained_control_strata')
    r=d[(d.gene=='MS4A12')&~d.degenerate_stratum]
    primary=r[r.cohort=='GSE178341']
    assert round(primary.rel_change_median.min(),3)==-.971
    assert round(primary.rel_change_median.max(),3)==-.951
    assert round(r.baseline_cp10k_normal.median(),2)==2.56
    assert round(r.cp10k_tumour.median(),3)==.104
    assert d[(d.gene=='MS4A12')&d.degenerate_stratum].n_tumour_mature_median.median()==1
    epithelial=r[r.label_selects_nothing]
    assert round(epithelial.rel_change_median.min(),3)==-1.000
    assert round(epithelial.rel_change_median.max(),3)==-.967
    bulk=table('retained_control_bulk').set_index('gene')
    assert round(bulk.loc['MS4A12','log2_fold_change'],2)==-8.18
    assert '-0.971' in prose('calibration.tex')
