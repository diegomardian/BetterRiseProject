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
def prose(name):return (PAPER/'sections'/name).read_text()

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
    a=prose('appendix.tex')
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
    assert 'Thirteen generator--estimator combinations' in prose('appendix.tex')


def test_model_reference_uncertainty_is_only_monte_carlo_for_mixtures():
    d=table('trial_learned_generator_primary');d=d[d.reference=='T_model']
    mc=d[d.T_model_is_monte_carlo]
    assert mc.generator.str.startswith('gmm').all()
    assert round(mc.T_model_mc_se.min(),3)==.019
    assert round(mc.T_model_mc_se.max(),3)==.027
    assert (d[~d.T_model_is_monte_carlo].T_model_mc_se==0).all()
    assert 'errors are 0.019--0.027' in prose('appendix.tex')


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
    for x in [r'75.5\%',r'81.0\%','two of eight',r'45.0--55.0\%']:
        assert x in prose('calibration.tex')


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


def test_clean_control_table_matches_every_saved_entry():
    d=table('residual_performance_clean_control')
    labels={'empirical-mean-calibrated':'Full-sample mean','same-point-narrow-interval':'Same point, narrow interval','known-bias-plus-0.5':'Mean plus $0.5$','independent-half-mean':'Half-sample mean'}
    text=prose('blind.tex')
    assert len(d)==8 and (d.n_valid==2000).all()
    for r in d.itertuples():
        residual='0' if r.max_residual_vs_reference==0 else f'{r.max_residual_vs_reference:.4f}'
        row=f'{labels[r.estimator]} & ${residual}$ & ${r.bias:.4f}$ & ${r.rmse:.4f}$ & ${100*r.interval_coverage:.2f}\\%$'
        assert row in text
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
    assert 'neither a representative prevalence estimate' in prose('refdesign.tex')
    assert 'already the norm' not in '\n'.join(x.read_text() for x in (PAPER/'sections').glob('*.tex'))
