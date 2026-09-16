"""Independent numerical checks for the added interval comparison."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from precision_followup import percentile_ci,welch_ci
from src.harness.interval import within_patient_intrinsic_ci
ROOT=Path(__file__).resolve().parent


def test_percentile_matches_original_including_sparse_samples():
    rng=np.random.default_rng(812)
    for n in (5,50,100,800):
        normal=rng.poisson(3,size=800).astype(float)
        for tumour in [np.zeros(n),rng.poisson(2,size=n).astype(float)]:
            for seed in [0,918]:
                expected=within_patient_intrinsic_ci(normal,tumour,frac_mature_normal=.4,frac_mature_tumour=n/2000,n_boot=200,seed=seed,weighting='normal')
                np.testing.assert_allclose(percentile_ci(normal,tumour,seed),expected,atol=1e-14)


def test_welch_matches_independent_scipy_implementation():
    rng=np.random.default_rng(441)
    for n in (5,50,100,800):
        normal=rng.gamma(2,3,size=800)
        for tumour in [np.zeros(n),rng.gamma(.2,4,size=n)]:
            expected=ttest_ind(tumour,normal,equal_var=False).confidence_interval()
            np.testing.assert_allclose(welch_ci(normal,tumour),.4*np.array(expected),rtol=1e-13)
    assert np.isnan(welch_ci(np.zeros(10),np.zeros(5))).all()


def test_followup_denominators_and_uncertainty():
    d=pd.read_csv(ROOT/'diagnostics/precision_followup.csv')
    assert len(d)==64
    assert d.n_attempted.eq(5000).all()
    assert (d.n_valid+d.n_excluded==d.n_attempted).all()
    assert d.n_attempted.sum()==320000
    for metric,count in [('coverage','n_covered'),('draw_coverage','n_draw_covered'),('rejection','n_rejected')]:
        np.testing.assert_allclose(d[metric],d[count]/d.n_valid)
        assert (d[metric+'_low']<=d[metric]).all() and (d[metric]<=d[metric+'_high']).all()
    null=d[d['shift']==1]
    assert (null.n_covered+null.n_rejected==null.n_valid).all()
    assert (null.n_covered==null.n_draw_covered).all()
    np.testing.assert_allclose(d.rejection,d.negative_rejection+d.positive_rejection)
    s=pd.read_csv(ROOT/'diagnostics/precision_seed_rates.csv')
    assert len(s)==320 and s.seed.nunique()==5
    for keys,g in s.groupby(['cohort','pool','n_cells_mature','shift','method']):
        row=d.set_index(['cohort','pool','n_cells_mature','shift','method']).loc[keys]
        assert g.n_valid.sum()==row.n_valid
        for metric in ['coverage','rejection']:
            assert np.isclose(np.average(g[metric],weights=g.n_valid),row[metric])


def test_followup_paired_differences_and_provenance():
    d=pd.read_csv(ROOT/'diagnostics/precision_followup.csv').set_index(['cohort','pool','n_cells_mature','shift','method'])
    paired=pd.read_csv(ROOT/'diagnostics/precision_paired.csv')
    for row in paired.itertuples():
        key=(row.cohort,row.pool,row.n_cells_mature,row.shift)
        a=d.loc[key+('percentile200',)]
        if row.comparison=='welch_minus_percentile200':
            b=d.loc[key+('welch',)]
            if row.n_pairs==a.n_valid==b.n_valid:
                assert np.isclose(row.difference,b[row.metric]-a[row.metric])
        else:assert np.isclose(row.difference,a.draw_coverage-a.coverage)
        assert 0<=row.mc_se<=1/np.sqrt(row.n_pairs)
    p=json.loads((ROOT/'diagnostics/precision_provenance.json').read_text())
    assert p['study_draws']==160000 and p['source_pair_checks']==120
    for path,hash_ in p['input_sha256'].items():
        assert hashlib.sha256((ROOT.parents[1]/path).read_bytes()).hexdigest()==hash_
    for name,hash_ in p['output_sha256'].items():
        assert hashlib.sha256((ROOT/'diagnostics'/name).read_bytes()).hexdigest()==hash_


def test_generated_table_and_new_headline_claims_match_results():
    from make_precision_results import render_table
    d=pd.read_csv(ROOT/'diagnostics/precision_followup.csv')
    assert (ROOT/'sections/precision_table.tex').read_text()==render_table(d)
    a=d[(d.n_cells_mature==50)&(d['shift']==.5)&(d.method=='percentile200')]
    assert (a.coverage_high<.9).sum()==3
    passing=a[a.coverage_low>.9]
    assert len(passing)==1 and passing.iloc[0].cohort=='smc' and passing.iloc[0].pool=='reference'
    for method,low,high in [('percentile200',9.8,32.2),('welch',8.9,30.8)]:
        null=d[(d.n_cells_mature==50)&(d['shift']==1)&(d.method==method)]
        assert round(100*null.rejection.min(),1)==low
        assert round(100*null.rejection.max(),1)==high
        assert (null.rejection_low>.05).all()
        effect=d[(d.n_cells_mature==50)&(d['shift']==.5)&(d.method==method)]
        assert (effect.rejection_high<.8).all()
        assert f'{low:.1f}--{high:.1f}\\%' in (ROOT/'sections/abstract.tex').read_text()
    k=d[(d.cohort=='kul3')&(d.pool=='pooled')&(d['shift']==.5)&(d.method=='percentile200')].set_index('n_cells_mature')
    assert k.loc[100].rejection_high<k.loc[50].rejection_low
    assert k.loc[800].rejection_low>k.loc[100].rejection_high
    assert k.loc[100].all_zero_fraction==0
