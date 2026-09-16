"""Independent checks for the fixed-pool skew-aware interval follow-up."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import pytest
from scipy.stats import bootstrap
from repair_followup import bootstrap_moments,skew_intervals,METHODS
from precision_followup import percentile_ci,welch_ci
ROOT=Path(__file__).resolve().parent

@pytest.mark.parametrize('kind',['continuous','counts','sparse','zero_tumor'])
def test_bca_matches_scipy_two_sample_jackknife(kind):
    rng=np.random.default_rng(338)
    if kind=='continuous':normal=rng.gamma(2,3,800);tumor=rng.gamma(.7,2,30)
    elif kind=='counts':normal=rng.poisson(4,800).astype(float);tumor=rng.poisson(2,50).astype(float)
    else:
        normal=np.where(rng.random(800)<.65,0,rng.poisson(20,800)).astype(float)
        tumor=np.zeros(30) if kind=='zero_tumor' else np.where(rng.random(30)<.7,0,rng.poisson(8,30)).astype(float)
    seed=23
    actual,_=skew_intervals(normal,tumor,bootstrap_moments(normal,tumor,seed))
    expected=bootstrap((normal,tumor),lambda x,y,axis:.4*(y.mean(axis=axis)-x.mean(axis=axis)),n_resamples=2000,method='BCa',rng=np.random.default_rng(seed))
    np.testing.assert_allclose(actual['bca2000'],expected.confidence_interval,rtol=1e-11,atol=1e-11)


def test_original_percentile_and_welch_equivalence():
    rng=np.random.default_rng(27);normal=rng.poisson(4,800).astype(float);tumor=rng.poisson(3,30).astype(float)
    got,_=skew_intervals(normal,tumor,bootstrap_moments(normal,tumor,51,n_boot=200))
    np.testing.assert_allclose(got['percentile200'],percentile_ci(normal,tumor,51))
    np.testing.assert_allclose(got['welch'],welch_ci(normal,tumor))


def test_studentized_interval_uses_reversed_pivot_quantiles():
    normal=np.array([0.,1.,3.,8.]);tumor=np.array([1.,2.,10.])
    # Independently form every pair of leave-one-out resamples.
    ns=[np.delete(normal,i) for i in range(len(normal))];ts=[np.delete(tumor,i) for i in range(len(tumor))]
    mn=[];vn=[];mt=[];vt=[];pivots=[]
    point=.4*(tumor.mean()-normal.mean())
    se=.4*np.sqrt(normal.var(ddof=1)/len(normal)+tumor.var(ddof=1)/len(tumor))
    for x in ns:
        for y in ts:
            mn.append(x.mean());vn.append(x.var(ddof=1));mt.append(y.mean());vt.append(y.var(ddof=1))
            sb=.4*np.sqrt(x.var(ddof=1)/len(normal)+y.var(ddof=1)/len(tumor))
            pivots.append((.4*(y.mean()-x.mean())-point)/sb)
    moments=[tuple(map(np.array,(mn,vn))),tuple(map(np.array,(mt,vt)))]
    got,_=skew_intervals(normal,tumor,moments)
    lo,hi=np.quantile(pivots,[.025,.975]);np.testing.assert_allclose(got['bootstrap_t2000'],[point-hi*se,point-lo*se])


def test_affine_invariance_and_exchange_of_arms():
    rng=np.random.default_rng(911);normal=rng.gamma(1,3,50);tumor=rng.gamma(.3,2,30)
    moments=bootstrap_moments(normal,tumor,7);base,_=skew_intervals(normal,tumor,moments)
    transformed=[(3*m+4,9*v) for m,v in moments]
    scaled,_=skew_intervals(3*normal+4,3*tumor+4,transformed)
    swapped,_=skew_intervals(tumor,normal,moments[::-1])
    for method in METHODS:
        np.testing.assert_allclose(scaled[method],3*base[method],rtol=1e-10,atol=1e-10)
        np.testing.assert_allclose(swapped[method],-base[method][::-1],rtol=1e-10,atol=1e-10)


def test_degenerate_methods_are_explicitly_unavailable():
    normal=np.zeros(800);tumor=np.zeros(30)
    with np.errstate(invalid='ignore',divide='ignore'):
        got,_=skew_intervals(normal,tumor,bootstrap_moments(normal,tumor,4))
    for method in ['bca2000','bootstrap_t2000','welch']:assert np.isnan(got[method]).all()
    np.testing.assert_array_equal(got['percentile2000'],[0,0])


def test_saved_repair_results_and_provenance():
    if not (ROOT/'diagnostics/repair_followup.csv').exists():pytest.skip('Run study after implementation checks')
    d=pd.read_csv(ROOT/'diagnostics/repair_followup.csv');assert len(d)==80
    assert set(d.method)==set(METHODS);assert d.n_attempted.eq(2000).all()
    assert (d.n_valid+d.n_excluded==d.n_attempted).all()
    assert d.n_attempted.sum()==160000
    for metric,count in [('coverage','n_covered'),('rejection','n_rejected')]:np.testing.assert_allclose(d[metric],d[count]/d.n_valid)
    null=d[d['shift']==1];assert (null.n_covered+null.n_rejected==null.n_valid).all()
    seeds=pd.read_csv(ROOT/'diagnostics/repair_seed_rates.csv');assert len(seeds)==400
    counts=seeds.groupby(['cohort','pool','n_cells_mature','shift','method'])[['n_attempted','n_valid','n_covered','n_rejected']].sum()
    pd.testing.assert_frame_equal(counts,d.set_index(counts.index.names)[counts.columns].sort_index())
    prov=json.loads((ROOT/'diagnostics/repair_provenance.json').read_text())
    for file,digest in prov['input_sha256'].items():assert hashlib.sha256((ROOT.parents[1]/file).read_bytes()).hexdigest()==digest
    for file,digest in prov['output_sha256'].items():assert hashlib.sha256((ROOT/'diagnostics'/file).read_bytes()).hexdigest()==digest


def test_rendered_table_and_reported_repair_claims():
    from make_repair_results import render_table
    d=pd.read_csv(ROOT/'diagnostics/repair_followup.csv')
    assert (ROOT/'sections/repair_table.tex').read_text()==render_table(d)
    null=d[(d.n_cells_mature==50)&(d['shift']==1)]
    assert (null.rejection_low>.05).all()
    pair=pd.read_csv(ROOT/'diagnostics/repair_paired.csv')
    budget=pair[(pair.n_cells_mature==50)&(pair['shift']==1)&(pair.metric=='rejection')&(pair.comparison=='percentile200_minus_percentile2000')]
    np.testing.assert_allclose(sorted(budget.difference),[.005,.006,.0085,.01])
    a=d[(d['shift']==.5)&(d.method=='percentile200')&(d.n_cells_mature==30)]
    np.testing.assert_allclose([a.coverage.min(),a.coverage.max(),a.rejection.min(),a.rejection.max()],[.781,.8905,.485,.6785])
    numerical=pd.read_csv(ROOT/'diagnostics/repair_numerical.csv')
    rate=numerical.bca_tail_under_10/numerical.n_studies
    np.testing.assert_allclose([rate.min(),rate.max()],[.052,.563])
