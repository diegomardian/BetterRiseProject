"""Independent checks on the selected incoming evidence and its interpretation."""
from pathlib import Path
import hashlib,json,itertools
import numpy as np
import pandas as pd
from audit_incoming_interval import check_aggregates
ROOT=Path(__file__).resolve().parent


def read(name):
    d=pd.read_csv(ROOT/'diagnostics'/f'incoming_{name}.csv')
    return d.rename(columns={'draw_reference_inclusion':'coverage','draw_reference_inclusion_mcse':'coverage_mcse','source_verdict':'verdict'})


def test_incoming_count_reaggregation_and_provenance():
    names=['interval_diagnosis_rates','interval_repair_rates','interval_repair_pooled','interval_repair_by_seed']
    report=check_aggregates({n:read(n) for n in names})
    assert report['source_null_tolerance_passes']==46
    prov=json.loads((ROOT/'diagnostics/incoming_interval_provenance.json').read_text())
    for name,sha in prov['aggregate_exports_sha256'].items():
        assert hashlib.sha256((ROOT/'diagnostics'/name).read_bytes()).hexdigest()==sha


def test_reported_substitution_rates_and_nonadditivity():
    d=read('interval_diagnosis_rates');n=d[d.n_cells_mature==50]
    ranges={('empirical','fixed_fraction'):(.097,.318),('zeros_removed_mean_matched','fixed_fraction'):(.071,.1265),('gaussian_matched','fixed_fraction'):(.0545,.0545),('empirical','balanced'):(.0745,.1145)}
    for (family,design),expected in ranges.items():
        x=n[(n.family==family)&(n.design==design)]
        np.testing.assert_allclose([x.null_rejection.min(),x.null_rejection.max()],expected)
    a=read('interval_diagnosis_attribution')
    assert (a[a.n_cells_mature==50].closure_verdict=='over_explained').all()
    assert not a[a.n_cells_mature==50].diagnosis_closes.any()
    assert not n[n.family=='skew_matched_no_zeros'].skew_matched.all()


def test_patient_method_tradeoff_is_not_universal_abstention():
    d=read('interval_repair_rates');x=d[(d.n_cells_mature==50)&(d.n_held_out==2)&(d.candidate=='pseudobulk_patient_t')]
    assert len(x)==4
    np.testing.assert_allclose([x.null_rejection.min(),x.null_rejection.max()],[55/1985,84/1335])
    np.testing.assert_allclose([x.discrimination.min(),x.discrimination.max()],[84/1993,105/1309])
    np.testing.assert_allclose([x.alt_abstention_rate.min(),x.alt_abstention_rate.max()],[.0035,.3455])
    assert (d[(d.n_cells_mature==50)&(d.n_held_out==5)].null_abstention_rate==0).all()


def test_equal_arm_symmetry_does_not_require_symmetric_source():
    source=np.array([0.,0.,0.,4.])
    assert np.mean((source-source.mean())**3)>0
    means=np.array([np.mean(x) for x in itertools.product(source,repeat=2)])
    differences=(means[:,None]-means).ravel()
    np.testing.assert_array_equal(np.sort(differences),np.sort(-differences))
    other=np.array([np.mean(x) for x in itertools.product(source,repeat=3)])
    asymmetric=(means[:,None]-other).ravel()
    assert np.mean(asymmetric**3)>0


def test_exact_stratum_balance_collapses_reference_weights():
    g=np.repeat([0,1,2],[4,6,8]);a=np.concatenate([np.tile([0,1],n) for n in [2,3,4]])
    y=3*g+(1+g)*a+np.sin(np.arange(len(g)))
    contrasts=np.array([y[(g==k)&(a==1)].mean()-y[(g==k)&(a==0)].mean() for k in range(3)])
    standardized=np.dot(np.bincount(g)/len(g),contrasts)
    design=np.column_stack([a,*[g==k for k in range(3)]])
    ols=np.linalg.lstsq(design,y,rcond=None)[0][0]
    np.testing.assert_allclose([ols,y[a==1].mean()-y[a==0].mean()],standardized,atol=1e-12)
