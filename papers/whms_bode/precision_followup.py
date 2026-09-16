"""Fresh fixed-design precision/interval study. Exports derived summaries only.

The sampled mature counts have the exact marginal law of generate_pseudobulk
at the fixed fractions; no non-mature cells contribute to the contrast.
"""
from pathlib import Path
import sys, itertools, json, hashlib, platform, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.stats import t
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from src.harness.calibration_gap import load_cohort_arrays, _pool_mask
from src.harness.interval import within_patient_intrinsic_ci
from make_fixed_pool_tables import wilson
SEEDS=tuple(range(202609160,202609165))
COUNTS=(5,50,100,800)
REPS=1000


def percentile_ci(normal,tumour,seed):
    rng=np.random.default_rng(seed)
    bn=normal[rng.integers(len(normal),size=(200,len(normal)))].mean(axis=1)
    bt=tumour[rng.integers(len(tumour),size=(200,len(tumour)))].mean(axis=1)
    return np.quantile(.4*(bt-bn),[.025,.975])


def welch_ci(normal,tumour):
    vn=normal.var(ddof=1)/len(normal);vt=tumour.var(ddof=1)/len(tumour)
    variance=vn+vt
    if not np.isfinite(variance) or variance<=0:return np.array([np.nan,np.nan])
    df=variance**2/(vn**2/(len(normal)-1)+vt**2/(len(tumour)-1))
    half=.4*t.ppf(.975,df)*np.sqrt(variance)
    estimate=.4*(tumour.mean()-normal.mean())
    return np.array([estimate-half,estimate+half])


def run_condition(job):
    cohort,pool,n,shift,pairs=job
    rows=[]
    for stream in SEEDS:
        # Count, shift and stream identifiers define independent draws per setting.
        rng=np.random.default_rng([stream,0 if cohort=='smc' else 1,0 if pool=='pooled' else 1,n,int(shift*100)])
        for rep in range(REPS):
            pair=int(rng.integers(len(pairs)));source=pairs[pair]
            normal=rng.choice(source,size=800).astype(float)
            tumour=rng.choice(source,size=n)
            if shift!=1: tumour=rng.binomial(tumour,shift)
            tumour=tumour.astype(float)
            boot_seed=int(rng.integers(0,2**63))
            mu=source.mean();variance=source.var()
            truth=.4*(shift-1)*mu;draw_truth=.4*(shift-1)*normal.mean()
            true_se=.4*np.sqrt(variance/800+(shift**2*variance+shift*(1-shift)*mu)/n)
            estimated_se=.4*np.sqrt(normal.var(ddof=1)/800+tumour.var(ddof=1)/n)
            base=dict(seed=stream,replicate=rep,pair=pair,truth=truth,draw_truth=draw_truth,
                      estimate=.4*(tumour.mean()-normal.mean()),
                      tumour_all_zero=bool((tumour==0).all()),se_ratio=estimated_se/true_se)
            for method,bounds in [('percentile200',percentile_ci(normal,tumour,boot_seed)),('welch',welch_ci(normal,tumour))]:
                lo,hi=bounds;valid=bool(np.isfinite(bounds).all())
                rows.append(base|dict(method=method,valid=valid,width=hi-lo,
                    covered=bool(valid and lo<=truth<=hi),draw_covered=bool(valid and lo<=draw_truth<=hi),
                    reject_negative=bool(valid and hi<0),reject_positive=bool(valid and lo>0)))
    d=pd.DataFrame(rows);summary=[];seed_summary=[];paired=[]
    keys=dict(cohort=cohort,pool=pool,n_cells_mature=n,shift=shift)
    for method,g in d.groupby('method'):
        v=g[g.valid];reject=v.reject_negative|v.reject_positive
        record=keys|dict(method=method,n_attempted=len(g),n_valid=len(v),n_excluded=len(g)-len(v),
            n_covered=int(v.covered.sum()),n_draw_covered=int(v.draw_covered.sum()),n_rejected=int(reject.sum()),
            coverage=float(v.covered.mean()),draw_coverage=float(v.draw_covered.mean()),rejection=float(reject.mean()),
            negative_rejection=float(v.reject_negative.mean()),positive_rejection=float(v.reject_positive.mean()),
            median_width=float(v.width.median()),median_se_ratio=float(v.se_ratio.median()),
            all_zero_fraction=float(v.tumour_all_zero.mean()),
            rejection_se_under_half=float(reject[v.se_ratio<.5].mean()) if (v.se_ratio<.5).any() else None,
            fraction_se_under_half=float((v.se_ratio<.5).mean()))
        for metric,col in [('coverage','covered'),('draw_coverage','draw_covered')]:
            low,high=wilson(int(v[col].sum()),len(v));record[metric+'_low']=low;record[metric+'_high']=high
        low,high=wilson(int(reject.sum()),len(v));record['rejection_low']=low;record['rejection_high']=high
        summary.append(record)
        for seed,h in g.groupby('seed'):
            h=h[h.valid]
            seed_summary.append(keys|dict(method=method,seed=int(seed),n_valid=len(h),coverage=float(h.covered.mean()),rejection=float((h.reject_negative|h.reject_positive).mean())))
    a=d[d.method=='percentile200'].set_index(['seed','replicate']);b=d[d.method=='welch'].set_index(['seed','replicate'])
    valid=a.valid&b.valid
    for name,x,y in [('coverage',a.covered,b.covered),('rejection',a.reject_negative|a.reject_positive,b.reject_negative|b.reject_positive)]:
        delta=y[valid].astype(int)-x[valid].astype(int)
        paired.append(keys|dict(comparison='welch_minus_percentile200',metric=name,n_pairs=len(delta),difference=float(delta.mean()),mc_se=float(delta.std(ddof=1)/np.sqrt(len(delta)))))
    delta=a.draw_covered[a.valid].astype(int)-a.covered[a.valid].astype(int)
    paired.append(keys|dict(comparison='draw_minus_fixed_percentile200',metric='coverage',n_pairs=len(delta),difference=float(delta.mean()),mc_se=float(delta.std(ddof=1)/np.sqrt(len(delta)))))
    print(cohort,pool,n,shift,'complete',flush=True)
    return summary,seed_summary,paired


def main():
    start=time.monotonic();jobs=[];pool_summaries=[]
    for cohort in ('smc','kul3'):
        a=load_cohort_arrays(cohort)
        for pool in ('pooled','reference'):
            mask=_pool_mask(pd.Series(a['tissue']),pool)
            counts=np.asarray(a['counts'])[mask,a['genes'].index('GUCA2A')]
            patient=np.asarray(a['patient_id'])[mask];mature=np.asarray(a['cell_type'])[mask]=='differentiated'
            pairs=[]
            for number,pair in enumerate(itertools.combinations(sorted(set(patient)),2),1):
                source=counts[np.isin(patient,pair)&mature];assert len(source)>0 and source.mean()>0
                pairs.append(source)
                pool_summaries.append(dict(cohort=cohort,pool=pool,pair_number=number,n_source_cells=len(source),mean_expression=source.mean(),variance_expression=source.var()))
            for n in COUNTS:
                for shift in (.5,1.):jobs.append((cohort,pool,n,shift,pairs))
    old=pd.read_csv(ROOT/'diagnostics/source_pool_properties.csv');new=pd.DataFrame(pool_summaries)
    for col in new:
        if col in ['cohort','pool']:assert new[col].tolist()==old[col].tolist()
        else:np.testing.assert_allclose(new[col],old[col],rtol=1e-12)
    # Exact interval implementation check using actual source distributions.
    for ci,job in enumerate(jobs):
        _,_,n,s,pairs=job;rng=np.random.default_rng(ci)
        normal=rng.choice(pairs[0],800).astype(float);tumour=rng.binomial(rng.choice(pairs[0],n),s).astype(float)
        np.testing.assert_allclose(percentile_ci(normal,tumour,ci),within_patient_intrinsic_ci(normal,tumour,frac_mature_normal=.4,frac_mature_tumour=n/2000,n_boot=200,seed=ci,weighting='normal'),rtol=1e-13,atol=1e-13)
    print('PASS: 120 source pools and 32 original-interval comparisons',flush=True)
    summaries=[];streams=[];paired=[]
    with ProcessPoolExecutor(max_workers=4) as executor:
        for result in executor.map(run_condition,jobs):
            summaries.extend(result[0]);streams.extend(result[1]);paired.extend(result[2])
    out=ROOT/'diagnostics'
    for name,records in [('precision_followup',summaries),('precision_seed_rates',streams),('precision_paired',paired)]:
        pd.DataFrame(records).to_csv(out/f'{name}.csv',index=False)
    paths=[Path(__file__),ROOT/'FOLLOWUP_DESIGN.md',REPO/'src/harness/pseudobulk.py',REPO/'src/harness/interval.py',REPO/'src/harness/calibration_gap.py',REPO/'src/reference/labels.py',ROOT/'diagnostics/source_pool_properties.csv']
    prov=dict(design_revision='after ccb6264',purpose='Exploratory precision and paired interval comparison on existing empirical pools',seed_streams=SEEDS,replicates_per_seed=REPS,outer_per_setting=5000,counts=COUNTS,shifts=[.5,1.],settings=32,study_draws=160000,interval_evaluations=320000,inner_bootstrap=200,original_interval_checks=32,source_pair_checks=120,sampling='Exact marginal mature-count law; fresh draws, not bitwise replay',python=platform.python_version(),numpy=np.__version__,elapsed_seconds=time.monotonic()-start,input_sha256={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('precision_*.csv')})
    (out/'precision_provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
    print('DONE',round(time.monotonic()-start,1),'seconds',flush=True)

if __name__=='__main__':main()
