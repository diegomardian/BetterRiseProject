"""Fixed-pool interval comparison at the reporting boundary and flag-wide band."""
from pathlib import Path
import sys,itertools,json,hashlib,platform,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import pandas as pd
from scipy.special import ndtr,ndtri
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from precision_followup import welch_ci
from make_fixed_pool_tables import wilson
from src.harness.calibration_gap import load_cohort_arrays,_pool_mask
SEEDS=tuple(range(202609170,202609175)); REPS=400; B=2000
COUNTS=(30,50); METHODS=('percentile200','percentile2000','bca2000','bootstrap_t2000','welch')


def bootstrap_moments(normal,tumor,seed,n_boot=B):
    rng=np.random.default_rng(seed);out=[]
    for values in (normal,tumor):
        draws=values[rng.integers(len(values),size=(n_boot,len(values)))]
        out.append((draws.mean(axis=1),draws.var(axis=1,ddof=1)))
    return out


def skew_intervals(normal,tumor,moments):
    (mn,vn),(mt,vt)=moments
    boot=.4*(mt-mn);point=.4*(tumor.mean()-normal.mean())
    invalid=np.array([np.nan,np.nan])
    intervals={'percentile200':np.quantile(boot[:200],[.025,.975]),
               'percentile2000':np.quantile(boot,[.025,.975]),'welch':welch_ci(normal,tumor)}
    u=np.concatenate((-(normal-normal.mean())/len(normal),(tumor-tumor.mean())/len(tumor)))
    denom=6*np.sum(u*u)**1.5
    a=np.sum(u**3)/denom if denom>0 else np.nan
    prob=(np.count_nonzero(boot<point)+np.count_nonzero(boot<=point))/(2*len(boot))
    z0=ndtri(prob);z=ndtri(np.array([.025,.975]));div=1-a*(z0+z)
    with np.errstate(invalid='ignore',divide='ignore'):
        adjusted=ndtr(z0+(z0+z)/div)
    bca_ok=np.isfinite(adjusted).all() and np.isfinite(z0) and np.isfinite(a) and (div>0).all() and adjusted[0]<=adjusted[1]
    intervals['bca2000']=np.quantile(boot,adjusted) if bca_ok else invalid.copy()
    se=.4*np.sqrt(normal.var(ddof=1)/len(normal)+tumor.var(ddof=1)/len(tumor))
    se_star=.4*np.sqrt(vn/len(normal)+vt/len(tumor))
    bad_pivots=int(np.count_nonzero((se_star<=0)|~np.isfinite(se_star)))
    if np.isfinite(se) and se>0 and bad_pivots==0:
        quantiles=np.quantile((boot-point)/se_star,[.025,.975])
        intervals['bootstrap_t2000']=point-se*quantiles[::-1]
    else:intervals['bootstrap_t2000']=invalid.copy()
    return intervals,dict(bca_acceleration=a,bca_min_tail_count=float(len(boot)*min(adjusted[0],1-adjusted[1])) if bca_ok else np.nan,bad_pivots=bad_pivots)


def run_condition(job):
    cohort,pool,n,shift,pairs=job;rows=[];diagnostics=[]
    for stream in SEEDS:
        rng=np.random.default_rng([stream,int(cohort=='kul3'),int(pool=='reference'),n,int(100*shift)])
        for rep in range(REPS):
            source=pairs[int(rng.integers(len(pairs)))];normal=rng.choice(source,800).astype(float)
            tumor=rng.choice(source,n)
            if shift!=1:tumor=rng.binomial(tumor,shift)
            tumor=tumor.astype(float);seed=int(rng.integers(2**63));truth=.4*(shift-1)*source.mean()
            intervals,diag=skew_intervals(normal,tumor,bootstrap_moments(normal,tumor,seed));diagnostics.append(diag)
            for method in METHODS:
                lo,hi=intervals[method];valid=bool(np.isfinite([lo,hi]).all() and lo<=hi)
                rows.append(dict(seed=stream,replicate=rep,method=method,valid=valid,
                    covered=bool(valid and lo<=truth<=hi),rejected=bool(valid and (hi<0 or lo>0)),width=hi-lo))
        print(cohort,pool,n,shift,'stream',stream,'complete',flush=True)
    d=pd.DataFrame(rows);keys=dict(cohort=cohort,pool=pool,n_cells_mature=n,shift=shift)
    summary=[];streams=[];paired=[]
    for method,g in d.groupby('method',sort=False):
        v=g[g.valid];row=keys|dict(method=method,n_attempted=len(g),n_valid=len(v),n_excluded=len(g)-len(v),median_width=float(v.width.median()))
        for metric,flag in [('coverage','covered'),('rejection','rejected')]:
            k=int(v[flag].sum());row['n_'+flag]=k;row[metric]=k/len(v) if len(v) else np.nan
            low,high=wilson(k,len(v)) if len(v) else (np.nan,np.nan);row[metric+'_low']=low;row[metric+'_high']=high
        summary.append(row)
        for seed,h in g.groupby('seed'):
            valid=h[h.valid];streams.append(keys|dict(method=method,seed=int(seed),n_attempted=len(h),n_valid=len(valid),n_covered=int(valid.covered.sum()),n_rejected=int(valid.rejected.sum())))
    baseline=d[d.method=='percentile2000'].set_index(['seed','replicate'])
    for method in METHODS:
        if method=='percentile2000':continue
        other=d[d.method==method].set_index(['seed','replicate']);mask=baseline.valid&other.valid
        for metric,flag in [('coverage','covered'),('rejection','rejected')]:
            delta=other.loc[mask,flag].astype(int)-baseline.loc[mask,flag].astype(int)
            paired.append(keys|dict(comparison=method+'_minus_percentile2000',metric=metric,n_pairs=len(delta),difference=float(delta.mean()),mc_se=float(delta.std(ddof=1)/np.sqrt(len(delta)))))
    diag=pd.DataFrame(diagnostics)
    numerical=keys|dict(n_studies=len(diag),bca_undefined=int(diag.bca_min_tail_count.isna().sum()),bca_tail_under_10=int((diag.bca_min_tail_count<10).sum()),median_bca_min_tail_count=float(diag.bca_min_tail_count.median()),studies_with_bad_pivots=int((diag.bad_pivots>0).sum()))
    return summary,streams,paired,numerical


def main():
    start=time.monotonic();jobs=[];pools=[]
    for cohort in ('smc','kul3'):
        a=load_cohort_arrays(cohort)
        for pool in ('pooled','reference'):
            mask=_pool_mask(pd.Series(a['tissue']),pool);counts=np.asarray(a['counts'])[mask,a['genes'].index('GUCA2A')]
            patient=np.asarray(a['patient_id'])[mask];mature=np.asarray(a['cell_type'])[mask]=='differentiated';pairs=[]
            for number,pair in enumerate(itertools.combinations(sorted(set(patient)),2),1):
                source=counts[np.isin(patient,pair)&mature];assert len(source)>0 and source.mean()>0;pairs.append(source)
                pools.append(dict(cohort=cohort,pool=pool,pair_number=number,n_source_cells=len(source),mean_expression=source.mean(),variance_expression=source.var()))
            for n in COUNTS:
                for shift in (.5,1.):jobs.append((cohort,pool,n,shift,pairs))
    old=pd.read_csv(ROOT/'diagnostics/source_pool_properties.csv');new=pd.DataFrame(pools)
    for col in new:
        if col in ('cohort','pool'):assert new[col].tolist()==old[col].tolist()
        else:np.testing.assert_allclose(new[col],old[col],rtol=1e-12)
    print('PASS: all 120 source-pool summaries match',flush=True)
    results=[[],[],[],[]]
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures=[executor.submit(run_condition,j) for j in jobs]
        for f in as_completed(futures):
            result=f.result()
            for i,r in enumerate(result):results[i].extend(r if isinstance(r,list) else [r])
    paths=[]
    for name,records in zip(('repair_followup','repair_seed_rates','repair_paired','repair_numerical'),results):
        path=ROOT/'diagnostics'/f'{name}.csv';pd.DataFrame(records).sort_values(['cohort','pool','n_cells_mature','shift']).to_csv(path,index=False);paths.append(path)
    inputs=[Path(__file__),ROOT/'REPAIR_FOLLOWUP_DESIGN.md',ROOT/'precision_followup.py',ROOT/'diagnostics/source_pool_properties.csv',REPO/'src/harness/calibration_gap.py',REPO/'src/harness/pseudobulk.py',REPO/'src/reference/labels.py']
    import scipy
    provenance=dict(purpose='Exploratory fixed-pool reporting-band and skew-aware interval comparison',source_revision='bf6b825',counts=COUNTS,shifts=[.5,1.],methods=METHODS,seed_streams=SEEDS,replicates_per_stream=REPS,outer_per_setting=2000,settings=16,study_draws=32000,interval_evaluations=160000,inner_bootstrap=B,source_pair_checks=120,elapsed_seconds=time.monotonic()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,input_sha256={str(f.relative_to(REPO)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},output_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in paths})
    (ROOT/'diagnostics/repair_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('DONE',round(time.monotonic()-start,1),'seconds',flush=True)
if __name__=='__main__':main()
