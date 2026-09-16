"""Rescore the unchanged incoming patient-t intervals against two fixed targets."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,hashlib,json,subprocess,sys,platform
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
SEEDS=(20260831,1,2,3,4)


def fixed_targets(values,patients,shift):
    labels=np.unique(patients)
    cell=.4*(shift-1)*values.mean()
    patient=.4*(shift-1)*np.mean([values[patients==p].mean() for p in labels])
    return cell,patient,len(labels)


def wilson(k,n):
    if n==0:return np.nan,np.nan
    z=1.959963984540054;p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-half,mid+half


def run(job):
    from src.harness.run_interval_repair import _load_cached,_replicate_seed,_composition
    from src.harness.pseudobulk import patient_holdout,generate_pseudobulk
    from src.harness.interval_repair import pseudobulk_patient_t_interval
    cohort,pool,holdout,count,shift,cache=job;a=_load_cached(Path(cache),pool)
    values=a['counts'][:,a['genes'].index('GUCA2A')];types=sorted(set(a['cell_type']))
    comp_n=_composition(.4,types,'differentiated');comp_t=_composition(count/2000,types,'differentiated')
    rows=[]
    for stream in SEEDS:
        for rep in range(400):
            seed=_replicate_seed(stream,count,shift,rep)
            _,held=patient_holdout(a['patient_id'],n_held_out=holdout,seed=seed)
            mask=(a['cell_type']=='differentiated')&np.isin(a['patient_id'],list(held))
            tc,tp,n_source=fixed_targets(values[mask],a['patient_id'][mask],shift)
            sample=generate_pseudobulk(a['counts'],a['cell_type'],a['patient_id'],a['genes'],composition_normal=comp_n,composition_tumour=comp_t,shift={'GUCA2A':shift},held_out_patients=held,n_cells=2000,seed=seed,mature_label='differentiated')
            cells=sample.mature_expression['GUCA2A'];normal,tumor=cells['normal'],cells['tumour']
            pn=sample.drawn_patient_id['normal'][sample.drawn_is_mature['normal']]
            pt=sample.drawn_patient_id['tumour'][sample.drawn_is_mature['tumour']]
            bounds=pseudobulk_patient_t_interval(normal,tumor,pn,pt,frac_mature_normal=.4)
            drawn=sample.truth.parametric['GUCA2A']['normal']['intrinsic']
            assert np.isclose(drawn,.4*(shift-1)*normal.mean())
            valid=bounds[0] is not None and bounds[1] is not None
            lo,hi=bounds if valid else (np.nan,np.nan)
            rows.append(dict(seed=stream,valid=valid,rejected=bool(valid and (hi<0 or lo>0)),cell_covered=bool(valid and lo<=tc<=hi),patient_covered=bool(valid and lo<=tp<=hi),draw_covered=bool(valid and lo<=drawn<=hi),width=hi-lo,n_source_patients=n_source,abs_target_gap=abs(tc-tp)))
    d=pd.DataFrame(rows);out=[]
    keys=dict(cohort=cohort,pool=pool,n_held_out=holdout,n_cells_mature=count,shift=shift)
    for stream,g in [('all',d),*list(d.groupby('seed'))]:
        v=g[g.valid];n=len(v);r=keys|dict(seed=stream,n_attempted=len(g),n_valid=n,n_abstained=len(g)-n,abstention=1-n/len(g),median_width=float(v.width.median()),source_patients_min=int(g.n_source_patients.min()),source_patients_max=int(g.n_source_patients.max()),median_abs_target_gap=float(g.abs_target_gap.median()))
        for metric,flag in [('rejection','rejected'),('cell_coverage','cell_covered'),('patient_coverage','patient_covered'),('draw_inclusion','draw_covered')]:
            k=int(v[flag].sum());r['n_'+flag]=k;r[metric]=k/n if n else np.nan;r[metric+'_low'],r[metric+'_high']=wilson(k,n)
        r['unconditional_rejection']=int(v.rejected.sum())/len(g);out.append(r)
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source-dir',type=Path,required=True);parser.add_argument('--cache-dir',type=Path,required=True);args=parser.parse_args()
    sys.path.insert(0,str(args.source_dir.resolve()))
    old=json.loads((ROOT/'diagnostics/incoming_interval_provenance.json').read_text())
    source_inputs={}
    for path,sha in old['source_paths_sha256'].items():
        if path.startswith('src/'):
            actual=hashlib.sha256((args.source_dir/path).read_bytes()).hexdigest();assert actual==sha;source_inputs[path]=actual
    jobs=[(c,p,h,n,s,str(args.cache_dir/f'lee_{c}.npz')) for c in ['smc','kul3'] for p in ['pooled','reference'] for h in [2,5] for n in [50,800] for s in [.5,1.]];rows=[]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for f in as_completed([ex.submit(run,j) for j in jobs]):
            rows.extend(f.result());print('completed',len(rows)//6,'of 32',flush=True)
    d=pd.DataFrame(rows).sort_values(['cohort','pool','n_held_out','n_cells_mature','shift','seed'],key=lambda s:s.astype(str))
    summary=d[d.seed=='all'].drop(columns='seed');streams=d[d.seed!='all']
    original=pd.read_csv(ROOT/'diagnostics/incoming_interval_repair_pooled.csv');original=original[(original.candidate=='pseudobulk_patient_t')&original.n_cells_mature.isin([50,800])]
    keys=['cohort','pool','n_held_out','n_cells_mature','shift'];a=summary.set_index(keys).sort_index();b=original.set_index(keys).sort_index()
    for new,prior in [('n_valid','n_scored'),('n_rejected','n_excludes_zero'),('n_draw_covered','n_covered'),('n_abstained','n_abstained')]:np.testing.assert_array_equal(a[new],b[prior])
    print('PASS exact incoming counts on all 32 settings',flush=True)
    outputs={}
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    from datetime import date
    result_dir=ROOT.parents[1]/'results'/f'{date.today().isoformat()}_{revision[:7]}';result_dir.mkdir(exist_ok=True)
    for name,frame in [('patient_target_followup',summary),('patient_target_seed_rates',streams)]:
        path=ROOT/'diagnostics'/f'{name}.csv';frame.to_csv(path,index=False);outputs[path.name]=hashlib.sha256(path.read_bytes()).hexdigest();frame.to_parquet(result_dir/f'{name}.parquet',index=False)
    prov=dict(source_revision=revision,incoming_revision=old['source_revision'],settings=32,study_draws=64000,seed_streams=SEEDS,replicates_per_stream=400,original_counts_reproduced=True,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,source_input_sha256=source_inputs,cache_sha256={c:hashlib.sha256((args.cache_dir/f'lee_{c}.npz').read_bytes()).hexdigest() for c in ['smc','kul3']},input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'PATIENT_TARGET_DESIGN.md',ROOT/'diagnostics/incoming_interval_repair_pooled.csv']},output_sha256=outputs)
    (ROOT/'diagnostics/patient_target_provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
    for name in ['patient_target_followup','patient_target_seed_rates']:(result_dir/f'{name}.meta.json').write_text(json.dumps(prov,indent=2)+'\n')
if __name__=='__main__':main()
