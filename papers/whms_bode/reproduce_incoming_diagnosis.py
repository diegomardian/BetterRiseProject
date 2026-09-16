"""Reproduce pinned diagnosis rates using an extracted incoming source tree.

Extract ae8c109's src/ into a temporary directory, then run this script with
--source-dir pointing there, --raw-dir pointing to local Lee data, and
--work-dir pointing to a scratch directory. Only aggregates belong in the repo.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,sys,time
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-dir',required=True,type=Path)
    parser.add_argument('--raw-dir',required=True,type=Path)
    parser.add_argument('--work-dir',required=True,type=Path)
    args=parser.parse_args();sys.path.insert(0,str(args.source_dir.resolve()))
    from src.harness.calibration_gap import load_cohort_arrays
    from src.harness.run_interval_repair import diagnosis_cell,_jobs
    start=time.monotonic();args.work_dir.mkdir(parents=True,exist_ok=True);cache={}
    for cohort in ['smc','kul3']:
        a=load_cohort_arrays(cohort,raw_dir=args.raw_dir.resolve())
        path=args.work_dir/f'lee_{cohort}.npz'
        np.savez_compressed(path,counts=a['counts'],cell_type=np.asarray(a['cell_type']),patient_id=np.asarray(a['patient_id']),tissue=np.asarray(a['tissue']),genes=np.asarray(a['genes']),study_id=np.asarray([a['study_id']]),n_patients=np.asarray([a['n_patients']]))
        cache[cohort]=path
    jobs=_jobs(cache,replicates=400,n_boot=2000,task='diagnosis');rows=[]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for f in as_completed([ex.submit(diagnosis_cell,j) for j in jobs]):
            rows.extend(f.result());print('completed',len(rows)//8,'of',len(jobs),'jobs',flush=True)
    d=pd.DataFrame(rows);d.to_csv(args.work_dir/'diagnosis_reproduced_streams.csv',index=False)
    keys=['cohort','pool','n_cells_mature','family','design']
    d['rejects']=(d.null_rejection*d.n_scored).round().astype(int)
    a=d.groupby(keys).agg(n_scored=('n_scored','sum'),rejects=('rejects','sum'))
    a['null_rejection']=a.rejects/a.n_scored
    b=pd.read_csv(ROOT/'diagnostics/incoming_interval_diagnosis_rates.csv').set_index(keys).sort_index()
    assert a.index.equals(b.index)
    np.testing.assert_array_equal(a.n_scored,b.n_scored)
    np.testing.assert_allclose(a.null_rejection,b.null_rejection,atol=1e-15,rtol=0)
    print('PASS all 64 diagnosis rates reproduced from raw source data',round(time.monotonic()-start,1),'seconds')
if __name__=='__main__':main()
