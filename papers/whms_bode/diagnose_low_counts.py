"""Exploratory bootstrap diagnostic; aggregate outputs only, no cell-level export.

Replays selected settings of the first controlled-grid seed and adds its null.
Preserves the original generator, interval, sample size and configuration seeds.
Not a new calibration, threshold search or patient-population validation.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from src.harness.calibration_gap import load_cohort_arrays, _pool_mask
from src.harness.attenuation import _composition, _configuration_seed
from src.harness.pseudobulk import patient_holdout, generate_pseudobulk
from src.harness.interval import within_patient_intrinsic_ci

COUNTS=(5,50,100,800)
SEED=20260831
REPLICATES=200

def main():
    summaries=[]
    for cohort in ('smc','kul3'):
        a=load_cohort_arrays(cohort)
        for pool in ('pooled','reference'):
            mask=_pool_mask(pd.Series(a['tissue']),pool)
            counts=np.asarray(a['counts'])[mask]
            celltype=np.asarray(a['cell_type'])[mask]
            patient=np.asarray(a['patient_id'])[mask]
            genes=a['genes'];target='GUCA2A';mature='differentiated'
            types=sorted(set(celltype));cn=_composition(.4,types,mature)
            for n in COUNTS:
                ct=_composition(n/2000,types,mature)
                for shift in (.5,1.):
                    rows=[]
                    for rep in range(REPLICATES):
                        seed=_configuration_seed(SEED,n_cells=2000,mature_fraction=n/2000,shift=shift,replicate=rep)
                        _,held=patient_holdout(patient,n_held_out=2,seed=seed)
                        s=generate_pseudobulk(counts,celltype,patient,genes,composition_normal=cn,composition_tumour=ct,shift={target:shift},held_out_patients=held,n_cells=2000,seed=seed,mature_label=mature)
                        cells=s.mature_expression[target];normal=cells['normal'];tumour=cells['tumour']
                        lo,hi=within_patient_intrinsic_ci(normal,tumour,frac_mature_normal=.4,frac_mature_tumour=n/2000,n_boot=200,seed=seed,weighting='normal')
                        truth=s.truth.parametric[target]['normal']['intrinsic']
                        assert lo is not None and hi is not None
                        rows.append(dict(width=hi-lo,covered=lo<=truth<=hi,exclude_zero=(lo>0 or hi<0),tumour_all_zero=bool((tumour==0).all()),normal_all_zero=bool((normal==0).all()),zero_width=hi==lo))
                    d=pd.DataFrame(rows);z=d.tumour_all_zero
                    result=dict(cohort=cohort,pool=pool,n_cells_mature=n,shift=shift,seed=SEED,n_attempted=len(d),n_valid=len(d),n_excludes_zero=int(d.exclude_zero.sum()),exclusion_rate=d.exclude_zero.mean(),coverage=d.covered.mean(),median_ci_width=d.width.median(),n_zero_width=int(d.zero_width.sum()),n_tumour_all_zero=int(z.sum()),n_excludes_zero_and_tumour_all_zero=int((d.exclude_zero&z).sum()),n_normal_all_zero=int(d.normal_all_zero.sum()))
                    summaries.append(result)
                    print(cohort,pool,n,shift,'excludes',round(result['exclusion_rate'],3),'all-zero T',result['n_tumour_all_zero'],flush=True)
    out=ROOT/'diagnostics';out.mkdir(exist_ok=True)
    result=pd.DataFrame(summaries)
    manifest=json.loads((ROOT/'results_manifest.json').read_text())
    original=pd.read_parquet(REPO/manifest['controlled_grid_rates_r200_b200'])
    original=original[(original.binning=='per_count')&(original.seed==SEED)].drop_duplicates(['cohort','pool','n_cells_mature'])
    both=result[result['shift']==.5].merge(original,on=['cohort','pool','n_cells_mature'],suffixes=('_diagnostic','_original'))
    assert len(both)==16
    for lhs,rhs in [('coverage_diagnostic','coverage_original'),('exclusion_rate','discrimination'),('median_ci_width_diagnostic','median_ci_width_original')]:
        np.testing.assert_allclose(both[lhs],both[rhs],rtol=1e-12,atol=1e-12)
    result.to_csv(out/'low_count_diagnostic.csv',index=False)
    provenance={'purpose':'Exploratory response to low-count anomaly; not independent calibration','seed':SEED,'counts':COUNTS,'replicates':REPLICATES,'inner_bootstrap':200,'original_alternative_settings_reproduced':16,'input_source_files':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [REPO/'src/harness/pseudobulk.py',REPO/'src/harness/interval.py',REPO/'src/harness/attenuation.py',Path(__file__)]},'output_sha256':hashlib.sha256((out/'low_count_diagnostic.csv').read_bytes()).hexdigest()}
    (out/'low_count_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('PASS: all 16 alternative settings reproduce the original first-seed summaries',flush=True)

if __name__=='__main__':main()
