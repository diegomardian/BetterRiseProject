"""Export and audit pinned incoming results without merging either branch."""
from pathlib import Path
from io import BytesIO
import hashlib,json,subprocess
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
REV='ae8c109'
RESULTS='results/2026-09-15_effd3b2'
NAMES=['interval_diagnosis_rates','interval_diagnosis_attribution','interval_repair_rates','interval_repair_pooled','interval_repair_by_seed']


def blob(path):
    return subprocess.check_output(['git','show',f'{REV}:{path}'],cwd=ROOT)


def check_aggregates(tables):
    d=tables['interval_diagnosis_rates'];r=tables['interval_repair_rates']
    assert len(d)==64 and (d.n_scored==2000).all() and (d.n_pool_undefined==0).all()
    np.testing.assert_allclose(d.null_rejection_mcse,np.sqrt(d.null_rejection*(1-d.null_rejection)/d.n_scored))
    q=d.pivot(index=['cohort','pool','n_cells_mature','family'],columns='design',values='null_rejection')
    np.testing.assert_array_equal(q.xs(800,level='n_cells_mature').balanced,q.xs(800,level='n_cells_mature').fixed_fraction)
    keys=['cohort','pool','n_held_out','n_cells_mature','candidate','shift']
    by=tables['interval_repair_by_seed'];po=tables['interval_repair_pooled'].set_index(keys).sort_index()
    counts=['n_attempted','n_scored','n_abstained','n_excludes_zero','n_covered']
    summed=by.groupby(keys)[counts].sum().sort_index()
    np.testing.assert_array_equal(summed,po[counts])
    assert (summed.n_attempted==summed.n_scored+summed.n_abstained).all()
    np.testing.assert_allclose(po.excludes_zero_rate,po.n_excludes_zero/po.n_scored,equal_nan=True)
    np.testing.assert_allclose(po.coverage,po.n_covered/po.n_scored,equal_nan=True)
    null=po.xs(1.,level='shift');alt=po.xs(.5,level='shift')
    np.testing.assert_array_equal(null.n_covered+null.n_excludes_zero,null.n_scored)
    idx=r.set_index(keys[:-1]).sort_index()
    np.testing.assert_allclose(idx.null_rejection,null.excludes_zero_rate,equal_nan=True)
    np.testing.assert_allclose(idx.discrimination,alt.excludes_zero_rate,equal_nan=True)
    assert (r.verdict=='pass').sum()==11
    assert (r.null_rejection<=.05+2*r.null_rejection_mcse).sum()==46
    assert set(r[r.n_held_out==5].candidate)=={'percentile','patient_cluster','pseudobulk_patient_t'}
    return {'diagnosis_rows':len(d),'repair_rows':len(r),'repair_seed_rows':len(by),'source_joint_pass_labels':11,'source_null_tolerance_passes':46,'repair_count_reaggregation':'passed','null_coverage_complementarity':'passed'}


def main():
    inputs={};outputs={};tables={}
    for n in NAMES:
        raw=blob(f'{RESULTS}/{n}.parquet');tables[n]=pd.read_parquet(BytesIO(raw))
        for ext in ['parquet','meta.json']:
            path=f'{RESULTS}/{n}.{ext}';inputs[path]=hashlib.sha256(blob(path)).hexdigest()
        target=ROOT/'diagnostics'/f'incoming_{n}.csv'
        exported=tables[n]
        if n=='interval_repair_rates':
            exported=exported.rename(columns={'coverage':'draw_reference_inclusion','coverage_mcse':'draw_reference_inclusion_mcse','verdict':'source_verdict'})
        exported.to_csv(target,index=False);outputs[target.name]=hashlib.sha256(target.read_bytes()).hexdigest()
    report=check_aggregates(tables)
    for path in ['src/harness/interval_diagnosis.py','src/harness/interval_repair.py','src/harness/run_interval_repair.py','src/harness/pseudobulk.py','docs/prereg_repaired_interval.md']:
        inputs[path]=hashlib.sha256(blob(path)).hexdigest()
    path=ROOT/'diagnostics/incoming_interval_provenance.json';prov=json.loads(path.read_text())
    prov.update(source_paths_sha256=inputs,aggregate_exports_sha256=outputs,aggregate_audit=report)
    path.write_text(json.dumps(prov,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
