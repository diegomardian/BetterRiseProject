"""Write a compact, reproducible audit of pinned inputs and corrected counts."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import pandas as pd
from _tables import REPO_ROOT, MANIFEST

PAPER=Path(__file__).resolve().parent

def main():
    manifest=json.loads(MANIFEST.read_text())
    def read(k):return pd.read_parquet(REPO_ROOT/manifest[k])
    learned=read('trial_learned_generator_primary')
    draw=learned[learned.reference=='T_draw']
    a=draw[draw.n_patients==100].set_index(['generator','estimator'])
    b=draw[draw.n_patients==5000].set_index(['generator','estimator'])
    switches=(~a.allclose_residual_vs_zero)&(~a.allclose_estimate_vs_reference)&(~b.allclose_residual_vs_zero)&b.allclose_estimate_vs_reference
    rates=read('controlled_grid_rates_r200_b200')
    at50=rates[(rates.binning=='per_count')&(rates.grid=='dense')&(rates.n_cells_mature==50)&(rates.pool=='reference')]
    notes={}
    for cohort,g in at50.groupby('cohort'):
        notes[cohort]={'seeds':len(g),'coverage_range':[float(g.coverage.min()),float(g.coverage.max())],
          'discrimination_range':[float(g.discrimination.min()),float(g.discrimination.max())],
          'seeds_meeting_both':int(((g.coverage>=.9)&(g.discrimination>=.8)).sum())}
    excluded=int((draw.n_draws-draw.n_valid_pairs).sum())
    nr=int(draw.n_nonfinite_reference.sum());ne=int(draw.n_nonfinite_estimate.sum())
    out={
      'source_revision':'b503971',
      'external_control_source_revision':'a336594',
      'scope':'Saved-table verification plus paired fixed-pool coverage and null evaluation on the same 6,400 exploratory draws; source-pool variance summaries; adopted synthetic external-control study; fresh 160,000-study precision and paired interval comparison.',
      'single_cell_target':'f_N * (s - 1) * drawn_normal_mature_mean; draw-dependent reference, not fixed pool mean',
      'learned_generator':{
        'scheduled_estimator_evaluations':7*10*6*6*30,
        'scheduled_cohorts':7*6*6*30,
        'recorded_estimator_evaluations':int(draw.n_draws.sum()),
        'generated_cohorts':int(draw.groupby(['generator','n_patients']).n_draws.first().sum()),
        'valid_draw_reference_pairs':int(draw.n_valid_pairs.sum()),
        'excluded_pairs':excluded,'nonfinite_reference':nr,'nonfinite_estimate':ne,'overlap':nr+ne-excluded,
        'settings':len(draw),'minimum_valid_pairs_per_setting':int(draw.n_valid_pairs.min()),
        'settings_passing_only_relative_plus_absolute_comparison':int(((~draw.allclose_residual_vs_zero)&draw.allclose_estimate_vs_reference).sum()),
        'actual_tolerance_switches_100_to_5000':int(switches.sum()),
        'switching_generator_estimator_pairs':[list(x) for x in switches[switches].index],
      },
      'direct_reference_pool_performance_at_50':notes,
      'external_audit':{'repositories':len(read('prevalence_audit_verdicts')),'verdicts':read('prevalence_audit_verdicts').verdict.value_counts().to_dict()},
      'diagnostic_artifact_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((PAPER/'diagnostics').glob('*')) if p.is_file()},
      'precision_followup':json.loads((PAPER/'diagnostics/precision_provenance.json').read_text()),
      'additional_diagnostic':json.loads((PAPER/'diagnostics/low_count_provenance.json').read_text()),
      'input_sha256':{name:hashlib.sha256((REPO_ROOT/path).read_bytes()).hexdigest() for name,path in manifest.items()},
      'review_environment':{'python':platform.python_version(),**{name:importlib.metadata.version(name) for name in ['numpy','pandas','pyarrow','matplotlib','pypdf','pytest']}},
      'verification':{'dedicated_result_tests':29,'external_control_harness_tests':35,'historical_original_manuscript_result_tests':55,'learned_table_numeric_entries':30},
    }
    (PAPER/'results_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Wrote results_audit.json; hashed',len(manifest),'pinned input tables')

if __name__=='__main__':main()
