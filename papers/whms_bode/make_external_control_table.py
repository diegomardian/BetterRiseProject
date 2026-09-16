"""Render the target-population comparison from version-pinned synthetic results."""
from pathlib import Path
import hashlib
import json
import pandas as pd
from _tables import pinned, REPO_ROOT

ROOT = Path(__file__).resolve().parent


def render(primary):
    selected = primary[primary.n_trial == 1600].set_index(['reference', 'estimator'])
    rows = []
    for ref, label in [('obs-pooled', 'Observed contrast, pooled'),
                       ('po-pooled', 'Potential outcomes, pooled'),
                       ('po-trial', 'Potential outcomes, enrolled')]:
        values = []
        for estimator in ['ate-standardisation', 'att-standardisation']:
            row = selected.loc[(ref, estimator)]
            values.extend([f'{row.bias:+.3f}' if row.bias else '0.000', f'{row.rmse:.3f}'])
        rows.append(label + ' & ' + ' & '.join(values) + r' \\')
    return r'''\begin{table}[htbp]
\centering\small
\caption{Reference choice reverses the estimator ranking in the synthetic
external-control example. Bias and RMSE use 200 studies with 1{,}600 enrolled
records and 4{,}800 controls. Columns identify the population used to average
observed stratum contrasts; only enrolled standardization matches the intended
target. All entries are in abstract outcome units.}
\label{tab:eca}
\begin{tabular}{@{}lrrrr@{}}
\toprule
& \multicolumn{2}{c}{Pooled standardization} & \multicolumn{2}{c}{Enrolled standardization} \\
\cmidrule(lr){2-3}\cmidrule(l){4-5}
Validation reference & Bias & RMSE & Bias & RMSE \\
\midrule
''' + '\n'.join(rows) + r'''
\bottomrule
\end{tabular}
\end{table}
'''


def main():
    inputs = {k: pinned(f'external_control_{k}') for k in ['primary','shift','balance','falsifiers']}
    primary = pd.read_parquet(inputs['primary'])
    (ROOT/'sections/external_control_table.tex').write_text(render(primary))
    # Tables contain only synthetic, aggregate results. Sidecars below retain
    # exact local provenance; the anonymous archive gets the CSVs only.
    hashes = {}
    for name, path in inputs.items():
        pd.read_parquet(path).to_csv(ROOT/f'diagnostics/external_control_{name}.csv', index=False)
        for source in [path, path.with_suffix('.meta.json')]:
            hashes[str(source.relative_to(REPO_ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
    for rel in ['src/harness/external_control_demo.py','docs/prereg_external_control_demo.md']:
        hashes[rel] = hashlib.sha256((REPO_ROOT/rel).read_bytes()).hexdigest()
    (ROOT/'diagnostics/external_control_provenance.json').write_text(json.dumps({
        'source_revision':'a336594', 'experiment_revision':'6b4a31a',
        'synthetic':True, 'seed':20260915, 'seed_streams':1,
        'primary_study_draws':1200, 'comparisons_per_study':40,
        'target':'Mean individual treatment effect over the enrolled records in each draw',
        'input_sha256':hashes,
    },indent=2)+'\n')
    print('Rendered external-control table and exported aggregate evidence.')


if __name__ == '__main__':
    main()
