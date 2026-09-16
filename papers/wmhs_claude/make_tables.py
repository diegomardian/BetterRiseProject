"""Generate the complete learned-estimator table; --check rejects stale values."""
from pathlib import Path
import argparse
import pandas as pd
from _tables import pinned

LABELS = {
    'gcomp-saturated':'G-computation, saturated',
    'ipw-saturated':'IPW, saturated',
    'gcomp-gmm-outcome-K1':'G-comp., mixture outcome $K=1$',
    'gcomp-gmm-outcome-K4':'G-comp., mixture outcome $K=4$',
    'gcomp-ridge-outcome-a1e-5':r'G-comp., ridge $\alpha=10^{-5}$',
    'gcomp-ridge-outcome-a1e-3':r'G-comp., ridge $\alpha=10^{-3}$',
    'gcomp-ridge-outcome-a1e-1':r'G-comp., ridge $\alpha=10^{-1}$',
    'gcomp-mlp-outcome':'G-comp., MLP outcome',
    'ols-stratum-dummies':'OLS, stratum dummies',
    'unadjusted':'Unadjusted difference',
}

def number(v):
    if v == 0: return '$0$'
    if abs(v) < .001:
        s,e=f'{v:.1e}'.split('e')
        return '$'+s+r'\!\times\!10^{'+str(int(e))+'}$'
    return f'${v:.3g}$'

def render():
    d=pd.read_parquet(pinned('trial_learned_generator_primary'))
    d=d[d.reference=='T_draw'].copy()
    assert set(d.estimator)==set(LABELS)
    d['family']=d.generator.map(lambda x:'Mixture' if x.startswith('gmm') else 'MLP' if x.startswith('mlp') else 'Parametric')
    agg=d.groupby(['estimator','family']).max_residual.max()
    rows=[label+' & '+' & '.join(number(agg.loc[key,fam]) for fam in ['Mixture','MLP','Parametric'])+r' \\' for key,label in LABELS.items()]
    return r'''\begin{table}[htbp]
\small\centering
\caption{Maximum absolute residual against the drawn reference, over cohort
sizes and generators within each family, conditional on finite pairs. All ten
estimators in the expanded run are shown. The saturated zero row is exact;
other entries are rounded from the pinned primary table. Changing generator
capacity preserves the saturated identities but can change other residuals.}
\label{tab:learned}
\begin{tabular}{@{}lccc@{}}
\toprule
Estimator & Mixtures & MLP generator & Parametric controls \\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
'''

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
    target=Path(__file__).resolve().parent/'sections/learned_table.tex'
    expected=render()
    if args.check:
        assert target.read_text()==expected,'Learned-generator table is stale; rerun make_tables.py'
        print('PASS: all 30 learned-table entries match the pinned results')
    else: target.write_text(expected)
