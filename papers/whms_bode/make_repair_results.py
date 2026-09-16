"""Render the complete fixed-pool interval comparison from saved summaries."""
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP
import pandas as pd
ROOT=Path(__file__).resolve().parent
METHODS=[('percentile200','Percentile, 200'),('percentile2000','Percentile, 2000'),('bca2000','BCa, 2000'),('bootstrap_t2000','Bootstrap-$t$, 2000'),('welch','Welch')]

def render_table(d):
    idx=d.set_index(['cohort','pool','n_cells_mature','method','shift']);rows=[]
    for cohort,pool in [('smc','pooled'),('smc','reference'),('kul3','pooled'),('kul3','reference')]:
        for j,(method,label) in enumerate(METHODS):
            vals=[]
            for n in (30,50):
                alt=idx.loc[(cohort,pool,n,method,.5)];null=idx.loc[(cohort,pool,n,method,1.)]
                vals += [alt.coverage,null.rejection,alt.rejection]
            name=f'{cohort.upper()} / {pool}' if j==0 else ''
            percentages = [(Decimal(str(v))*100).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP) for v in vals]
            rows.append(name+' & '+label+' & '+' & '.join(str(v) for v in percentages)+r' \\')
        rows.append(r'\addlinespace')
    return r'''\begin{table}[htbp]
\centering\small
\caption{Fixed-pool evaluation inside the flag-wide band and at the reporting
boundary. Each cohort/pool/count/shift uses 2{,}000 studies; all methods share
samples. Bootstrap budgets follow method names. $C_{.5}$ is coverage at the
halving, $R_1$ null rejection (one minus null coverage), and $D_{.5}$ detection
at the halving, all in percent. Nominal coverage is 95\%; the historical
alternative-coverage minimum is 90\% and detection requirement 80\%.
Maximum Monte Carlo SE is 1.12 percentage points; Wilson intervals and paired
comparisons are supplied with the aggregate results. All intervals are finite.}
\label{tab:repair}
\begin{tabular}{@{}llrrrrrr@{}}
\toprule
 & & \multicolumn{3}{c}{30 cells: flag-wide} & \multicolumn{3}{c}{50 cells: report} \\
\cmidrule(lr){3-5}\cmidrule(l){6-8}
Cohort / pool & Interval & $C_{.5}$ & $R_1$ & $D_{.5}$ & $C_{.5}$ & $R_1$ & $D_{.5}$ \\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
if __name__=='__main__':
    (ROOT/'sections/repair_table.tex').write_text(render_table(pd.read_csv(ROOT/'diagnostics/repair_followup.csv')))
