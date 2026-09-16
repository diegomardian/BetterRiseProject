"""Paired target comparison, Monte Carlo precision and source-pool summaries."""
from pathlib import Path
import json
import math
import pandas as pd
ROOT=Path(__file__).resolve().parent

def wilson(k,n,z=1.959963984540054):
    p=k/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den
    radius=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-radius,center+radius

def signal_to_noise(mu,variance,n_t,n_n=800,s=.5):
    # Independent replacement draws followed by binomial thinning of the tumour arm.
    return (1-s)*mu/math.sqrt((s*s*variance+s*(1-s)*mu)/n_t+variance/n_n)

def fixed_table(d):
    q=d[d['shift']==.5].set_index(['cohort','pool','n_cells_mature'])
    out=[r'\begin{table}[htbp]',r'\centering\small',r'\caption{Paired evaluation of the same nominal 95\% intervals at $s=0.5$.',r'Entries are drawn-reference inclusion (D) and fixed-pool coverage (F), in \%.',r'Each entry uses 200 finite replicates; only the evaluation target changes.}',r'\label{tab:fixedpool}',r'\begin{tabular}{@{}lrrrrrrrr@{}}',r'\toprule',r'& \multicolumn{2}{c}{5 cells} & \multicolumn{2}{c}{50 cells} & \multicolumn{2}{c}{100 cells} & \multicolumn{2}{c}{800 cells} \\',r'cohort / pool & D & F & D & F & D & F & D & F \\',r'\midrule']
    for c,p in [('smc','pooled'),('smc','reference'),('kul3','pooled'),('kul3','reference')]:
        nums=[f'{100*q.loc[(c,p,n),col]:.1f}' for n in [5,50,100,800] for col in ['coverage','fixed_pool_coverage']]
        out.append(c.upper()+' / '+p+' & '+' & '.join(nums)+r' \\')
    return '\n'.join(out+[r'\bottomrule',r'\end{tabular}',r'\end{table}'])+'\n'

def reversal_table():
    rows=json.loads((ROOT/'diagnostics/candidate_reversals.json').read_text())
    out=[r'\begin{table}[htbp]',r'\centering\small',r'\caption{The three SMC/reference seed streams with a later failing count.',r'The last column requires all subsequent evaluated counts through 800 to pass;',r'it is a descriptive alternative to the first-crossing rule.}',r'\label{tab:reversals}',r'\begin{tabular}{@{}lrrr@{}}',r'\toprule',r'seed stream & first qualifying & later failing & all subsequent pass \\',r'\midrule']
    for r in rows:
        if r['later_failures']:
            out.append(f"{r['seed']} & {r['first_candidate']:.0f} & "+', '.join(f'{v:.0f}' for v in r['later_failures'])+f" & {r['all_subsequent_evaluated_candidate']:.0f}"+r' \\')
    return '\n'.join(out+[r'\bottomrule',r'\end{tabular}',r'\end{table}'])+'\n'

def main():
    d=pd.read_csv(ROOT/'diagnostics/low_count_diagnostic.csv')
    (ROOT/'sections/fixed_pool_table.tex').write_text(fixed_table(d))
    (ROOT/'sections/reversal_table.tex').write_text(reversal_table())
    precision=[]
    for r in d.itertuples():
        precision.append(dict(cohort=r.cohort,pool=r.pool,count=r.n_cells_mature,shift=r.shift,fixed_pool_wilson95=wilson(r.n_fixed_pool_covered,r.n_valid),draw_benchmark_wilson95=wilson(r.n_draw_benchmark_covered,r.n_valid),null_or_effect_exclusion_wilson95=wilson(r.n_excludes_zero,r.n_valid),paired_difference=r.paired_coverage_difference,paired_difference_mc_se=r.paired_difference_mc_se))
    (ROOT/'diagnostics/coverage_precision.json').write_text(json.dumps(precision,indent=2)+'\n')
    p=pd.read_csv(ROOT/'diagnostics/source_pool_properties.csv')
    p['snr_at_50']=[signal_to_noise(r.mean_expression,r.variance_expression,50) for r in p.itertuples()]
    summaries=[]
    for (cohort,pool),g in p.groupby(['cohort','pool']):
        n_patients=10 if cohort=='smc' else 6
        summaries.append(dict(cohort=cohort,pool=pool,pairs=len(g),total_source_cells=int(g.n_source_cells.sum()/(n_patients-1)),median_source_cells=float(g.n_source_cells.median()),median_mean=float(g.mean_expression.median()),median_squared_cv=float(g.squared_cv.median()),median_zero_fraction=float(g.zero_fraction.median()),median_snr_at_50=float(g.snr_at_50.median())))
    (ROOT/'diagnostics/source_pool_summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    print(json.dumps([s for s in summaries if s['pool']=='reference'],indent=2))
    print('50-cell fixed-pool Wilson intervals:')
    for x in precision:
        if x['count']==50 and x['shift']==.5: print(x['cohort'],x['pool'],[round(v*100,1) for v in x['fixed_pool_wilson95']])

if __name__=='__main__':main()
