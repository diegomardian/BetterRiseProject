"""Render diagnostic rates and record reversals from pinned aggregate results."""
from pathlib import Path
import json
import pandas as pd
from _tables import pinned
ROOT=Path(__file__).resolve().parent

def render():
    d=pd.read_csv(ROOT/'diagnostics/low_count_diagnostic.csv').set_index(['cohort','pool','n_cells_mature','shift'])
    out=[r'\begin{table}[htbp]',r'\centering\small',r'\caption{Exclusion of zero (\%) with the same percentile interval under a halving',r'(effect, $s=0.5$) and the null ($s=1$). Each entry uses 200 finite draws from',r'one seed stream; nominal null rejection is 5\%. Diagnostic runs are exploratory.',r'Five-cell intervals are evaluated before the reporting gate, which abstains',r'below 20 cells.}',r'\label{tab:lowcounts}',r'\begin{tabular}{@{}lrrrrrrrr@{}}',r'\toprule',r'& \multicolumn{2}{c}{5 cells} & \multicolumn{2}{c}{50 cells} & \multicolumn{2}{c}{100 cells} & \multicolumn{2}{c}{800 cells} \\',r'cohort / pool & effect & null & effect & null & effect & null & effect & null \\',r'\midrule']
    for c,p in [('smc','pooled'),('smc','reference'),('kul3','pooled'),('kul3','reference')]:
        nums=[f'{100*d.loc[(c,p,n,s),"exclusion_rate"]:.1f}' for n in [5,50,100,800] for s in [.5,1.]]
        out.append(c.upper()+' / '+p+' & '+' & '.join(nums)+r' \\')
    return '\n'.join(out+[r'\bottomrule',r'\end{tabular}',r'\end{table}'])+'\n'

def reversals():
    d=pd.read_parquet(pinned('controlled_grid_rates_r200_b200'))
    d=d[d.binning=='per_count'].drop_duplicates(['cohort','pool','seed','n_cells_mature'])
    rows=[]
    for (c,p,seed),q in d.groupby(['cohort','pool','seed']):
        q=q.sort_values('n_cells_mature');v=(q.coverage>=.9)&(q.discrimination>=.8)
        first=float(q.loc[v,'n_cells_mature'].min()) if v.any() else None
        tail=[float(n) for n in q.n_cells_mature if ((q[q.n_cells_mature>=n].coverage>=.9)&(q[q.n_cells_mature>=n].discrimination>=.8)).all()]
        later=q[q.n_cells_mature>=first] if first is not None else q.iloc[:0]
        fails=later[(later.coverage<.9)|(later.discrimination<.8)]
        rows.append(dict(cohort=c,pool=p,seed=int(seed),first_candidate=first,later_failures=fails.n_cells_mature.tolist(),all_subsequent_evaluated_candidate=min(tail) if tail else None))
    return rows

if __name__=='__main__':
    (ROOT/'sections/low_count_table.tex').write_text(render())
    (ROOT/'diagnostics/candidate_reversals.json').write_text(json.dumps(reversals(),indent=2)+'\n')
