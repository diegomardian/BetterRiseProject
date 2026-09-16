"""Generate the precision table and pooled rejection plot from saved aggregates."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent


def render_table(d):
    q=d[d.n_cells_mature==50].set_index(['cohort','pool','method','shift'])
    rows=[]
    def rate(row,metric):
        return f'{100*row[metric]:.1f} [{100*row[metric+"_low"]:.1f}, {100*row[metric+"_high"]:.1f}]'
    for cohort,pool in [('smc','pooled'),('smc','reference'),('kul3','pooled'),('kul3','reference')]:
        for method,label in [('percentile200','Percentile'),('welch','Welch')]:
            a=q.loc[(cohort,pool,method,.5)];n=q.loc[(cohort,pool,method,1.)]
            rows.append(f'{cohort.upper()} / {pool} & {label} & {rate(a,"coverage")} & {rate(n,"rejection")} & {rate(a,"rejection")}' + r' \\')
    return r'''\begin{table}[htbp]
\centering\small
\caption{At the 50-cell reporting boundary, changing intervals improves coverage
but does not restore the required performance. Rates (\%) and 95\% Wilson
Monte Carlo intervals use 5{,}000 studies per cohort, pool and effect. Methods
share the same samples. Coverage targets the fixed-pool effect at $s=0.5$;
null rejection uses $s=1$; detection uses $s=0.5$. Targets are at least 90\%
coverage, nominal 5\% null rejection and at least 80\% detection.}
\label{tab:precision}
\begin{tabular}{@{}llrrr@{}}
\toprule
Cohort / pool & Interval & Coverage & Null rejection & Detection \\
\midrule
'''+'\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
'''


def main():
    d=pd.read_csv(ROOT/'diagnostics/precision_followup.csv')
    (ROOT/'sections/precision_table.tex').write_text(render_table(d))
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(6.4,2.6),sharey=True)
    for ax,cohort in zip(axes,['smc','kul3']):
        for method,color,label in [('percentile200','#00758a','Percentile'),('welch','#a44800','Welch')]:
            for shift,style,kind in [(.5,'-','halving'),(1.,'--','null')]:
                q=d[(d.cohort==cohort)&(d.pool=='pooled')&(d.method==method)&(d['shift']==shift)].sort_values('n_cells_mature')
                ax.errorbar(q.n_cells_mature,100*q.rejection,yerr=[100*(q.rejection-q.rejection_low),100*(q.rejection_high-q.rejection)],color=color,linestyle=style,marker='o',markersize=3,capsize=2,label=f'{label}, {kind}')
        ax.set_title(cohort.upper()+' / pooled tissues');ax.set_xscale('log');ax.set_xticks([5,50,100,800],[5,50,100,800]);ax.set_xlabel('Mature cells in diseased arm');ax.set_ylim(0,90)
        ax.axhline(5,color='.4',linewidth=.8,alpha=.6)
    axes[0].set_ylabel('Exclusion of zero (%)')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,1.08))
    fig.tight_layout(rect=[0,0,1,.91])
    fig.savefig(ROOT/'figures/precision_rejection.pdf',bbox_inches='tight',metadata={'Creator':'','Producer':'','CreationDate':None})
    plt.close(fig)

if __name__=='__main__':main()
