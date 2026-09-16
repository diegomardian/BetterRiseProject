"""Regenerate manuscript figures from the pinned, unchanged result tables."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from _tables import pinned

OUT=Path(__file__).resolve().parent/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
                     'axes.titlesize':9,'legend.fontsize':7,'xtick.labelsize':7,
                     'ytick.labelsize':7,'pdf.fonttype':42})

def save(fig,name):
    for ax in fig.axes:
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.savefig(OUT/name,bbox_inches='tight',metadata={'Author':'','Creator':'Matplotlib'})
    plt.close(fig)

def controlled():
    raw=pd.read_parquet(pinned('controlled_grid_rates_r200_b200'))
    raw=raw[raw.binning=='per_count'].copy()
    keys=['cohort','pool','seed','n_cells_mature']
    # A shared count must be identical across grid views of the same draws.
    for col in ['coverage','discrimination','n_estimated']:
        assert raw.groupby(keys)[col].nunique(dropna=False).max()==1
    raw=raw.drop_duplicates(keys)
    raw=raw[(raw.n_cells_mature>0)&(raw.n_estimated>0)]
    fig,axs=plt.subplots(2,2,figsize=(5.5,3.55),sharex=True,sharey=True)
    colors=['#146b83','#b45309']
    for i,cohort in enumerate(['smc','kul3']):
        for j,pool in enumerate(['pooled','reference']):
            ax=axs[i,j];d=raw[(raw.cohort==cohort)&(raw.pool==pool)]
            for col,color,target in zip(['coverage','discrimination'],colors,[.90,.80]):
                s=d.groupby('n_cells_mature')[col].agg(['median','min','max']).sort_index()
                ax.plot(s.index,s['median'],color=color,lw=1.25,label=('Benchmark inclusion' if col=='coverage' else 'Exclusion of zero'))
                ax.fill_between(s.index,s['min'],s['max'],color=color,alpha=.13,lw=0)
                ax.axhline(target,color=color,ls='--',lw=.7,alpha=.8)
            ax.axvline(50,color='#777777',ls=':',lw=.9)
            ax.set_xscale('log');ax.set_ylim(0,1.03);ax.set_xlim(5,850)
            ax.set_xticks([10,50,100,400,800]);ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_title(f'{cohort.upper()} / '+('pooled tissues' if pool=='pooled' else 'reference tissue'))
            for side in ['top','right']:ax.spines[side].set_visible(False)
            if j==0:ax.set_ylabel('Probability')
            if i==1:ax.set_xlabel('Mature cells in diseased arm')
    axs[0,0].legend(loc='lower right',frameon=False,fontsize=6.5)
    fig.tight_layout(h_pad=1.0,w_pad=1.0)
    save(fig,'controlled_calibration.pdf')

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    controlled()
