"""Regenerate historical local figures (not compiled in the submission) from the pinned, unchanged result tables."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from _tables import pinned

OUT=Path(__file__).resolve().parent/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,
                     'axes.titlesize':10,'legend.fontsize':9,'xtick.labelsize':8,
                     'ytick.labelsize':8,'pdf.fonttype':42})

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
    fig,axs=plt.subplots(2,2,figsize=(5.8,5.0),sharex=True,sharey=True)
    colors=['#146b83','#b45309']
    for i,cohort in enumerate(['smc','kul3']):
        for j,pool in enumerate(['pooled','reference']):
            ax=axs[i,j];d=raw[(raw.cohort==cohort)&(raw.pool==pool)]
            for col,color,target in zip(['coverage','discrimination'],colors,[.90,.80]):
                s=d.groupby('n_cells_mature')[col].agg(['median','min','max']).sort_index()
                ax.plot(s.index,s['median'],color=color,lw=1.6,label=('Benchmark inclusion' if col=='coverage' else 'Exclusion of zero'))
                ax.fill_between(s.index,s['min'],s['max'],color=color,alpha=.13,lw=0)
                ax.axhline(target,color=color,ls='--',lw=1.1,alpha=.9)
            ax.axvline(50,color='#222222',ls=':',lw=1.3)
            ax.set_xscale('log');ax.set_ylim(0,1.03);ax.set_xlim(5,850)
            ax.set_xticks([5,20,50,200,800]);ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_title(f'{cohort.upper()} / '+('pooled tissues' if pool=='pooled' else 'reference tissue'))
            for side in ['top','right']:ax.spines[side].set_visible(False)
            if j==0:ax.set_ylabel('Probability')
            if i==1:ax.set_xlabel('Mature cells in diseased arm')
    h,l=axs[0,0].get_legend_handles_labels()
    fig.legend(h,l,loc='upper center',ncol=2,frameon=False,fontsize=9)
    for ax in axs.flat:
        ax.text(.97,.91,'90%',transform=ax.transAxes,ha='right',fontsize=7,color=colors[0])
        ax.text(.97,.73,'80%',transform=ax.transAxes,ha='right',fontsize=7,color=colors[1])
        ax.text(50,.04,'50',ha='center',fontsize=8,color='#222222')
    fig.tight_layout(rect=(0,0,1,.94),h_pad=1.4,w_pad=1.5)
    save(fig,'controlled_calibration.pdf')

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    controlled()
