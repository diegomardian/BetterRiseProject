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
                ax.plot(s.index,s['median'],color=color,lw=1.25,label=col.capitalize())
                ax.fill_between(s.index,s['min'],s['max'],color=color,alpha=.13,lw=0)
                ax.axhline(target,color=color,ls='--',lw=.7,alpha=.8)
            ax.axvline(50,color='#777777',ls=':',lw=.9)
            ax.set_xscale('log');ax.set_ylim(0,1.03);ax.set_xlim(5,850)
            ax.set_xticks([10,50,100,400,800]);ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_title(f'{cohort.upper()} / '+('pooled tissues' if pool=='pooled' else 'reference tissue'))
            for side in ['top','right']:ax.spines[side].set_visible(False)
            if j==0:ax.set_ylabel('Probability')
            if i==1:ax.set_xlabel('Mature cells in diseased arm')
    axs[0,0].legend(loc='lower right',frameon=False)
    fig.tight_layout(h_pad=1.0,w_pad=1.0)
    save(fig,'controlled_calibration.pdf')

def recovery():
    d=pd.read_parquet(pinned('calibration_gap_recovery'))
    d=d[(d['shift']==.5)&(d.grid=='extended')]
    fig,ax=plt.subplots(figsize=(5.5,1.8))
    for pool,color,label in [('pooled','#146b83','Pooled tissues'),('reference','#b45309','Reference tissue')]:
        s=d[d.pool==pool].groupby('median_n_cells_mature')[['ratio_median','ratio_q25','ratio_q75']].median().sort_index()
        s=s[s.index>0]
        ax.plot(s.index,s.ratio_median,'o-',ms=2.5,lw=1.2,color=color,label=label)
        ax.fill_between(s.index,s.ratio_q25,s.ratio_q75,color=color,alpha=.14,lw=0)
    ax.axhline(1,color='#555555',ls='--',lw=.8)
    ax.axvline(50,color='#777777',ls=':',lw=.8)
    ax.set_xscale('log');ax.set_xlim(15,850);ax.set_ylim(.4,1.9)
    ax.set_xticks([20,50,100,200,400,800]);ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel('Mature cells in diseased arm');ax.set_ylabel('Estimate / requested reference')
    ax.legend(frameon=False,loc='upper right');ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();save(fig,'recovery.pdf')

def trial():
    d=pd.read_parquet(pinned('trial_recovery'))
    styles=[('gcomp-from-generator','#146b83','o','Saturated G-comp.'),('ipw-saturated','#56a7ba','s','Saturated IPW'),('ipw-cross-fitted','#2e7d32','^','Cross-fitted IPW'),('ols-stratum-dummies','#7955a0','D','OLS'),('unadjusted','#b45309','v','Unadjusted')]
    fig,axs=plt.subplots(1,2,figsize=(5.5,2.6))
    for name,color,marker,label in styles:
        r=d[d.estimator==name].sort_values('n_patients');x=r.n_patients
        axs[0].plot(x,r.ratio_median,color=color,marker=marker,ms=2.5,lw=1,label=label)
        axs[0].fill_between(x,r.ratio_q25,r.ratio_q75,color=color,alpha=.10,lw=0)
        axs[1].plot(x,r.max_residual_vs_realised.clip(lower=1e-16),color=color,marker=marker,ms=2.5,lw=1)
    axs[0].axhline(1,color='#555555',ls=':',lw=.8);axs[0].set_ylim(.3,2.8)
    axs[0].set_title('Requested-effect recovery');axs[0].set_ylabel('Estimated / requested effect')
    axs[1].set_title('Empirical-reference residual');axs[1].set_yscale('log');axs[1].set_ylim(2e-17,100)
    axs[1].set_ylabel(r'Maximum $|\hat\theta-T(D)|$');axs[1].axhspan(2e-17,1e-12,color='#ededed',zorder=0)
    axs[1].text(120,2e-12,'Numerical equality',fontsize=6.5,color='#555555')
    for ax in axs:
        ax.set_xscale('log');ax.set_xlabel('Patients per cohort');ax.spines[['top','right']].set_visible(False)
    h,l=axs[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,-.01))
    fig.tight_layout(rect=(0,.17,1,1));save(fig,'trial_recovery.pdf')

def historical():
    d=pd.read_parquet(pinned('calibration_gap_bins_r500'))
    d=d[d.grid=='extended']
    fig,axs=plt.subplots(1,2,figsize=(5.5,2.3),sharey=True)
    for ax,pool in zip(axs,['pooled','reference']):
        q=d[d.pool==pool]
        ax.axvspan(1,20,color='#eeeeee');ax.axvspan(20,50,color='#f6f6f6')
        for col,color,target in [('coverage','#146b83',.9),('discrimination','#b45309',.8)]:
            r=q.groupby('n_cells_mature')[col].agg(['median','min','max']).sort_index()
            ax.plot(r.index,r['median'],color=color,marker='o',ms=2.5,lw=1,label=col.capitalize())
            ax.fill_between(r.index,r['min'],r['max'],color=color,alpha=.13,lw=0)
            ax.axhline(target,color=color,ls='--',lw=.7)
        ax.axvline(50,color='#777777',ls=':',lw=.8)
        ax.set_xscale('log');ax.set_xlim(8,850);ax.set_ylim(0,1.05)
        ax.set_xticks([10,50,100,400,800]);ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xlabel('Median cell count within bin');ax.spines[['top','right']].set_visible(False)
    axs[0].set_title('Pooled tissues: no qualifying bin');axs[0].set_ylabel('Probability')
    axs[1].set_title('Reference tissue: candidate 90 (7/8 seeds)')
    axs[1].axvline(90,color='#333333',ls=':',lw=.9)
    h,l=axs[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=2,frameon=False)
    fig.tight_layout(rect=(0,.10,1,1));save(fig,'historical_calibration.pdf')

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    controlled();trial()
