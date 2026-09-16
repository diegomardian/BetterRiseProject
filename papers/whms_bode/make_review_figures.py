"""Mechanism schematic and complete interval comparison with pointwise uncertainty."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parent
META={'Creator':'','Producer':'','CreationDate':None}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})


def mechanism():
    fig=plt.figure(figsize=(6.5,2.2));ax=fig.add_axes([.01,.04,.61,.93]);ax.axis('off')
    ax.text(.5,.97,'Same simulated draw; same reported interval',ha='center',va='top',weight='bold')
    ax.text(.5,.80,r'$e_T=\bar Y_T-s\mu_P$     $e_N=\bar Y_N-\mu_P$',ha='center')
    for x,title,formula,color in [(0,'Score the fixed parameter',r'$\hat\theta-\theta_P=f_N(e_T-e_N)$','#e5f1f3'),(.53,'Score the sampled reference',r'$\hat\theta-i_{\rm req}=f_N(e_T-s e_N)$','#f8eadb')]:
        ax.add_patch(FancyBboxPatch((x,.29),.46,.31,boxstyle='round,pad=0.015',facecolor=color,edgecolor='.6'))
        ax.text(x+.23,.52,title,ha='center',fontsize=8.2);ax.text(x+.23,.37,formula,ha='center',fontsize=9)
        ax.annotate('',xy=(x+.23,.63),xytext=(.5,.74),arrowprops={'arrowstyle':'->','color':'.4'})
    ax.text(.5,.10,'Tumor-arm error is unchanged.\nReference-arm error is multiplied by s.',ha='center',fontsize=9)
    b=fig.add_axes([.76,.20,.22,.58]);b.barh([1,0],[1,.25],height=.45,color=['#087c91','#b25a15']);b.set_xlim(0,1.15);b.set_yticks([1,0],['Fixed','Sampled']);b.set_xticks([0,.25,1],['0','25%','100%']);b.set_title('Reference-arm variance\nat a halving (s = 0.5)',fontsize=9);b.spines['left'].set_visible(False);b.spines['bottom'].set_visible(False);b.tick_params(length=0)
    fig.savefig(ROOT/'figures/reference_mechanism.pdf',bbox_inches='tight',metadata=META);plt.close(fig)


def repair():
    d=pd.read_csv(ROOT/'diagnostics/repair_followup.csv').set_index(['cohort','pool','n_cells_mature','method','shift'])
    methods=[('percentile200','Percentile 200'),('percentile2000','Percentile 2000'),('bca2000','BCa 2000'),('bootstrap_t2000','Bootstrap-t 2000'),('welch','Welch')]
    groups=[('smc','pooled'),('smc','reference'),('kul3','pooled'),('kul3','reference')]
    fig,axes=plt.subplots(1,3,figsize=(6.5,5.1),sharey=True)
    ticks=[];labels=[]
    for g,(c,p) in enumerate(groups):
        for j,(m,label) in enumerate(methods):
            y=g*6+j;ticks.append(y);labels.append(label)
            for n,off,color,marker in [(30,-.13,'#087c91','o'),(50,.13,'#b25a15','s')]:
                for ax,(shift,metric) in zip(axes,[(.5,'coverage'),(1.,'rejection'),(.5,'rejection')]):
                    r=d.loc[(c,p,n,m,shift)];v=100*r[metric]
                    ax.errorbar(v,y+off,xerr=[[v-100*r[metric+'_low']],[100*r[metric+'_high']-v]],fmt=marker,color=color,markersize=2.5,capsize=1,linewidth=.7)
        axes[0].text(-.99,g*6-1.05,c.upper()+' / '+p,transform=axes[0].get_yaxis_transform(),weight='bold',fontsize=8.5,clip_on=False)
    for ax,title,limits,target in zip(axes,['Alternative coverage','Null rejection','Detection'],[(65,100),(0,46),(30,85)],[95,5,80]):
        ax.set_title(title,fontsize=9);ax.set_xlim(*limits);ax.axvline(target,color='.35',linestyle='--',linewidth=.8);ax.set_xlabel('Percent');ax.set_ylim(23.2,-1.7);ax.grid(axis='x',alpha=.15)
    axes[0].axvline(90,color='.65',linestyle=':',linewidth=.8)
    axes[0].set_yticks(ticks,labels,fontsize=8)
    for n,col,mark in [(30,'#087c91','o'),(50,'#b25a15','s')]:axes[0].plot([],[],marker=mark,color=col,label=f'{n} cells',markersize=4,linestyle='none')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.61,1.03))
    fig.subplots_adjust(left=.31,right=.99,bottom=.1,top=.91,wspace=.23)
    fig.savefig(ROOT/'figures/repair_comparison.pdf',bbox_inches='tight',metadata=META);plt.close(fig)
if __name__=='__main__':mechanism();repair()
