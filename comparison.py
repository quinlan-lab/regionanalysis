import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import toolshed as ts
import seaborn as sns
sns.set_style('whitegrid')
from cyvcf2 import VCF

def read_values(path='/uufs/chpc.utah.edu/common/home/u1021864/analysis/scoredregions.bed'):
    gerp,phast,cadd,topgerp,topphast,topcadd,toptengerp,toptenphast,toptencadd,ccrpct,topccr=[],[],[],[],[],[],[],[],[],[],[]
    for i, region in enumerate(ts.reader(path, header="ordered")):
        g = np.mean(map(float,region['GERP'].split(",")))
        p = np.mean(map(float,region['phastCons'].split(",")))
        c = np.mean(map(float,region['CADD'].split(",")))
        pct=float(region['weighted_pct']); gene=region['gene']
        gerp.append(g);phast.append(p);cadd.append(c);ccrpct.append(pct)
        if pct>=90:
            toptengerp.append(g);toptenphast.append(p);toptencadd.append(c);topccr.append((gene,pct))
            if pct>=99:
                topgerp.append(g);topphast.append(p);topcadd.append(c)
    return gerp,phast,cadd,topgerp,topphast,topcadd,toptengerp,toptenphast,toptencadd,ccrpct,topccr

def barplot(dataset):
    mi,ma=min(dataset),max(dataset)
    rng = (mi,ma)
    t,t_edges=np.histogram(dataset, bins=40, range=rng)
    t=map(lambda x: float(x)/sum(t), t)
    width_t = (t_edges[1]-t_edges[0])
    return t, t_edges, width_t

def clinvar(path='patho.vcf'):
    clin=VCF(path)
    genes=[]
    for variant in clin:
        info = variant.INFO
        gene=info.get('GENEINFO')
        if gene:
            genes.append(gene.split(':')[0])
    return genes

gerp,phast,cadd,topgerp,topphast,topcadd,toptengerp,toptenphast,toptencadd,ccrpct,topccr=read_values('scoredregions.bed')
genes=clinvar('patho.vcf') #patho.vcf, vcf file after clinvarfilter.py for pathogenic variants
#TODO ADD top 1% and shrink point size for right size plots!!
fig, axarr = plt.subplots(3,2)
fig.tight_layout()
t,t_edges,width_t=barplot(gerp)
axarr[0,0].bar(t_edges[:-1], t, width = width_t, color = 'magenta', label = 'all', alpha = 0.5)
t,t_edges,width_t=barplot(topgerp)
axarr[0,0].bar(t_edges[:-1], t, width = width_t, color = 'cyan', label = 'top 1%', alpha = 0.5)
t,t_edges,width_t=barplot(toptengerp)
axarr[0,0].bar(t_edges[:-1], t, width = width_t, color = 'yellow', label = 'top 10%', alpha = 0.5)
axarr[0,0].set_xlabel('GERP');axarr[0,0].set_ylabel('Frequency')
axarr[0,0].legend(loc='best')
axarr[0,0].axvline(np.nanmean(gerp), color='magenta', linestyle='dashed', linewidth=2)
axarr[0,0].axvline(np.nanmean(topgerp), color='cyan', linestyle='dashed', linewidth=2)
axarr[0,0].axvline(np.nanmean(toptengerp), color='yellow', linestyle='dashed', linewidth=2)
t,t_edges,width_t=barplot(phast)
axarr[1,0].bar(t_edges[:-1], t, width = width_t, color = 'magenta', label = 'all', alpha = 0.5)
t,t_edges,width_t=barplot(topphast)
axarr[1,0].bar(t_edges[:-1], t, width = width_t, color = 'cyan', label = 'top 1%', alpha = 0.5)
t,t_edges,width_t=barplot(toptenphast)
axarr[1,0].bar(t_edges[:-1], t, width = width_t, color = 'yellow', label = 'top 10%', alpha = 0.5)
axarr[1,0].set_xlabel('phastCons');axarr[1,0].set_ylabel('Frequency')
axarr[1,0].legend(loc='best')
axarr[1,0].axvline(np.nanmean(phast), color='magenta', linestyle='dashed', linewidth=2)
axarr[1,0].axvline(np.nanmean(topphast), color='cyan', linestyle='dashed', linewidth=2)
axarr[1,0].axvline(np.nanmean(toptenphast), color='yellow', linestyle='dashed', linewidth=2)
t,t_edges,width_t=barplot(cadd)
axarr[2,0].bar(t_edges[:-1], t, width = width_t, color = 'magenta', label = 'all', alpha = 0.5)
t,t_edges,width_t=barplot(topcadd)
axarr[2,0].bar(t_edges[:-1], t, width = width_t, color = 'cyan', label = 'top 1%', alpha = 0.5)
t,t_edges,width_t=barplot(toptencadd)
axarr[2,0].bar(t_edges[:-1], t, width = width_t, color = 'yellow', label = 'top 10%', alpha = 0.5)
axarr[2,0].set_xlabel('CADD');axarr[2,0].set_ylabel('Frequency')
axarr[2,0].legend(loc='best')
axarr[2,0].axvline(np.nanmean(cadd), color='magenta', linestyle='dashed', linewidth=2)
axarr[2,0].axvline(np.nanmean(topcadd), color='cyan', linestyle='dashed', linewidth=2)
axarr[2,0].axvline(np.nanmean(toptencadd), color='yellow', linestyle='dashed', linewidth=2)
print np.nanmean(gerp), np.nanmean(phast), np.nanmean(cadd)
#color = ["white" if x[0] in genes else "red" for x in topccr]
ind = [i for i,x in enumerate(topccr) if x[0] not in genes]
topccr = [x[1] for x in topccr]
topccr = [topccr[i] for i in ind]; toptengerp=[toptengerp[i] for i in ind]; toptenphast=[toptenphast[i] for i in ind]; toptencadd=[toptencadd[i] for i in ind]
axarr[0,1].set_xlim(90,100)
axarr[1,1].set_xlim(90,100)
axarr[2,1].set_xlim(90,100)
axarr[0,1].scatter(topccr, toptengerp, c='magenta', marker='.', label='ccr_vs_gerp', alpha=0.2, s=1) #c=color
axarr[0,1].set_xlabel('CCR Percentile');axarr[0,1].set_ylabel('Score')
axarr[0,1].legend(loc='best')
axarr[0,1].axhline(np.nanmean(toptengerp), c='k', ls='dashed', lw=2)
axarr[1,1].scatter(topccr, toptenphast, c='lawngreen', marker='.', label='ccr_vs_phastcons', alpha=0.2, s=1)
axarr[1,1].set_xlabel('CCR Percentile');axarr[1,1].set_ylabel('Score')
axarr[1,1].legend(loc='best')
axarr[1,1].axhline(np.nanmean(toptenphast), c='k', ls='dashed', lw=2)
axarr[2,1].scatter(topccr, toptencadd, c='orange', marker='.', label='ccr_vs_cadd', alpha=0.2, s=1)
axarr[2,1].set_xlabel('CCR Percentile');axarr[2,1].set_ylabel('Score')
axarr[2,1].legend(loc='best')
axarr[2,1].axhline(np.nanmean(toptencadd), c='k', ls='dashed', lw=2)
plt.savefig('conservation.png',bbox_inches='tight')
