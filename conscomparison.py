import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as patches
import numpy as np
import toolshed as ts
import seaborn as sns
sns.set_style('whitegrid')
from cyvcf2 import VCF
import cPickle as pickle
import os.path

def read_values(path='/uufs/chpc.utah.edu/common/home/u1021864/analysis/scoredregions.bed'):
    gerp,phast,cadd,topgerp,topphast,topcadd,toptengerp,toptenphast,toptencadd,ccrs=[],[],[],[],[],[],[],[],[],[]
    for i, region in enumerate(ts.reader(path, header="ordered")):
        g = np.mean(map(float,region['GERP'].split(",")))
        p = np.mean(map(float,region['phastCons'].split(",")))
        c = np.mean(map(float,region['CADD'].split(",")))
        pct=float(region['weighted_pct']); gene=region['gene']; chrom=region['chrom']; ranges=region['ranges']
        gerp.append(g);phast.append(p);cadd.append(c);ccrs.append((chrom,ranges,gene,pct))
        if pct>=90:
            toptengerp.append(g);toptenphast.append(p);toptencadd.append(c)
            if pct>=99:
                topgerp.append(g);topphast.append(p);topcadd.append(c)
    return gerp,phast,cadd,topgerp,topphast,topcadd,toptengerp,toptenphast,toptencadd,ccrs

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


if not os.path.isfile("var.pickle"):
    gerp,phast,cadd,topgerp,topphast,topcadd,toptengerp,toptenphast,toptencadd,ccrs=read_values('scoredregions.bed')
    #genes=clinvar('patho.vcf') #patho.vcf, vcf file after clinvarfilter.py for pathogenic variants
    var={}
    var['gerp']=gerp; var['phast']=phast; var['cadd']=cadd; var['toptengerp']=toptengerp; var['toptenphast']=toptenphast; var['toptencadd']=toptencadd; var['topgerp']=topgerp; var['topphast']=topphast; var['topcadd']=topcadd; var['ccrs']=ccrs
    pickle.dump(var, open("var.pickle", "wb"))

var = pickle.load(open("var.pickle", "rb"))
gerp=var['gerp']; phast=var['phast']; cadd=var['cadd']; toptengerp=var['toptengerp']; toptenphast=var['toptenphast']; toptencadd=var['toptencadd']; topgerp=var['topgerp']; tophphast=var['topphast']; topcadd=var['topcadd']; ccrs=var['ccrs']

fig, axarr = plt.subplots(2,3)
fig.tight_layout()
a,a_edges,width_a=barplot(gerp)
axarr[0,0].bar(a_edges[:-1], a, width = width_a, color = 'magenta', label = 'all', alpha = 0.5)
t,t_edges,width_t=barplot(toptengerp)
axarr[1,0].bar(t_edges[:-1], t, width = width_t, color = 'cyan', label = 'top 10%', alpha = 0.5)
#ma=max(max(t),max(a))
axarr[0,0].set_ylim(0,0.20)#axarr[0,0].set_ylim(0,ma)
axarr[1,0].set_ylim(0,0.20)#axarr[1,0].set_ylim(0,ma)
#t,t_edges,width_t=barplot(topgerp)
#axarr[2,0].bar(t_edges[:-1], t, width = width_t, color = 'yellow', label = 'top 1%', alpha = 0.5)
axarr[0,0].set_xlabel('Mean GERP');axarr[0,0].set_ylabel('Frequency')
axarr[1,0].set_xlabel('Mean GERP');axarr[1,0].set_ylabel('Frequency')
axarr[0,0].legend(loc='upper left')
axarr[1,0].legend(loc='upper left')
axarr[0,0].axvline(np.nanmean(gerp), color='magenta', linestyle='dashed', linewidth=2)
axarr[1,0].axvline(np.nanmean(gerp), color='magenta', linestyle='dashed', linewidth=2)

a,a_edges,width_a=barplot(phast)
axarr[0,1].bar(a_edges[:-1], a, width = width_a, color = 'magenta', label = 'all', alpha = 0.5)
t,t_edges,width_t=barplot(toptenphast)
axarr[1,1].bar(t_edges[:-1], t, width = width_t, color = 'cyan', label = 'top 10%', alpha = 0.5)
#ma=max(max(t),max(a))
axarr[0,1].set_ylim(0,0.35)#axarr[0,1].set_ylim(0,ma)
axarr[1,1].set_ylim(0,0.35)#axarr[1,1].set_ylim(0,ma)
#t,t_edges,width_t=barplot(topphast)
#axarr[2,1].bar(t_edges[:-1], t, width = width_t, color = 'yellow', label = 'top 1%', alpha = 0.5)
axarr[0,1].set_xlabel('Mean phastCons');axarr[0,1].set_ylabel('Frequency')
axarr[1,1].set_xlabel('Mean phastCons');axarr[1,1].set_ylabel('Frequency')
axarr[0,1].legend(loc='upper left')
axarr[1,1].legend(loc='upper left')
axarr[0,1].axvline(np.nanmean(phast), color='magenta', linestyle='dashed', linewidth=2)
axarr[1,1].axvline(np.nanmean(phast), color='magenta', linestyle='dashed', linewidth=2)

a,a_edges,width_a=barplot(cadd)
axarr[0,2].bar(a_edges[:-1], a, width = width_a, color = 'magenta', label = 'all', alpha = 0.5)
t,t_edges,width_t=barplot(toptencadd)
axarr[1,2].bar(t_edges[:-1], t, width = width_t, color = 'cyan', label = 'top 10%', alpha = 0.5)
#ma=max(max(t),max(a))
axarr[0,2].set_ylim(0,0.14)#axarr[0,2].set_ylim(0,ma)
axarr[1,2].set_ylim(0,0.14)#axarr[1,2].set_ylim(0,ma)
#t,t_edges,width_t=barplot(topcadd)
#axarr[2,2].bar(t_edges[:-1], t, width = width_t, color = 'yellow', label = 'top 1%', alpha = 0.5)
axarr[0,2].set_xlabel('Mean CADD');axarr[0,2].set_ylabel('Frequency')
axarr[1,2].set_xlabel('Mean CADD');axarr[1,2].set_ylabel('Frequency')
axarr[0,2].legend(loc='best')
axarr[1,2].legend(loc='best')
axarr[0,2].axvline(np.nanmean(cadd), color='magenta', linestyle='dashed', linewidth=2)
axarr[1,2].axvline(np.nanmean(cadd), color='magenta', linestyle='dashed', linewidth=2)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=None)
plt.savefig('consdist.pdf', format='pdf', bbox_inches='tight')
print np.nanmean(gerp), np.nanmean(phast), np.nanmean(cadd)


fig, axarr = plt.subplots(3)
fig.tight_layout()
#color = ["white" if x[0] in genes else "red" for x in topccr]
#color = ["red" if x < np.nanmean(ccrpct) else "blue" if x > 90 else "black" for x in ccrpct]
#ind = [i for i,x in enumerate(ccrpct) if x[0] not in genes]
ccrpct = [x[3] for x in ccrs]
#ccrpct = [ccrpct[i] for i in ind]; gerp=[gerp[i] for i in ind]; phast=[phast[i] for i in ind]; cadd=[cadd[i] for i in ind]

key=None;gerpsum=0;phastsum=0;caddsum=0
for i,j,k,l in zip(gerp,phast,cadd,ccrs):
    if key == None or key != l[0]+l[1]+l[2]:
        if i > np.nanmean(gerp) and l[3] < np.nanmean(ccrpct):
            gerpsum+=1
        if j > np.nanmean(phast) and l[3] < np.nanmean(ccrpct):
            phastsum+=1
        if k > np.nanmean(cadd) and l[3] < np.nanmean(ccrpct):
            caddsum+=1
    key=l[0]+l[1]+l[2]
print gerpsum,phastsum,caddsum
        

axarr[0].set_xlim(0,100)
axarr[1].set_xlim(0,100)
axarr[2].set_xlim(0,100)
g=axarr[0].hexbin(ccrpct, gerp, cmap='rainbow', bins='log')
axarr[0].set_xlabel('CCR Percentile');axarr[0].set_ylabel('GERP')
cbar = fig.colorbar(g, ax=axarr[0], orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
axarr[0].legend(loc='best')
axarr[0].add_patch(patches.Rectangle((90, min(gerp)), 10, np.nanmean(gerp)-min(gerp), fill=False, edgecolor='red', linewidth=2))
axarr[0].add_patch(patches.Rectangle((0, np.nanmean(gerp)), np.nanmean(ccrpct), max(gerp)-np.nanmean(gerp), fill=False, edgecolor='black', linewidth=2))

#TODO find # of regions in each box!!
p=axarr[1].hexbin(ccrpct, phast, cmap='rainbow', bins='log')
axarr[1].set_xlabel('CCR Percentile');axarr[1].set_ylabel('phastCons')
cbar = fig.colorbar(p, ax=axarr[1], orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
axarr[1].legend(loc='best')
axarr[1].add_patch(patches.Rectangle((90, min(phast)), 10, np.nanmean(phast)-min(phast), fill=False, edgecolor='red', linewidth=2))
axarr[1].add_patch(patches.Rectangle((0, np.nanmean(phast)), np.nanmean(ccrpct), max(phast)-np.nanmean(phast), fill=False, edgecolor='black', linewidth=2))

c=axarr[2].hexbin(ccrpct, cadd, cmap='rainbow', bins='log')
axarr[2].set_xlabel('CCR Percentile');axarr[2].set_ylabel('CADD')
cbar = fig.colorbar(c, ax=axarr[2], orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
axarr[2].legend(loc='best')
axarr[2].add_patch(patches.Rectangle((90, min(cadd)), 10, np.nanmean(cadd)-min(cadd), fill=False, edgecolor='red', linewidth=2))
axarr[2].add_patch(patches.Rectangle((0, np.nanmean(cadd)), np.nanmean(ccrpct), max(cadd)-np.nanmean(cadd), fill=False, edgecolor='black', linewidth=2))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)
plt.savefig('consvsconst.pdf', format='pdf', bbox_inches='tight')
