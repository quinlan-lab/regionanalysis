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
from collections import defaultdict
#import cProfile

def read_values(path='/uufs/chpc.utah.edu/common/home/u1021864/analysis/scoredregions.bed'):
    var=defaultdict(defaultdict)
    ccrs=defaultdict(list); genes=defaultdict(list)
    for i, region in enumerate(ts.reader(path, header="ordered")):
        ccrs['gerp'].append(np.mean(map(float,region['GERP'].split(","))))
        ccrs['phast'].append(np.mean(map(float,region['phastCons'].split(","))))
        ccrs['cadd'].append(np.mean(map(float,region['CADD'].split(","))))
        length=sum([int(i.split("-")[1])-int(i.split("-")[0]) for i in region['ranges'].split(',')])
        ccrs['pct'].append(float(region['weighted_pct'])); ccrs['gene'].append(region['gene']); ccrs['chrom'].append(region['chrom']); ccrs['ranges'].append(region['ranges']); ccrs['length'].append(length)
        if genes[region['gene']]:
            genes[region['gene']][0]+=1
            genes[region['gene']][1]+=length
        else:
            genes[region['gene']]=[1,length]
    var['ccrs']=ccrs
    var['genes']=genes
    return var

def barplot(dataset):
    mi,ma=np.nanmin(dataset),np.nanmax(dataset)
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
    var=read_values('scoredregions.bed')
    #clingenes=clinvar('patho.vcf') #patho.vcf, vcf file after clinvarfilter.py for pathogenic variants
    pickle.dump(var, open("var.pickle", "wb"))

var = pickle.load(open("var.pickle", "rb"))
ccrs = var['ccrs']; genes = var['genes']
gerp=ccrs['gerp']; phast=ccrs['phast']; cadd=ccrs['cadd']; ccrpct=ccrs['pct']
toptengerp,toptenphast,toptencadd,topgerp,topphast,topcadd=[],[],[],[],[],[]
for region in ccrs:
    if ccrs['pct'] >= 90:
        toptengerp.append(ccrs['gerp'])
        toptenphast.append(ccrs['phast'])
        toptencadd.append(ccrs['cadd'])
        if ccrs['pct'] >= 99:
            topgerp.append(ccrs['gerp'])
            topphast.append(ccrs['phast'])
            topcadd.append(ccrs['cadd'])

meangerp=np.nanmean(gerp)
meanphast=np.nanmean(phast)
meancadd=np.nanmean(cadd)
meanccr=np.nanmean(ccrpct)

prop=[]
f1=open('geneexlength.txt', 'r') #sed 's/\"//g' exacresiduals/flatexome.bed | sed 's/;//g' | awk '{a[$4]+=($3-$2)} END {for (i in a) print i"\t"a[i]}' > geneexlength.txt
for line in f1:
    line=line.strip().split('\t')
    if line[0] in genes.keys():
        genes[line[0]].append(line[1])
        
    
f, ax = plt.subplots()
for i in genes:
    nums=genes[i]
    if nums[0]/float(nums[2])>.26:
        print i, nums
    prop.append(nums[0]/float(nums[2]))
mi=0; ma = max(prop); rng = (mi,ma)
p,p_edges=np.histogram(prop, bins=40, range=rng) #bins=400
p=map(lambda x: float(x)/sum(p), p)
width_p = (p_edges[1]-p_edges[0])
ax.bar(p_edges[:-1], p, width = width_p, color = 'r', log = True)
ax.xaxis.set_ticks(p_edges)#[1::2])
#ax1.set_xlim(0,115)
labs=ax.xaxis.set_ticklabels(p_edges, rotation=45, ha='right', fontname='Cmr10')
ticks=ax.get_xticks()
ax.xaxis.set_ticklabels(['{:.1f}%'.format(x*100) for x in ticks], fontname='Cmr10')
ticks=ax.get_yticks()
ax.yaxis.set_ticklabels(['{:f}%'.format(x*100) for x in ticks], fontname='Cmr10')
ax.set_xlabel('Number of regions in gene/size of gene')
ax.set_ylabel('Proportion of genes in bin')
plt.tight_layout()
plt.savefig('geneproportions.pdf', bbox_inches='tight')

def positive_selection(ccrs):
    key=None;gerpsum=0;phastsum=0;caddsum=0;gerpgene=set();phastgene=set();caddgene=set();
    for i in range(0,len(ccrs['pct'])):
        if key == None or key != [ccrs['pct'][i],ccrs['chrom'][i],ccrs['gene'][i],ccrs['ranges'][i]]:
            if ccrs['gerp'][i] > meangerp and ccrs['pct'][i] < meanccr and ccrs['length'][i] > 10:
                gerpsum+=1; gerpgene.add(ccrs['gene'][i])
            if ccrs['phast'][i] > meanphast and ccrs['pct'][i] < meanccr and ccrs['length'][i] > 10:
                phastsum+=1; phastgene.add(ccrs['gene'][i])
            if ccrs['cadd'][i] > meancadd and ccrs['pct'][i] < meanccr and ccrs['length'][i] > 10:
                caddsum+=1; caddgene.add(ccrs['gene'][i])
        key=[ccrs['pct'][i],ccrs['chrom'][i],ccrs['gene'][i],ccrs['ranges'][i]]
    print gerpsum,len(gerpgene),phastsum,len(phastgene),caddsum,len(caddgene)
positive_selection(ccrs)

fig, axarr = plt.subplots(3)
fig.tight_layout()
#color = ["white" if x[0] in clingenes else "red" for x in topccr]
#color = ["red" if x < meanccr else "blue" if x > 90 else "black" for x in ccrpct]
#ind = [i for i,x in enumerate(ccrpct) if x[0] not in clingenes]
#ccrpct = [ccrpct[i] for i in ind]; gerp=[gerp[i] for i in ind]; phast=[phast[i] for i in ind]; cadd=[cadd[i] for i in ind]
        
axarr[0].set_xlim(0,100)
axarr[1].set_xlim(0,100)
axarr[2].set_xlim(0,100)
g=axarr[0].hexbin(ccrpct, gerp, cmap='rainbow', bins='log')
axarr[0].set_xlabel('CCR Percentile');axarr[0].set_ylabel('GERP')
cbar = fig.colorbar(g, ax=axarr[0], orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
axarr[0].legend(loc='best')
axarr[0].add_patch(patches.Rectangle((90, min(gerp)), 10, meangerp-min(gerp), fill=False, edgecolor='red', linewidth=2))
axarr[0].add_patch(patches.Rectangle((0, meangerp), meanccr, max(gerp)-meangerp, fill=False, edgecolor='black', linewidth=2))

p=axarr[1].hexbin(ccrpct, phast, cmap='rainbow', bins='log')
axarr[1].set_xlabel('CCR Percentile');axarr[1].set_ylabel('phastCons')
cbar = fig.colorbar(p, ax=axarr[1], orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
axarr[1].legend(loc='best')
axarr[1].add_patch(patches.Rectangle((90, min(phast)), 10, meanphast-min(phast), fill=False, edgecolor='red', linewidth=2))
axarr[1].add_patch(patches.Rectangle((0, meanphast), meanccr, max(phast)-meanphast, fill=False, edgecolor='black', linewidth=2))

c=axarr[2].hexbin(ccrpct, cadd, cmap='rainbow', bins='log')
axarr[2].set_xlabel('CCR Percentile');axarr[2].set_ylabel('CADD')
cbar = fig.colorbar(c, ax=axarr[2], orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
axarr[2].legend(loc='best')
axarr[2].add_patch(patches.Rectangle((90, min(cadd)), 10, meancadd-min(cadd), fill=False, edgecolor='red', linewidth=2))
axarr[2].add_patch(patches.Rectangle((0, meancadd), meanccr, max(cadd)-meancadd, fill=False, edgecolor='black', linewidth=2))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=None)
plt.savefig('consvsconst.pdf', format='pdf', bbox_inches='tight')

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
axarr[0,0].axvline(meangerp, color='magenta', linestyle='dashed', linewidth=2)
axarr[1,0].axvline(meangerp, color='magenta', linestyle='dashed', linewidth=2)

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
axarr[0,1].axvline(meanphast, color='magenta', linestyle='dashed', linewidth=2)
axarr[1,1].axvline(meanphast, color='magenta', linestyle='dashed', linewidth=2)

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
axarr[0,2].axvline(meancadd, color='magenta', linestyle='dashed', linewidth=2)
axarr[1,2].axvline(meancadd, color='magenta', linestyle='dashed', linewidth=2)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=None)
plt.savefig('consdist.pdf', format='pdf', bbox_inches='tight')
print meangerp, meanphast, meancadd
