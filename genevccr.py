import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import pearsonr

import gzip

with open('exacresiduals/flatexome.bed', 'r') as genefile:
    data = []
    for line in genefile:
        line = line.strip().split("\t")
        data.append(line)

data.sort(key=lambda s:s[3])
prevgene=None
genes={}
lengths=[]
for gene in data:
    if prevgene is None:
        pass
    elif prevgene[3] != gene[3]:
        genes[prevgene[3]]=[sum(lengths),0]
        lengths=[]
    else:
        lengths.append(float(prevgene[2])-float(prevgene[1]))
    prevgene=gene
genes[prevgene[3]]=[sum(lengths),0]
        
        
        
data=[]
with gzip.open('essentials/gnomadbased-ccrs.bed.gz', 'rb') as ccrfile:
    for line in ccrfile:
        line = line.strip().split("\t")
        ccr = float(line[-1])
        try:
            if genes[line[3]][1] < ccr:
                genes[line[3]][1] = ccr
        except KeyError:
            pass

lengths,ccrs=[],[]
for i in genes:
    if genes[i][1] != 0:
        lengths.append(genes[i][0])
        ccrs.append(genes[i][1])

fig, ax = plt.subplots()
r, pval = pearsonr(ccrs, lengths)
ax.scatter(ccrs, lengths, alpha = 0.3, s = 1, label = "r = " + str(r) + "\np-value = " + str(pval))
m, b = np.polyfit(ccrs, lengths, 1)
ccrs=np.array(ccrs)
ax.plot(ccrs, m*ccrs + b, 'k-')
plt.legend(loc='best')
#plt.ylim(0,20000)
# colors = [(0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (0, 0, 0)]
# cmap_name='mymap'
# cm = LinearSegmentedColormap.from_list(cmap_name, colors)
# g=ax.hexbin(ccrs, lengths, cmap=cm, alpha=0.5, mincnt=1) #, bins='log')
# def y_format(x,y):
   # return '{:.0f}'.format(round(10**x-1,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
# counts,edges=np.histogram(g.get_array(),bins=8)
# cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges) #, format=FuncFormatter(y_format))
# cbar.set_label('Number of Genes', rotation=270, labelpad=20)
plt.tight_layout()
sns.despine()
plt.xlabel('Maximum CCR Percentile for Gene')
plt.ylabel('Total Gene Coding Sequence')
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/genevccr.pdf', bbox_inches='tight')
