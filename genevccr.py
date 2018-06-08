import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import pearsonr
import sys
import gzip

with open(sys.argv[1], 'r') as genefile:
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
with gzip.open(sys.argv[2], 'rb') as ccrfile:
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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 21
fig, ax = plt.subplots()
r, pval = pearsonr(lengths, ccrs)
ax.scatter(lengths, ccrs, alpha = 0.3, s = 1, label = "r = " + '{:.3f}'.format(r) + "\np-value = " + '{:.3g}'.format(pval), color='#388aac')
m, b = np.polyfit(lengths, ccrs, 1)
lengths=np.array(lengths)
ax.plot(lengths, m*lengths + b, 'k-')
plt.legend(loc='best')
plt.xlim(0,20000)
plt.ylim(0,100)
plt.tight_layout()
sns.despine()
plt.xlabel('Total gene coding sequence (bp)')
plt.ylabel('Maximum CCR percentile for gene')
plt.savefig(sys.argv[3], bbox_inches='tight')
