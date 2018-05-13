# script for pLI v GERP comparison
from bw import BigWig
import numpy as np
import tabix
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import pearsonr

def read_gerp(gerp, region):
    chrom = region[0]
    if "chr" not in chrom:
        chrom="chr"+chrom
    start=int(region[1]); end=int(region[2])
    gerps=np.frombuffer(gerp.values(chrom, start, end), dtype='f')
    return np.nanmean(gerps), len(gerps)

gerppath = sys.argv[1] #'/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
plipath = sys.argv[2] #see pathoscore
outfile = 'pligerp.pdf'

def perchrom(pli_gerp_chrom):
    plis = []; gerps = []
    plipath, gerppath, chrom = pli_gerp_chrom
    pLI = tabix.open(plipath)
    gerp = BigWig(gerppath)

    gerpdict={}; lengths=[]; scores=[]; rangeprev = None
    for region in pLI.querys(chrom):
        gerpscore, overlap = read_gerp(gerp, region)
        gerps.append(float(gerpscore))
        plis.append(float(region[-1]))

    return gerps, plis

import multiprocessing as mp
p = mp.Pool(12)

gerp = []; pli = []
for gerps, plis in p.imap_unordered(perchrom, ((plipath, gerppath, str(chrom)) for chrom in range(1, 23))):
   for g, p in zip(gerps, plis):
        gerp.append(g)
        pli.append(p)

fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(pli, gerp, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(round(10**x) if round(10**x,-1) <= 10 else round(10**x,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Genes', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(pli, gerp)
print (pr, pval)
plt.tight_layout()
sns.despine()
print ("Pearson's r: " + str(pr) + "\np-value: " + str(pval))
plt.xlabel('pLI')
plt.ylabel('Mean GERP++')
plt.savefig(outfile, format='pdf', bbox_inches='tight')
