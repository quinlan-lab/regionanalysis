import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
from matplotlib import pyplot as plt
import cPickle as pickle
import numpy as np
import seaborn as sns
sns.set_style('white')
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

import toolshed as ts
import sys

it = ts.reader(sys.argv[1]) # exacresiduals/results/2016_12_10/weightedresiduals.txt
iterable = (i for i in it)
totcpg,totcov,totresid,totpct = [],[],[],[]
topcpg,topcov,topresid,toppct=[],[],[],[]
lengths=[]; cpgs=[]; covs=[]; rangeprev = None
for region in iterable:
    cpg=float(region['cpg']);
    cov=float(region['cov_score']);
    resid=float(region['resid']);
    pct=float(region['weighted_pct']);
    ranges=region['ranges']
    if rangeprev is None:
        pass # will be dealt with at next pass through loop
    elif rangeprev != ranges:
        totcpg.append(cpgprev); totcov.append(covprev); totresid.append(residprev); totpct.append(pctprev)
        #if pct > 99:
        #    topcpg.append(cpgprev); topcov.append(covprev); topresid.append(residprev); toppct.append(pctprev)
    rangeprev = ranges; covprev = cov; cpgprev = cpg; residprev = resid; pctprev = pct

totcpg.append(cpgprev)
totcov.append(covprev)
totresid.append(residprev)
totpct.append(pctprev)

X = {"CpG": [], "intercept": []}
X['intercept'] = np.ones(len(totcov))
X['CpG'] = totcpg
X = pd.DataFrame(X)
results = sm.OLS(totcov, X, hasconst=True).fit()
intercept=results.params['intercept']; cpgcoef=results.params['CpG']

fig = plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 2])
ax0 = plt.subplot(gs[0])
ax0.set_ylim(0,max(totcov))
colors = [(1, 1, 1), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax0.hexbin(totcpg, totcov, cmap=cm, bins='log', alpha=0.5, mincnt=1)
ax0.set_xlabel('CpG proportion');ax0.set_ylabel('Sum of coverage fractions')
def y_format(x,y):
    return '{:.0f}'.format(10**x-1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
    
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax0, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#ax0.plot(topcpg,topcov,'b.')
x=[0,np.max(totcpg)]
y=[intercept,cpgcoef*np.max(totcpg)]
ax0.plot(x,y,'k-')
ax0.set_ylabel('Sum of fraction of exomes in ExAC\n with 10x coverage for each bp')
#ax1 = plt.subplot(gs[1])
#ax1.plot(totcpg,totresid,'r.')
#ax1.plot(topcpg,topresid,'b.')
#ax1.set_xlabel('CpG fraction')
#ax1.set_ylabel('Studentized Residuals')

plt.tight_layout()
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/regression.pdf', bbox_inches='tight')
plt.close()
