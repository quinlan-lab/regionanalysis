import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
from matplotlib import pyplot as plt
import cPickle as pickle
import numpy as np
import seaborn as sns
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
for gene in iterable:
    cpg=float(gene['cpg']);
    cov=float(gene['cov_score']);
    resid=float(gene['cov_cpg_resid']);
    pct=float(gene['weighted_pct'])
    totcpg.append(cpg); totcov.append(cov); totresid.append(resid); totpct.append(pct)
    if pct > 99:
        topcpg.append(cpg); topcov.append(cov); topresid.append(resid); toppct.append(pct)
#variables = pickle.load(open("exacresiduals/var.pickle", "rb"))
#resid=variables['resid']
#raws=variables['rawresid']
#cov=variables['cov']
#cpg=variables['cpg']
#genes=variables['genes']
X = {"CpG": [], "intercept": []}
X['intercept'] = np.ones(len(totcov))
X['CpG'] = totcpg
X = pd.DataFrame(X)
results = sm.OLS(totcov, X, hasconst=True).fit()
intercept=results.params['intercept']; cpgcoef=results.params['CpG']

fig = plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 2])
ax0 = plt.subplot(gs[0])
#ax0.plot(totcpg,totcov,'g.')
ax0.set_ylim(0,max(totcov))
colors = [(1, 1, 1), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax0.hexbin(totcpg, totcov, cmap=cm, bins='log', alpha=0.5)
ax0.set_xlabel('CpG proportion');ax0.set_ylabel('Sum of coverage fractions')
def y_format(x,y):
    return '{:.0f}'.format(10**x-1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
    
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax0, orientation='vertical', extend='both', extendrect=True, drawedges=True, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270)
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
plt.savefig("plots/regression.pdf", bbox_inches='tight')
plt.close()
