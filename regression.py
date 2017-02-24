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
from matplotlib.ticker import FormatStrFormatter

import toolshed as ts
import sys

it = ts.reader(sys.argv[1])
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
g=ax0.hexbin(totcpg, totcov, cmap='rainbow', bins='log')
ax0.set_xlabel('CpG proportion');ax0.set_ylabel('Sum of coverage fractions')
cbar = fig.colorbar(g, ax=ax0, orientation='vertical', extend='both', format=FormatStrFormatter('$10^{%.1f}$'))
#ax0.plot(topcpg,topcov,'b.')
x=[0,np.max(totcpg)]
y=[intercept,cpgcoef*np.max(totcpg)]
ax0.plot(x,y,'r-')
ax0.set_ylabel('sum(Coverage fractions)')
#ax1 = plt.subplot(gs[1])
#ax1.plot(totcpg,totresid,'r.')
#ax1.plot(topcpg,topresid,'b.')
#ax1.set_xlabel('CpG fraction')
#ax1.set_ylabel('Studentized Residuals')

plt.tight_layout()
plt.savefig("plots/regression.pdf", bbox_inches='tight')
plt.close()
