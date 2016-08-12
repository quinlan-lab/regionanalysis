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

variables = pickle.load(open("exacresiduals/var.pickle", "rb"))
resid=variables['resid']
raws=variables['rawresid']
ys=variables['cov']
cpg=variables['cpg']
genes=variables['genes']
fig = plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 2])
ax0 = plt.subplot(gs[0])
ax0.plot(cpg,ys,'.')
x=[0,np.max(cpg)]
y=[variables['intercept'],variables['cpgcoef']*np.max(cpg)]
ax0.plot(x,y,'r-')
ax0.set_ylabel('log(1.0 + sum(Coverage fractions))')
ax1 = plt.subplot(gs[1])
ax1.plot(cpg,resid,'r.')
ax1.set_xlabel('CpG fraction')
ax1.set_ylabel('Studentized Residuals')

plt.savefig("/uufs/chpc.utah.edu/common/home/u1021864/analysis/plots/regression.png")
plt.close()

resid_pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))
cov_pctile = 100.0 * np.sort(ys).searchsorted(ys) / float(len(ys))

