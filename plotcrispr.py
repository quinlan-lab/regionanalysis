import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import subprocess as sp
import argparse
import scipy.stats as ss

parser = argparse.ArgumentParser()
parser.add_argument("regions")
parser.add_argument("sgrnas")
parser.add_argument("-s", "--study-name", help="name of your study, for plot label", required=True)
args=parser.parse_args()

regions=args.regions
sgrnas=args.sgrnas
study_name=args.study_name

import plotdistro as pltdist
import pandas as pd

cs=[]; pct=[]
ixn = pltdist.intersect(regions,sgrnas,True) #sys.argv[2]='patho.bed'
for line in ixn.split("\n"):
    fields = line.strip().split("\t")
    pct.append(float(fields[13])); cs.append(float(fields[17]))

fig = plt.figure()
ax = fig.add_subplot(111)
xmi = np.min(cs)
xma = np.max(cs)
ymi = np.min(pct)
yma = np.max(pct)
print xmi,xma,ymi,yma

title="CRISPR score vs Residual Percentile for "+study_name+"\n"
ax.plot(cs,pct,'.')
X = {"cs": cs, "intercept": np.ones(len(pct))}
X = pd.DataFrame(X)

import statsmodels.api as sm

results = sm.OLS(pct, X, hasconst=True).fit()
x=[np.min(cs),np.max(cs)]
y=[results.params['intercept'],results.params['cs']*np.max(cs)]
ax.plot(x,y,'r-')
plt.ticklabel_format(useOffset=False)
ax.set_xlabel("CS score")
ax.set_ylabel("Percentile")
ax.set_xlim(xmi,xma)
ax.set_ylim(ymi,yma)
plt.title(title)
plt.savefig("plots/"+study_name+"crispr.png")
plt.close()
