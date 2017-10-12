import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import seaborn as sns
sns.set_style('white')
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence
import pandas as pd
from collections import defaultdict

ccrint = open(sys.argv[1], 'r') # all variants ccr file intersected with doubleton<= CCR file
plotout = sys.argv[2] # pdf location

ccrs, dccrs = [], []
maxccrs, medccrs, maxdccrs = [], [], []
ccrprev, dccrprev, ccrcalc = None, None, []
for line in ccrint:
    fields = line.strip().split("\t")
    # may need to fix if I change the number of fields in CCR file
    dccr = float(fields[13]) # doubleton <= CCRs score
    ccr = float(fields[-1]) # all CCRs score
    gene = fields[3]
    ccrs.append(ccr); dccrs.append(dccr)
    if dccrprev is None:
        pass
    elif dccrprev == dccr:
        ccrcalc.append(ccrprev)
    elif dccrprev != dccr:
        ccrcalc.append(ccrprev)
        maxccrs.append(max(ccrcalc))
        medccrs.append(np.median(ccrcalc))
        maxdccrs.append(dccrprev)
        ccrcalc = []
    ccrprev = ccr; dccrprev = dccr; lineprev = line
ccrcalc.append(ccrprev)
maxccrs.append(max(ccrcalc))
medccrs.append(np.median(ccrcalc))
maxdccrs.append(dccrprev)

# all intersections of ccrs with doubleton ccrs

X = {"CCRs": [], "intercept": []}
X['intercept'] = np.ones(len(dccrs))
X['CCRs'] = ccrs
X = pd.DataFrame(X)
results = sm.OLS(dccrs, X.astype(float), hasconst=True).fit()
intercept=results.params['intercept']; cpgcoef=results.params['CCRs']

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
ax.set_ylim(0,max(ccrs))
colors = [(1, 1, 1), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(dccrs, ccrs, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(10**x-1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
    
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
x=[0,np.max(dccrs)]
y=[intercept,cpgcoef*np.max(dccrs)]
ax.plot(x,y,'k-')
ax.set_xlabel('Doubleton or greater CCRs');ax.set_ylabel('ALL CCRs')
plt.tight_layout()
plt.savefig(plotout, bbox_inches='tight')
plt.close()

# max ccr intersection with doubleton ccr

X = {"CCRs": [], "intercept": []}
X['intercept'] = np.ones(len(maxdccrs))
X['CCRs'] = maxccrs
X = pd.DataFrame(X)
results = sm.OLS(maxdccrs, X.astype(float), hasconst=True).fit()
intercept=results.params['intercept']; cpgcoef=results.params['CCRs']

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
ax.set_ylim(0,max(maxccrs))
colors = [(1, 1, 1), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(maxdccrs, maxccrs, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(10**x-1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
    
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
x=[0,np.max(maxdccrs)]
y=[intercept,cpgcoef*np.max(maxdccrs)]
ax.plot(x,y,'k-')
ax.set_xlabel('Doubleton or greater CCRs');ax.set_ylabel('Max of ALL CCRs')
plt.tight_layout()
outsplit=plotout.rpartition(".")
outfile="".join(outsplit[0:-2])+'max'+"".join(outsplit[-2:])
plt.savefig(outfile, bbox_inches='tight')
plt.close()

# max ccr intersection with doubleton ccr

X = {"CCRs": [], "intercept": []}
X['intercept'] = np.ones(len(maxdccrs))
X['CCRs'] = medccrs
X = pd.DataFrame(X)
results = sm.OLS(maxdccrs, X.astype(float), hasconst=True).fit()
intercept=results.params['intercept']; cpgcoef=results.params['CCRs']

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
ax.set_ylim(0,max(medccrs))
colors = [(1, 1, 1), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(maxdccrs, medccrs, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(10**x-1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
    
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
x=[0,np.max(maxdccrs)]
y=[intercept,cpgcoef*np.max(maxdccrs)]
ax.plot(x,y,'k-')
ax.set_xlabel('Doubleton or greater CCRs');ax.set_ylabel('Median of ALL CCRs')
plt.tight_layout()
outsplit=plotout.rpartition(".")
outfile="".join(outsplit[0:-2])+'med'+"".join(outsplit[-2:])
plt.savefig(outfile, bbox_inches='tight')
plt.close()

# plot of FDR by bin
bins=range(0,101,5)
#bindict=defaultdict(float)
bindict=defaultdict(lambda : [0,0,0.0])
for maxccr, dccr in zip(maxccrs, maxdccrs):
    for i in range(0,len(bins)):
        if i == 0:
            if maxccr == bins[0]:
                bindict[bins[0]][0] += 0.0 # no FDR for 0th percentile
            if dccr == bins[0]:
                bindict[bins[0]][1] += 1.0
        elif maxccr <= bins[i] and maxccr > bins[i-1]:
            if maxccr < dccr:
                bindict[bins[i]][0] += 1.0
        if dccr <= bins[i] and dccr > bins[i-1]:
            bindict[bins[i]][1] += 1.0

bins, vals = [], []
print "\t".join(["Bin", "doubleton CCRs > max full CCR score", "count of doubleton CCRs in bin"])
for bin in bindict:
    vals.append(bindict[bin][0]/bindict[bin][1])
    bins.append(bin)
    print "\t".join(map(str,[bin, bindict[bin][0], bindict[bin][1]]))
bins, vals = zip(*sorted(zip(bins, vals)))
matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
binar=[i-5 for i in bins]
width = (bins[1]-bins[0])+.1
ax.bar(bins, vals, width = width, color = 'b', alpha = 0.7)
ax.set_xlabel('Doubleton or greater CCR percentile');ax.set_ylabel('FDR (fraction of doubleton CCRs > max full CCR score/count of doubleton CCRs in bin)')
outsplit=plotout.rpartition(".")
outfile="".join(outsplit[0:-2])+'bar'+"".join(outsplit[-2:])
plt.savefig(outfile, bbox_inches='tight')

# max ccr intersection with doubleton ccr > 95%

X = {"CCRs": [], "intercept": []}
X['intercept'] = np.ones(len(maxdccrs))
X['CCRs'] = maxccrs
X = pd.DataFrame(X)
results = sm.OLS(maxdccrs, X.astype(float), hasconst=True).fit()
intercept=results.params['intercept']; cpgcoef=results.params['CCRs']

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
colors = [(1, 1, 1), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
cmap_name='mymap'
maxdccrs=np.array(maxdccrs); maxccrs=np.array(maxccrs)
ds=maxdccrs[np.where(maxdccrs>=95)]
cs=maxccrs[np.where(maxdccrs>=95)]
g=ax.plot(ds, cs, '.', alpha=0.5)
def y_format(x,y):
    return '{:.0f}'.format(10**x-1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
    
ax.set_xlabel('Doubleton or greater CCRs');ax.set_ylabel('Max of ALL CCRs')
plt.tight_layout()
outsplit=plotout.rpartition(".")
outfile="".join(outsplit[0:-2])+'zoomed'+"".join(outsplit[-2:])
plt.savefig(outfile, bbox_inches='tight')
plt.close()
