import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
from collections import defaultdict
import numpy as np
import seaborn as sns
sns.set_style('white')

f95=open(sys.argv[1],'r')
f99=open(sys.argv[2],'r')

def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % float(height),
                ha='center', va='bottom')

# 95th PCT

ccrs=defaultdict(int)
bins=defaultdict(int)
for line in f95:
    fields=line.strip().split("\t")
    ccr=fields[3]+fields[6]+fields[13] # gene,ranges,weighted_pct
    if fields[-1] != "0":
        ccrs[ccr] += 1
    else:
        ccrs[ccr] += 0
    
print "Last bin total (95):"
tot=0
for ccr in ccrs:
    if ccrs[ccr] == 0:
        bins[0] += 1
    if ccrs[ccr] == 1:
        bins[1] += 1
    if ccrs[ccr] == 2:
        bins[2] += 1
    if ccrs[ccr] == 3:
        bins[3] += 1
    if ccrs[ccr] == 4:
        bins[4] += 1
    if ccrs[ccr] == 5:
        bins[5] += 1
    if ccrs[ccr] == 6:
        bins[6] += 1
    if ccrs[ccr] == 7:
        bins[7] += 1
    if ccrs[ccr] == 8:
        bins[8] += 1
    if ccrs[ccr] == 9:
        bins[9] += 1
    if ccrs[ccr] >= 10:
        bins["10+"] += 1
        tot+=ccrs[ccr]

print tot
print bins

keys=bins.keys(); values=bins.values()
fig, ax = plt.subplots(1)
width=0.4
lefts=np.arange(0,1.2*len(keys),1.2)
keys, values = zip(*sorted(zip(keys,values)))

rects=ax.bar(left=lefts,height=values,width=width,tick_label=keys,color="#A1DAD7")
autolabel(rects, ax)

# 99th PCT

ccrs=defaultdict(int)
bins=defaultdict(int)
for line in f99:
    fields=line.strip().split("\t")
    ccr=fields[3]+fields[6]+fields[13] # gene,ranges,weighted_pct
    if fields[-1] != "0":
        ccrs[ccr] += 1
    else:
        ccrs[ccr] += 0
    
print "Last bin total (99):"
tot=0
for ccr in ccrs:
    if ccrs[ccr] == 0:
        bins[0] += 1
    if ccrs[ccr] == 1:
        bins[1] += 1
    if ccrs[ccr] == 2:
        bins[2] += 1
    if ccrs[ccr] == 3:
        bins[3] += 1
    if ccrs[ccr] == 4:
        bins[4] += 1
    if ccrs[ccr] == 5:
        bins[5] += 1
    if ccrs[ccr] == 6:
        bins[6] += 1
    if ccrs[ccr] == 7:
        bins[7] += 1
    if ccrs[ccr] == 8:
        bins[8] += 1
    if ccrs[ccr] == 9:
        bins[9] += 1
    if ccrs[ccr] >= 10:
        bins["10+"] += 1
        tot+=ccrs[ccr]

print tot
print bins

keys=bins.keys(); values=bins.values()
width=0.4
alefts=np.arange(0+width,1.2*len(keys)+width, 1.2)
keys, values = zip(*sorted(zip(keys,values)))

rects2=ax.bar(left=alefts,height=values,width=width,tick_label=keys,color="#388aac")
autolabel(rects2, ax)
ax.set_xlabel("Number of ClinVar pathogenic variants in 95th and 99th percentile CCR")
ax.set_ylabel("Frequency")
lims=ax.get_ylim()
#ax.set_ylim(lims[0]/1.5, lims[1]+.4)
ax.set_yscale("log",basey=10,nonposy="clip")
ax.set_xticks(lefts-width*.5)
sns.despine(bottom=True)
plt.tight_layout()
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/pathogenicsinccrs9599.pdf', bbox_inches='tight')
