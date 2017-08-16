import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from collections import defaultdict
import sys
import seaborn as sns
sns.set_style('white')
from scipy.stats import fisher_exact

##make contingency table for odds ratios like:
##|============|============|===============
##|            |   In bin   |  Not in bin
##|------------|------------|---------------
##| Pathogenic |    Pb      |     Pn
##|------------|------------|---------------
##| Benign     |    Bb      |     Bn
##|============|============|===============
#
#                 Pb/Pn
# "Odds Ratio" = -------
#                 Bb/Bn

def conf_interval(ratio, std_error):
    """
    Calculate 95% confidence interval for odds ratio and relative risk.
    """

    _lci = np.log(ratio) - 1.96*std_error
    _uci = np.log(ratio) + 1.96*std_error

    lci = round(np.exp(_lci), 2)
    uci = round(np.exp(_uci), 2)

    return (lci, uci)

patho = open(sys.argv[1], "r") #patho-ccr.txt or neurodev-ccr.txt
benign = open(sys.argv[2], "r") #benign-ccr.txt or control-ccr.txt
filename = sys.argv[3]

bins=[0,20,80,90,95,100]
#bins=[0,80,90,100]
presults=defaultdict(lambda:[0.0,0.0])
bresults=defaultdict(lambda:[0.0,0.0])
lens=defaultdict(lambda:0.0)
pset,bset=set(),set()
totlen=0.0
for line in patho:
    fields=line.strip().split("\t")
    info = fields[7]; ccr=fields[8:]
    if "_exclude" in info:
        continue
    ccrkey="\t".join(ccr)
    if ccrkey not in pset or bset:
        leng=int(ccr[2])-int(ccr[1])
        totlen+=leng
    else:
        leng=0.0
        pset.add(ccrkey)
    pct=float(ccr[-1])
    for i, bin in enumerate(bins):
        if i == len(bins)-1:
            break
        if pct >= bins[i] and pct <= bins[i+1]:
            presults[str(bins[i])+"-"+str(bins[i+1])][0]+=1.0
            lens[str(bins[i])+"-"+str(bins[i+1])]+=leng
        else:
            presults[str(bins[i])+"-"+str(bins[i+1])][1]+=1.0

for line in benign:
    fields=line.strip().split("\t")
    info = fields[7]; ccr=fields[8:]
    ccrkey = "\t".join(ccr)
    if ccrkey not in pset or bset:
        leng=int(ccr[2])-int(ccr[1])
        totlen+=leng
    else:
        bset.add(ccrkey)
        leng=0.0
    pct=float(ccr[-1])
    leng=int(ccr[2])-int(ccr[1])
    for i, bin in enumerate(bins):
        if i == len(bins)-1:
            break
        if pct >= bins[i] and pct <= bins[i+1]:
            bresults[str(bins[i])+"-"+str(bins[i+1])][0]+=1.0
            lens[str(bins[i])+"-"+str(bins[i+1])]+=leng
        else:
            bresults[str(bins[i])+"-"+str(bins[i+1])][1]+=1.0
        
keys,ors,props,cis=[],[],[],[]
for i in presults:
    a,b,c,d=presults[i][0],bresults[i][0],presults[i][1],bresults[i][1]
    try:    
        odds=(a*d)/(b*c)
    except ZeroDivisionError:
        odds=(a*d)/(1*c)
    or_se=np.sqrt((1/a)+(1/b)+(1/c)+(1/d))
    or_ci = conf_interval(odds, or_se)
    keys.append(i); ors.append(odds)
    props.append(lens[i]/totlen*100); cis.append(or_ci)
    print i, presults[i], bresults[i]
    print odds
    print or_ci
    print lens[i]/totlen*100

def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        if height==1: # 2
            height=rect.get_y()
            ax.text(rect.get_x() + rect.get_width()/2., .65*height,
                    '%.2f' % float(height),
                    ha='center', va='bottom')
            continue
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%.2f' % float(height),
                ha='center', va='bottom')

#fig, axarr = plt.subplots(2, sharex=True)
fig, ax = plt.subplots(1)
width=0.4
lefts=np.arange(0,.6*len(keys),.6)
ticks=[key+'\n'+'%.2f' % prop+"%" for key, prop in zip(keys, props)]
ticks, ors, cis = zip(*sorted(zip(ticks,ors,cis)))
y_r = [[ors[i]-cis[i][0] for i in range(len(cis))],[cis[i][1]-ors[i] for i in range(len(cis))]]
print cis, y_r
ax.errorbar([i for i in lefts], ors, yerr=y_r, color='k', capsize=4, fmt="none")
bottoms = [height if height<1 else 1 for height in ors] #2
ors2 = [1 if height<1 else height for height in ors] #2
rects=ax.bar(lefts,ors2,width=width,bottom=bottoms,tick_label=ticks)
#autolabel(rects, ax)
ax.axhline(y=1, color='k') # 2
ax.set_xlabel("CCR Percentile Bins"+'\n'+"Percent of CCRs that score these variants"+'\n'+"(not whole exome, 95-100 is ~4.1% of all CCRs)")
ax.set_ylabel("Odds Ratio")
ax.set_title(filename[0].upper()+filename[1:])
ax.set_yscale("log",basey=10,nonposy="clip") # 2
lims=ax.get_ylim()
ax.set_ylim(lims[0]/1.5, lims[1]+.4)
sns.despine()
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/oddsratio'+filename+'.pdf', bbox_inches='tight')
