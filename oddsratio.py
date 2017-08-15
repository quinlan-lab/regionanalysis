import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from collections import defaultdict
import sys

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

patho = open(sys.argv[1], "r") #patho-ccr.txt or neurodev-ccr.txt
benign = open(sys.argv[2], "r") #benign-ccr.txt or control-ccr.txt
filename = sys.argv[3]

bins=[0,80,90,95,99,100]
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
        
keys,ors,props=[],[],[]
for i in presults:
    try:
        odds=(presults[i][0]/presults[i][1])/(bresults[i][0]/bresults[i][1])
    except ZeroDivisionError:
        odds=(presults[i][0]/presults[i][1])/(1/bresults[i][1])
    keys.append(i); ors.append(odds)
    props.append(lens[i]/totlen*100)
    print i, presults[i], bresults[i]
    print odds
    print lens[i]/totlen*100

#fig, axarr = plt.subplots(2, sharex=True)
fig, ax = plt.subplots(1)
width=0.6
lefts=np.arange(0,.8*len(keys),.8)
ticks=[key+'\n'+'%.2f' % prop+"%" for key, prop in zip(keys, props)]
ticks, ors = zip(*sorted(zip(ticks,ors)))
ax.bar(lefts,ors,width=width,tick_label=ticks,alpha=0.5)
ax.set_xlabel("CCR Percentile Bins"+'\n'+"Percent of CCRs that score these variants"+'\n'+"(not whole exome, 95-100 is ~4.1% of all CCRs)")
ax.set_ylabel("Odds Ratio")
ax.set_title(filename[0].upper()+filename[1:])
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/oddsratio'+filename+'.pdf', bbox_inches='tight')
