import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from collections import defaultdict
import sys
import seaborn as sns
sns.set_style('white')
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--varfiles", help = "variant intersection files for plot", nargs = 2)
parser.add_argument("-a", "--adgenefiles", help = "AD gene variant intersection files for AD gene part of plot, optional", nargs = 2)
parser.add_argument("-o", "--output", help = "output file name for plot")
parser.add_argument("-l", "--labels", nargs = 2, help = "labels for plot")
args = parser.parse_args()
varfiles = args.varfiles
adgenefiles = args.adgenefiles
filename = args.output
labels = args.labels
if labels:
    label=labels[0]
    adlabel=labels[-1]
else:
    label="All Genes"
    adlabel="AD Genes"
patho = open(varfiles[0], "r")
benign = open(varfiles[1], "r")
if adgenefiles:
    pathoad = open(adgenefiles[0], "r")
    benignad = open(adgenefiles[1], "r")

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

    lci = round(np.exp(_lci), 3)
    uci = round(np.exp(_uci), 3)

    return (lci, uci)

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

def generateresults(infile, results, rset, lens, totlen, bins, pathogenic):
    for line in infile:
        fields=line.strip().split("\t")
        info = fields[7]; ccr=fields[8:]
        if pathogenic and "_exclude" in info:
            continue
        ccrkey="\t".join(ccr)
        if ccrkey not in rset:
            leng=int(ccr[2])-int(ccr[1])
            totlen+=leng
        else:
            leng=0.0
            rset.add(ccrkey)
        pct=float(ccr[-1])
        for i, bin in enumerate(bins):
            if i == len(bins)-1:
                break
            if pct >= bins[i] and pct <= bins[i+1]:
                results[str(bins[i])+" - "+str(bins[i+1])][0]+=1.0
                lens[str(bins[i])+" - "+str(bins[i+1])]+=leng
            else:
                results[str(bins[i])+" - "+str(bins[i+1])][1]+=1.0
    return results, rset, lens, totlen

bins=[0,20,80,90,95,100]
presults=defaultdict(lambda:[0.0,0.0])
bresults=defaultdict(lambda:[0.0,0.0])
lens=defaultdict(lambda:0.0)
rset=set()
totlen=0.0
pathogenic = True
presults, rset, lens, totlen = generateresults(patho, presults, rset, lens, totlen, bins, pathogenic)
pathogenic = False
bresults, rset, lens, totlen = generateresults(benign, bresults, rset, lens, totlen, bins, pathogenic)

keys,ors,props,cis=[],[],[],[]
for i in presults:
    a,b,c,d=presults[i][0],bresults[i][0],presults[i][1],bresults[i][1]
    if 0 == b: b=1
    if 0 == c: c=1
    odds=(a*d)/(b*c)
    if 0 == a: a=1
    if 0 == d: d=1
    or_se=np.sqrt((1/a)+(1/b)+(1/c)+(1/d))
    or_ci = conf_interval(odds, or_se)
    keys.append(i); ors.append(odds)
    props.append(lens[i]/totlen*100); cis.append(or_ci)
    print i, presults[i], bresults[i]
    print odds
    print or_ci
    print lens[i]/totlen*100

if adgenefiles:
    padresults=defaultdict(lambda:[0.0,0.0])
    badresults=defaultdict(lambda:[0.0,0.0])
    lens=defaultdict(lambda:0.0)
    rset=set()
    totlen=0.0
    pathogenic = True
    padresults, rset, lens, totlen = generateresults(pathoad, padresults, rset, lens, totlen, bins, pathogenic)
    pathogenic = False
    badresults, rset, lens, totlen = generateresults(benignad, badresults, rset, lens, totlen, bins, pathogenic)

    akeys,aors,aprops,acis=[],[],[],[]
    for i in padresults:
        a,b,c,d=padresults[i][0],badresults[i][0],padresults[i][1],badresults[i][1]
        if 0 == b: b=1
        if 0 == c: c=1
        odds=(a*d)/(b*c)
        if 0 == a: a=1
        if 0 == d: d=1
        or_se=np.sqrt((1/a)+(1/b)+(1/c)+(1/d))
        or_ci = conf_interval(odds, or_se)
        akeys.append(i); aors.append(odds)
        aprops.append(lens[i]/totlen*100); acis.append(or_ci)
        print i, padresults[i], badresults[i]
        print odds
        print or_ci
        print lens[i]/totlen*100

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(1)

width=0.4
if adgenefiles:
    lefts=np.arange(0,1.2*len(keys),1.2)
else:
    lefts=np.arange(0,0.8*len(keys),0.8)
ticks=[key+'\n'+'(%.2f' % prop+"%)" for key, prop in zip(keys, props)]
ticks, ors, cis = zip(*sorted(zip(ticks,ors,cis)))
y_r = [[ors[i]-cis[i][0] for i in range(len(cis))],[cis[i][1]-ors[i] for i in range(len(cis))]]
print cis, y_r
ax.errorbar([i for i in lefts], ors, yerr=y_r, color='k', capsize=4, fmt="none")
bottoms = [height if height<1 else 1 for height in ors] #2 for log2
ors2 = [1-height if height<1 else height-1 for height in ors] #2 for log2
print ors2, bottoms, ticks
rects=ax.bar(left=lefts,height=ors2,width=width,bottom=bottoms,tick_label=ticks,color=(161/255.0,218/255.0,215/255.0), edgecolor=(96/255.0, 133/255.0, 131/255.0),label=label)
#autolabel(rects, ax)
ax.axhline(y=1, color='k') # 2 for log2
ax.set_title(filename[0].upper()+filename[1:])
if max(ors) > 10:
    ax.set_yscale("log",basey=10,nonposy="clip") # 2 for log2
for i,o in enumerate(ors2):
    ax.text(lefts[i]-width*.25, 0.14 if o<1 else 5, '%.3f' % ors[i] if o<1 else '%.1f' % ors[i], color='k', size=8) #2 for log2

if adgenefiles:
    alefts=np.arange(0+width,1.2*len(akeys)+width, 1.2)
    aticks=[key+'\n'+'(%.2f' % prop+"%)" for key, prop in zip(akeys, aprops)]
    aticks, aors, acis = zip(*sorted(zip(aticks,aors,acis)))
    y_r = [[aors[i]-acis[i][0] for i in range(len(acis))],[acis[i][1]-aors[i] for i in range(len(acis))]]
    print acis, y_r
    ax.errorbar([i for i in alefts], aors, yerr=y_r, color='k', capsize=4, fmt="none")
    abottoms = [height if height<1 else 1 for height in aors] #2 for log2
    aors2 = [1-height if height<1 else height-1 for height in aors] #2 for log2
    print aors2, abottoms, aticks
    rects=ax.bar(left=alefts,height=aors2,width=width,bottom=abottoms,tick_label=aticks,color=(56/255.0,138/255.0,172/255.0),edgecolor=(96/255.0, 133/255.0, 131/255.0),label=adlabel)
    ax.set_xticks(alefts-width*.5)
    #autolabel(rects, ax)
    ax.set_title(filename[0].upper()+filename[1:])
    ax.legend(loc='upper left')
    if max(aors) > 10:
        ax.set_yscale("log",basey=10,nonposy="clip") # 2 for log2
        for i,o in enumerate(aors2):
            ax.text(alefts[i]-width*.25, 0.14 if o<1 else 5, '%.3f' % aors[i] if o<1 else '%.1f' % aors[i], color='k', size=8) #2 for log2

lims=ax.get_ylim()
ax.set_ylim(lims[0]/1.5, lims[1]+.4)
ax.set_xlabel("CCR Percentile Bins"+'\n'+"Percent of CCRs that score these variants"+'\n'+"(not whole exome, 95-100 is ~4.1% of all CCRs)")
ax.set_ylabel("Odds Ratio\n(pathogenic versus benign)")
sns.despine()
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/oddsratio'+filename+'.pdf', bbox_inches='tight')
