import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import numpy as np

f=open(sys.argv[1], "r")

gos={}
totalbp=0.0

for line in f:
    fields=line.strip().split("\t")
    ccr=float(fields[13]); bp=int(fields[2])-int(fields[1])
    goterm=fields[-1].split(";")[-2]
    try:
        if ccr > 90:
            gos[goterm][0]+=bp
        else:
            gos[goterm][1]+=bp
    except KeyError or IndexError:
        if ccr > 90:
            gos[goterm]=[bp,0]
        else:
            gos[goterm]=[0,bp]
    totalbp+=bp

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
#import seaborn
N = len(gos.keys())
ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
labels=[]; lows=[]; highs=[]
for key in gos:
    labels.append(key)
    highs.append(gos[key][0])
    lows.append(gos[key][1])

highs=map(lambda x: float(x)/sum(highs), highs)
lows=map(lambda x: float(x)/sum(lows), lows)
rects1 = ax.bar(ind, lows, width, color = 'b')
rects2 = ax.bar(ind + width, highs, width, color = 'r')
ax.set_ylabel("Total bp")
ax.set_title('GO Term Enrichment in CCRs overlapping Pfam domains of unknown function')
ax.set_xticks(ind + width/2)
ax.set_xticklabels(tuple(labels))
ax.legend((rects1[0], rects2[0]), ('Low CCRs (< 90%)', 'High CCRs (> 90%)'))
plt.savefig('unknowns.pdf',bbox_inches='tight')
