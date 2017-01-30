import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import toolshed as ts

lens=[]
it = ts.reader(sys.argv[1]) #$HOME/analysis/exacresiduals/results/2016_12_10/unfilteredresiduals.txt
iterable = (i for i in it)
for region in iterable:
    lens.append(int(region['end'])-int(region['start']))

fig=plt.figure()
ax1 = fig.add_subplot(3, 1, 1)
ax2 = fig.add_subplot(3, 1, 2)
ax3 = fig.add_subplot(3, 1, 3)
ax1.hist(lens,bins=[1,2,3,4,5,6,7,8,9,10,20,50,100],facecolor='red')
ax2.hist(lens,bins=[50,100,500],facecolor='green')
ax3.hist(lens,bins=[500,1000,10000],facecolor='blue')
ax1.set_ylabel('Frequency')
ax2.set_ylabel('Frequency')
ax3.set_ylabel('Frequency')
ax2.set_xlim(50,500)
ax3.set_xlim(500,max(lens))
ax3.set_xlabel('Distance between functional variants')
plt.savefig('funchist.png',bbox_inches='tight')

totlen=0.0; ordered=[]; exomecov=[]; wgtdpct=[]
it = ts.reader(sys.argv[2]) #$HOME/analysis/exacresiduals/results/2016_12_10/weightedresiduals.txt
iterable = (i for i in it)
for region in iterable:
    length=int(region['end'])-int(region['start'])
    totlen+=length
    ordered.append((length,float(region['weighted_pct']),region['ranges']+region['gene']))
ordered=sorted(ordered,key=lambda x:x[1])
length=0.0
prevregion=None
for region in ordered:
    length+=region[0]
    if ordered[2]!=prevregion and prevregion!=None:
        exomecov.append(prevvalues[0])
        wgtdpct.append(prevvalues[1])
    prevregion=region[2]
    prevvalues=(length/totlen,region[1])
plt.figure()
plt.plot(wgtdpct,exomecov)
plt.xlabel('CCR Percentile (highest = most constrained)')
plt.ylabel('Proportion of exome covered in total')
plt.savefig('exomecov.png',bbox_inches='tight')
