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

plt.hist(lens,500,facecolor='green')
plt.savefig('funchist.png',bbox_inches='tight')

totlen=0.0; ordered=[]; exomecov=[]; wgtdpct=[]
it = ts.reader(sys.argv[2]) #$HOME/analysis/exacresiduals/results/2016_12_10/weightedresiduals.txt (sorted by weighted_pct)
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
plt.savefig('exomecov.png',bbox_inches='tight')
