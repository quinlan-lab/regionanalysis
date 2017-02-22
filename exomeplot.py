import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import font_manager
import seaborn as sns
import numpy as np
import toolshed as ts

lens=[]
it = ts.reader('$HOME/analysis/exacresiduals/results/2016_12_10/unfilteredresiduals.txt')
iterable = (i for i in it)
for region in iterable:
    lens.append(int(region['end'])-int(region['start']))

fig=plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
mi=0; ma = max(lens); rng = (mi,ma)
p,p_edges=np.histogram(lens, bins=40, range=rng) #bins=400
p=map(lambda x: float(x)/sum(p), p)
width_p = (p_edges[1]-p_edges[0])
ax1.bar(p_edges[:-1], p, width = width_p, color = 'r', log = True)
ax1.xaxis.set_ticks(p_edges)#[1::2])
#ax1.set_xlim(0,115)
labs=ax1.xaxis.set_ticklabels(p_edges, rotation=45, ha='right', fontname='Cmr10')
ax1.xaxis.set_major_formatter(FormatStrFormatter('${%.0f}$'))
ticks=ax1.get_yticks()
ax1.yaxis.set_ticklabels(['{:f}%'.format(x*100) for x in ticks], fontname='Cmr10')
ax1.axvline(np.nanmean(lens), color='black', linestyle='dashed', linewidth=2)
fig.tight_layout()
ax1.set_ylabel('Proportion of regions represented')
ax1.set_xlabel('Distance between functional variants')
plt.savefig('funchist.pdf', format='pdf', bbox_inches='tight')

totlen=0.0; ordered=[]; exomecov=[]; wgtdpct=[]
it = ts.reader('$HOME/analysis/exacresiduals/results/2016_12_10/weightedresiduals.txt')
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
