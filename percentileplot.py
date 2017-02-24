import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import FormatStrFormatter

pct=[]
for region in open(sys.argv[1],'r'):
    region=region.strip().split('\t')
    try:
        pct.append(float(region[11]))
    except ValueError:
        print region
fig=plt.figure()
ax1=fig.add_subplot(1, 1, 1)
mi=0; ma = max(pct); rng = (mi,ma)
p,p_edges=np.histogram(pct, bins=40, range=rng) #bins=400
p=map(lambda x: float(x)/sum(p), p)
width_p = (p_edges[1]-p_edges[0])
ax1.bar(p_edges[:-1], p, width = width_p, color = 'r', log = True)
ax1.xaxis.set_ticks(p_edges)#[1::2])
#ax1.set_xlim(0,115)
labs=ax1.xaxis.set_ticklabels(p_edges, rotation=45, ha='right', fontname='Cmr10')
ax1.xaxis.set_major_formatter(FormatStrFormatter('${%.0f}$'))
ticks=ax1.get_yticks()
ax1.yaxis.set_ticklabels(['{:f}%'.format(x*100) for x in ticks], fontname='Cmr10')
ax1.axvline(np.nanmean(pct), color='black', linestyle='dashed', linewidth=2)
ax1.set_xlabel('CCR Percentile')
ax1.set_ylabel('Proportion of Regions in LQTS genes')
fig.tight_layout()
plt.savefig('lqts.pdf', format='pdf', bbox_inches='tight')
