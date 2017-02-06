import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import toolshed as ts
import seaborn as sns
sns.set_style('whitegrid')

def read_values(path='/uufs/chpc.utah.edu/common/home/u1021864/analysis/scoredregions.bed'):
    gerp,phast,cadd,topgerp,topphast,topcadd,ccrpct=[],[],[],[],[],[],[]
    for i, region in enumerate(ts.reader(path, header="ordered")):
        g = np.mean(map(float,region['GERP'].split(",")))
        p = np.mean(map(float,region['phastCons'].split(",")))
        c = np.mean(map(float,region['CADD'].split(",")))
        pct=float(region['weighted_pct'])
        gerp.append(g);phast.append(p);cadd.append(c);ccrpct.append(pct)
        if pct>=90:
            topgerp.append(g);topphast.append(p);topcadd.append(c)
    return gerp,phast,cadd,topgerp,topphast,topcadd,ccrpct

fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows=4,ncols=1)
fig.tight_layout()
gerp,phast,cadd,topgerp,topphast,topcadd,ccrpct=read_values('test.bed')

mi,ma=np.amin(gerp),np.amax(gerp)
rng = (mi,ma)
m,m_edges=np.histogram(gerp, bins=40, range=rng)
m=map(lambda x: float(x)/sum(m), m)
width_m = (m_edges[1]-m_edges[0])
ax1.bar(m_edges[:-1], m, width = width_m, color = 'r', label = 'all', alpha = 0.5)
mi,ma=np.amin(topgerp),np.amax(topgerp)
rng = (mi,ma)
t,t_edges=np.histogram(topgerp, bins=40, range=rng)
t=map(lambda x: float(x)/sum(t), t)
width_t = (t_edges[1]-t_edges[0])
ax1.bar(t_edges[:-1], t, width = width_t, color = 'b', label = 'top 10%', alpha = 0.5)
ax1.set_xlabel('GERP');ax1.set_ylabel('Frequency')
mi,ma=np.amin(phast),np.amax(phast)
rng = (mi,ma)
m,m_edges=np.histogram(phast, bins=40, range=rng)
m=map(lambda x: float(x)/sum(m), m)
width_m = (m_edges[1]-m_edges[0])
ax2.bar(m_edges[:-1], m, width = width_m, color = 'r', label = 'all', alpha = 0.5)
mi,ma=np.amin(topphast),np.amax(topphast)
rng = (mi,ma)
t,t_edges=np.histogram(topphast, bins=40, range=rng)
t=map(lambda x: float(x)/sum(t), t)
width_t = (t_edges[1]-t_edges[0])
ax2.bar(t_edges[:-1], t, width = width_t, color = 'b', label = 'top 10%', alpha = 0.5)
ax2.set_xlabel('phastCons');ax2.set_ylabel('Frequency')
mi,ma=np.amin(cadd),np.amax(cadd)
rng = (mi,ma)
m,m_edges=np.histogram(cadd, bins=40, range=rng)
m=map(lambda x: float(x)/sum(m), m)
width_m = (m_edges[1]-m_edges[0])
ax3.bar(m_edges[:-1], m, width = width_m, color = 'r', label = 'all', alpha = 0.5)
mi,ma=np.amin(topcadd),np.amax(topcadd)
rng = (mi,ma)
t,t_edges=np.histogram(topcadd, bins=40, range=rng)
t=map(lambda x: float(x)/sum(t), t)
width_t = (t_edges[1]-t_edges[0])
ax3.bar(t_edges[:-1], t, width = width_t, color = 'b', label = 'top 10%', alpha = 0.5)
ax3.set_xlabel('CADD');ax3.set_ylabel('Frequency')
ax4.plot(ccrpct,gerp,label='ccr_vs_gerp',color='b', alpha=0.5)
ax4.plot(ccrpct,phast,label='ccr_vs_phastcons',color='r', alpha=0.5)
ax4.plot(ccrpct,cadd,label='ccr_vs_cadd',color='g', alpha=0.5)
ax4.set_xlabel('CCR Percentile');ax4.set_ylabel('Score')
plt.legend()
plt.savefig('conservation.png',bbox_inches='tight')
