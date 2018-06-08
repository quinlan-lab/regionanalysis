import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
from collections import defaultdict
from operator import itemgetter
from itertools import groupby
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from scipy.stats import pearsonr

miszccr = open(sys.argv[1], 'r') 
rvisccr = open(sys.argv[2], 'r')
pliccr = open(sys.argv[3], 'r')


miszlist, rvislist, plilist = [], [], []
for line in miszccr:
    fields = line.strip().split("\t")
    gene, val = fields[-2], float(fields[-1])
    miszlist.append((gene, val)) 
for line in rvisccr:
    fields = line.strip().split("\t")
    gene, val = fields[-2], float(fields[-1])
    rvislist.append((gene, val)) 
for line in pliccr:
    fields = line.strip().split("\t")
    gene, val = fields[-2], float(fields[-1])
    plilist.append((gene, val))

sorter = itemgetter(0)
grouper = itemgetter(0)
misz, rvis, pli, miszcounts, rviscounts, plicounts = [], [], [], [], [], []
for key, grp in groupby(sorted(miszlist, key = sorter), grouper): #sort by each gene value (whole line)
    grp=list(grp)
    misz.append(grp[0][-1])
    miszcounts.append(len(grp))
for key, grp in groupby(sorted(rvislist, key = sorter), grouper): #sort by each gene value (whole line)
    grp=list(grp)
    rvis.append(grp[0][-1])
    rviscounts.append(len(grp))
for key, grp in groupby(sorted(plilist, key = sorter), grouper): #sort by each gene value (whole line)
    grp=list(grp)
    pli.append(grp[0][-1])
    plicounts.append(len(grp))

host = host_subplot(111, axes_class=AA.Axes)
plt.subplots_adjust(bottom=0.25)

par1 = host.twiny()
par2 = host.twiny()

offset = 40
new_fixed_axis1 = par2.get_grid_helper().new_fixed_axis
new_fixed_axis2 = par2.get_grid_helper().new_fixed_axis
par1.axis["bottom"] = new_fixed_axis1(loc="bottom",
                                    axes=par1,
                                    offset=(0, -offset))
par2.axis["bottom"] = new_fixed_axis2(loc="bottom",
                                    axes=par2,
                                    offset=(0, -offset-40))

par1.axis["top"].toggle(all=False)
par2.axis["top"].toggle(all=False)

#host.set_xlim(0, 2)
#host.set_ylim(0, 2)

host.set_xlabel("Missense Z Constraint")
host.set_ylabel("Count of CCRs >= 95")
par1.set_xlabel("RVIS")
par2.set_xlabel("pLI")

mpr, mpval = pearsonr(misz, miszcounts)
rpr, rpval = pearsonr(rvis, rviscounts)
ppr, ppval = pearsonr(pli, plicounts)

p1 = host.scatter(misz, miszcounts, s=5, alpha=0.5, label="Missense Z Constraint, r=" + '{:.3f}'.format(mpr), c='b')
p2 = par1.scatter(rvis, rviscounts, s=5, alpha=0.5, label="RVIS, r=" + '{:.3f}'.format(rpr), c='g')
p3 = par2.scatter(pli, plicounts, s=5, alpha=0.5, label="pLI, r=" + '{:.3f}'.format(ppr), c='r')

host.set_xlim(-6.5, 14)
par1.set_xlim(0, 100)
par2.set_xlim(0, 1)

host.legend()

#host.axis["top"].label.set_color(p1.get_color())
#par1.axis["bottom"].label.set_color(p2.get_color())
#par2.axis["bottom"].label.set_color(p3.get_color())

plt.draw()
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/conscorr.pdf', bbox_inches='tight')
