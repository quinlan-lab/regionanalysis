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
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

miszccr = open(sys.argv[1], 'r') 
rvisccr = open(sys.argv[2], 'r')
pliccr = open(sys.argv[3], 'r')

miszlist, rvislist, plilist = [], [], []
misz0list, rvis0list, pli0list = [], [], []
for line in miszccr:
    fields = line.strip().split("\t")
    gene, val = fields[3], float(fields[4])
    ccrgene = fields[8]
    if fields[-1] == "0":
        misz0list.append((gene, val))
        continue
    if gene == ccrgene:
        miszlist.append((gene, val)) 
for line in rvisccr:
    fields = line.strip().split("\t")
    gene, val = fields[3], float(fields[4])
    ccrgene = fields[8]
    if fields[-1] == "0":
        rvis0list.append((gene, val))
        continue
    if gene == ccrgene:
        rvislist.append((gene, val)) 
for line in pliccr:
    fields = line.strip().split("\t")
    gene, val = fields[3], float(fields[4])
    ccrgene = fields[8]
    if fields[-1] == "0":
        pli0list.append((gene, val))
        continue
    if gene == ccrgene:
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
for entry in misz0list:
    misz.append(entry[-1])
    miszcounts.append(0)
for entry in rvis0list:
    rvis.append(entry[-1])
    rviscounts.append(0)
for entry in pli0list:
    pli.append(entry[-1])
    plicounts.append(0)

mpr, mpval = pearsonr(misz, miszcounts)
rpr, rpval = pearsonr(rvis, rviscounts)
ppr, ppval = pearsonr(pli, plicounts)
print mpr, rpr, ppr

# because we arbitrarily cutoff at the 30 count...
miszcounts = [30 if i > 30 else i for i in miszcounts]
rviscounts = [30 if i > 30 else i for i in rviscounts]
plicounts = [30 if i > 30 else i for i in plicounts]

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24,7))
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g1=ax1.hexbin(misz, miszcounts, cmap=cm, bins='log', alpha=0.5, mincnt=1, gridsize=40)
ax1.set_ylabel('Count of CCRs >= 95th percentile', fontsize=24)
ax1.set_xlabel('Missense Z constraint', fontsize=24)
ax1.set_ylim(-1,30) # because we arbitrarily cutoff at the 30 count
ax1.text(10, 25, 'r= ' + '{:.3f}'.format(mpr), fontsize=16) # because we arbitrarily cutoff at the 30 count
#ax1.text(10, 60, 'r= ' + '{:.3f}'.format(mpr), fontsize=16)
for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    label.set_fontsize(21)
g2=ax2.hexbin(rvis, rviscounts, cmap=cm, bins='log', alpha=0.5, mincnt=1, gridsize=40)
ax2.set_ylabel('Count of CCRs >= 95th percentile', fontsize=24)
ax2.set_xlabel('RVIS', fontsize = 24)
ax2.text(80, 25, 'r= ' + '{:.3f}'.format(rpr), fontsize=16) # because we arbitrarily cutoff at the 30 count
#ax2.text(80, 60, 'r= ' + '{:.3f}'.format(rpr), fontsize=16)
ax2.set_ylim(-1,30) # because we arbitrarily cutoff at the 30 count
ax2.set_xlim(100,0) # because we arbitrarily cutoff at the 30 count
for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
    label.set_fontsize(21)
g3=ax3.hexbin(pli, plicounts, cmap=cm, bins='log', alpha=0.5, mincnt=1, gridsize=40)
ax3.set_ylabel('Count of CCRs >= 95th percentile', fontsize=24)
ax3.set_xlabel('pLI', fontsize=24)
ax3.text(0.8, 25, 'r= ' + '{:.3f}'.format(ppr), fontsize=16) # because we arbitrarily cutoff at the 30 count
#ax3.text(0.8, 60, 'r= ' + '{:.3f}'.format(ppr), fontsize=16)
ax3.set_ylim(-1,30) # because we arbitrarily cutoff at the 30 count
for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
    label.set_fontsize(21)
def y_format(x,y):
    return '{:.0f}'.format(round(10**x)) # more round number binning
counts,edges=np.histogram(g1.get_array(),bins=8)
cbar1 = fig.colorbar(g1, ax=ax1, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar1.set_label('Number of genes', rotation=270, labelpad=20, fontsize=18)
cbar1.ax.tick_params(labelsize=18) 
counts,edges=np.histogram(g2.get_array(),bins=8)
cbar2 = fig.colorbar(g2, ax=ax2, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar2.set_label('Number of genes', rotation=270, labelpad=20, fontsize=18)
cbar2.ax.tick_params(labelsize=18) 
counts,edges=np.histogram(g3.get_array(),bins=8)
cbar3 = fig.colorbar(g3, ax=ax3, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar3.set_label('Number of genes', rotation=270, labelpad=20, fontsize=18)
cbar3.ax.tick_params(labelsize=18) 
plt.tight_layout()
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/conscorr.pdf', bbox_inches='tight')

# host = host_subplot(111, axes_class=AA.Axes)
# plt.subplots_adjust(bottom=0.25)

# par1 = host.twiny()
# par2 = host.twiny()

# offset = 40
# new_fixed_axis1 = par2.get_grid_helper().new_fixed_axis
# new_fixed_axis2 = par2.get_grid_helper().new_fixed_axis
# par1.axis["bottom"] = new_fixed_axis1(loc="bottom",
                                    # axes=par1,
                                    # offset=(0, -offset))
# par2.axis["bottom"] = new_fixed_axis2(loc="bottom",
                                    # axes=par2,
                                    # offset=(0, -offset-40))

# par1.axis["top"].toggle(all=False)
# par2.axis["top"].toggle(all=False)

# #host.set_xlim(0, 2)
#host.set_ylim(0, 2)

# host.set_xlabel("Missense Z Constraint")
# host.set_ylabel("Count of CCRs >= 95")
# par1.set_xlabel("RVIS")
# par2.set_xlabel("pLI")

# p1 = host.scatter(misz, miszcounts, s=5, alpha=0.5, label="Missense Z Constraint, r=" + '{:.3f}'.format(mpr), c='b')
# p2 = par1.scatter(rvis, rviscounts, s=5, alpha=0.5, label="RVIS, r=" + '{:.3f}'.format(rpr), c='g')
# p3 = par2.scatter(pli, plicounts, s=5, alpha=0.5, label="pLI, r=" + '{:.3f}'.format(ppr), c='r')

# host.set_xlim(-6.5, 14)
# par1.set_xlim(0, 100)
# par2.set_xlim(0, 1)

# host.legend()

#host.axis["top"].label.set_color(p1.get_color())
#par1.axis["bottom"].label.set_color(p2.get_color())
#par2.axis["bottom"].label.set_color(p3.get_color())

# plt.draw()
# plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/conscorr.pdf', bbox_inches='tight')
