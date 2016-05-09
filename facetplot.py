import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import sys
from matplotlib import gridspec
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np

f = open('nonpasstop.txt','r')
f.readline()
cpg=[];gerp=[];nonpass=[]
for line in f:
    fields = line.strip().split("\t")
    cpg.append(float(fields[8]))
    gerp.append(float(fields[9]))
    nonpass.append(float(fields[13]))

cmean=np.mean(cpg);cstd=np.std(cpg);gmean=np.mean(gerp);gstd=np.std(gerp);nmean=np.mean(nonpass);nstd=np.std(nonpass)
fig = plt.figure()
gs = gridspec.GridSpec(3, 1, height_ratios=[3, 3, 3])
ax0 = fig.add_subplot(gs[0])
ax0.plot(cpg,gerp,'.',zorder=-1)
ax0.set_xlabel('cpg')
ax0.set_ylabel('gerp')
ax0.add_patch(
    patches.Rectangle(
        (cmean-cstd, gmean-gstd),   # (x,y)
        2*cstd,         # width
        2*gstd,         # height
        fill=False,     # remove background
        ec='red',       # color of edge
        linewidth=2.0
    )
)
ax1 = fig.add_subplot(gs[1])
ax1.plot(cpg,nonpass,'.',zorder=-1)
ax1.set_xlabel('cpg')
ax1.set_ylabel('non-pass density')
ax1.add_patch(
    patches.Rectangle(
        (cmean-cstd, nmean-nstd),   # (x,y)
        2*cstd,         # width
        2*nstd,         # height
        fill=False,     # remove background
        ec='red',       # color of edge
        linewidth=2.0
    )
)
ax2 = fig.add_subplot(gs[2])
ax2.plot(gerp,nonpass,'.',zorder=-1)
ax2.set_xlabel('gerp')
ax2.set_ylabel('non-pass density')
ax2.add_patch(
    patches.Rectangle(
        (gmean-gstd, nmean-nstd),   # (x,y)
        2*gstd,         # width
        2*nstd,         # height
        fill=False,     # remove background
        ec='red',       # color of edge
        linewidth=2.0
    )
)
plt.tight_layout()
plt.savefig('facets.png',bbox_inches='tight')
#sns.jointplot(x, y, kind="hex", joint_kws={'gridsize':200})
#sns.plt.savefig('s.png',bbox_inches='tight')
#g = sns.JointGrid(x, y)
#g.ax_marg_x.hist(x, bins=np.arange(xmin, xmax, 100))
#g.ax_marg_y.hist(u, bins=np.arange(ymin, ymax), orientation="horizontal")
#g.plot_joint(plt.hexbin, gridsize=25, extent=[xmin, xmax, ymin, ymax], cmap="Blues")
#g.fig.savefig("s.png", bbox_inches="tight")
