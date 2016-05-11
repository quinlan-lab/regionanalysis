import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import sys
from matplotlib import pyplot as plt
import numpy as np

f = open('residualgist/residuals.txt','r')
f.readline()
cpg=[];lens=[]
for line in f:
    fields = line.strip().split("\t")
    lens.append(int(fields[2])-int(fields[1])+1)
    cpg.append(float(fields[8]))

x = np.array(lens); y = np.array(cpg)
mask = (x >= 0) & (x <= 3)
x, y = x[mask], y[mask]
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y)
plt.hexbin(x, y, cmap=plt.cm.YlOrRd_r)
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel("length")
plt.ylabel("CpG")
plt.title("Hexagon binning")
cb = plt.colorbar()
cb.set_label('counts')
plt.savefig('s.png',bbox_inches='tight')
#sns.jointplot(x, y, kind="hex", joint_kws={'gridsize':200})
#sns.plt.savefig('s.png',bbox_inches='tight')
#g = sns.JointGrid(x, y)
#g.ax_marg_x.hist(x, bins=np.arange(xmin, xmax, 100))
#g.ax_marg_y.hist(u, bins=np.arange(ymin, ymax), orientation="horizontal")
#g.plot_joint(plt.hexbin, gridsize=25, extent=[xmin, xmax, ymin, ymax], cmap="Blues")
#g.fig.savefig("s.png", bbox_inches="tight")
