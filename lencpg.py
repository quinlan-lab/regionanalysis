import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import sys
from matplotlib import pyplot as plt
import numpy as np

f = open('exacresiduals/residuals.txt','r')
f.readline()
cpg=[];lens=[]
for line in f:
    fields = line.strip().split("\t")
    lens.append(int(fields[2])-int(fields[1])+1)
    cpg.append(float(fields[8]))

'''
small length range
'''
x = np.array(lens); y = np.array(cpg)
mask = (x >= 0) & (x <= 20)
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
plt.close()
'''
mid length range
'''
x = np.array(lens); y = np.array(cpg)
mask = (x >= 20) & (x <= 200)
x, y = x[mask], y[mask]
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y)
plt.hexbin(x, y, cmap=plt.cm.YlOrRd_r)
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel("length")
plt.ylabel("CpG")
plt.title("Hexagon binning")
cb = plt.colorbar()
cb.set_label('counts')
plt.savefig('m.png',bbox_inches='tight')
plt.close()

'''
larger length range
'''
x = np.array(lens); y = np.array(cpg)
mask = (x >= 200) & (x <= 2000)
x, y = x[mask], y[mask]
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y)
plt.hexbin(x, y, cmap=plt.cm.YlOrRd_r)
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel("length")
plt.ylabel("CpG")
plt.title("Hexagon binning")
cb = plt.colorbar()
cb.set_label('counts')
plt.savefig('l.png',bbox_inches='tight')
plt.close()
