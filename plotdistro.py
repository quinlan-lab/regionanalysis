import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import subprocess as sp

#plots the distribution of residuals in total, particularly top vs mid

f = open('regions/topresid.txt','r')

topresid,midresid,totresid = [],[],[]

for line in f:
    fields = line.strip().split("\t")
    topresid.append(float(fields[7]))
f.close()
f = open('regions/midresid.txt','r')
for line in f:
    fields = line.strip().split("\t")
    midresid.append(float(fields[7]))
f.close()
f = open('residualgist/residuals.txt', 'r')
f.readline()
for line in f:
    fields = line.strip().split("\t")
    totresid.append(float(fields[7]))
f.close()

fig = plt.figure()
ax = fig.add_subplot(111)
mi = np.min(totresid)
ma = np.max(totresid)
rng = (mi,ma)
ax.hist(topresid,bins=40,range=rng,alpha=0.3,color='r',label='topresid')
ax.hist(midresid,bins=40,range=rng,alpha=0.3,color='g',label='midresid')
ax.hist(totresid,bins=40,range=rng,alpha=0.3,color='b',label='totresid')
ax.set_xlabel("Residuals")
ax.set_ylabel("Frequency")
plt.title("Residual Comparison")
ax.legend(loc='upper right')
plt.savefig("plots/residuals.png")
plt.close()

#plotting the distribution of residuals for pathogenic vs benign

def intersect(hotspot):
    def killproc(p):
        try:
            p.kill()
        except OSError:
            pass
    p1 = sp.Popen("sed '1d' residualgist/residuals.txt", shell = True, stdout = sp.PIPE)
    p2 = sp.Popen(['bedtools', 'intersect', '-a', 'stdin', '-b', hotspot], stdin = p1.stdout, stdout = sp.PIPE)
    output,error = p2.communicate()
    killproc(p1); killproc(p2)
    return output.strip()

benign = intersect('benign.bed')
patho = intersect('patho.bed')
bentile,pathtile = [],[]
for line in benign.split("\n"):
    fields = line.strip().split("\t")
    bentile.append(float(fields[8]))
for line in patho.split("\n"):
    fields = line.strip().split("\t")
    pathtile.append(float(fields[8]))

fig = plt.figure()
ax = fig.add_subplot(111)
mi = np.min(bentile+pathtile)
ma = np.max(bentile+pathtile)
rng = (mi,ma)
b,b_edges=np.histogram(bentile, bins=40, range=rng)
p,p_edges=np.histogram(pathtile, bins=40, range=rng)
b=map(lambda x: float(x)/sum(b), b)
p=map(lambda x: float(x)/sum(p), p)
width_b = (b_edges[1]-b_edges[0])
ax.bar(b_edges[:-1], b, width = width_b, color = 'b', label = 'benign', alpha = 0.3)
width_p = (p_edges[1]-p_edges[0])
ax.bar(p_edges[:-1], p, width = width_p, color = 'r', label = 'pathogenic', alpha = 0.3)
for i,j in zip(p, p_edges):
    print i,j
ax.set_xlabel("Percentile")
ax.set_ylabel("Frequency")
plt.title("Percentile Pathogenicity Comparison")
ax.legend(loc='upper right')
plt.savefig("plots/distribution.png")
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
#bentile = [i for i in bentile if i >= 90]
#pathtile = [i for i in pathtile if i >= 90]
mi = np.min(bentile+pathtile)
ma = np.max(bentile+pathtile)
rng = (mi,ma)
b,b_edges=np.histogram(bentile, bins=400, range=rng)
p,p_edges=np.histogram(pathtile, bins=400, range=rng)
b=map(lambda x: float(x)/sum(b), b)
p=map(lambda x: float(x)/sum(p), p)
width_b = (b_edges[1]-b_edges[0])
ax.bar(b_edges[:-1], b, width = width_b, color = 'b', label = 'benign', alpha = 0.3)
width_p = (p_edges[1]-p_edges[0])
ax.bar(p_edges[:-1], p, width = width_p, color = 'r', label = 'pathogenic', alpha = 0.3)
ax.set_xlim(90, 100)
ax.set_xlabel("Percentile")
ax.set_ylabel("Frequency")
plt.title("Percentile Pathogenicity Comparison: Top 10%")
ax.legend(loc='upper right')
plt.tight_layout()
plt.savefig("plots/top10distribution.png")
plt.close()
