import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import subprocess as sp

#plots the distribution of residuals in total, particularly top vs mid

f = open('regions/topresid.txt','r')

topresid,midresid,totresid = [],[],[]
lengths = []

for line in f:
    fields = line.strip().split("\t")
    topresid.append(float(fields[10])) # raw resids/studentized
f.close()
f = open('regions/midresid.txt','r')
for line in f:
    fields = line.strip().split("\t")
    midresid.append(float(fields[10]))
f.close()
f = open('exacresiduals/results/2016_06_16/exonicresiduals.txt', 'r')
f.readline()
for line in f:
    fields = line.strip().split("\t")
    totresid.append(float(fields[10]))
    lengths.append(int(fields[2])-int(fields[1]))
f.close()

'''
residual distribution
'''
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

'''
distribution of lengths
'''
fig = plt.figure()
ax = fig.add_subplot(111)
mi = np.min(lengths)
ma = np.max(lengths)
rng = (mi,ma)
l,l_edges=np.histogram(lengths, bins=50, range=rng)
l=map(lambda x: float(x)/sum(l), l)
width_l = (l_edges[1]-l_edges[0])
ax.bar(l_edges[:-1], l, width = width_l, color = 'b')
ax.set_xlabel("Length")
ax.set_ylabel("Frequency")
plt.title("Lengths of residuals, with introns excised")
plt.savefig("plots/lengths.png")

