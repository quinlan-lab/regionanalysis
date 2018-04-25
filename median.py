import numpy as np
import gzip
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
import toolshed as ts
import sys

plotout=sys.argv[3]

l95, l99 = [], []
with gzip.open(sys.argv[1], 'rb') as f:
    for line in f:
        fields = line.strip().split("\t")
        ranges = fields[6]; ccr = float(fields[-1])
        length = sum([float(i.split("-")[1]) - float(i.split("-")[0]) for i in ranges.split(",")])
        if ccr >= 95:
            l95.append(length)
        if ccr >= 99:
            l99.append(length)

lengths = []
for i, d in enumerate(ts.reader(sys.argv[2], header="ordered")):
    if d['varflag'] == "VARTRUE":
        continue
    length=sum([float(i.split("-")[1]) - float(i.split("-")[0]) for i in d['ranges'].split(",")])
    if length > 4000:
        print d['ranges']
    lengths.append(length)

matplotlib.rcParams.update({'font.size': 10})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots()
step=100
ma = max(lengths)
p,p_edges=np.histogram(lengths, bins=np.linspace(0, ma, step), range=(0,ma)) #bins=10
width_p = [p_edges[i+1] - p_edges[i] for i in range(0, len(p_edges)-1)]
ax.bar(p_edges[:-1], p, width = width_p, color = (161/255.0,218/255.0,215/255.0), edgecolor = (96/255.0, 133/255.0, 131/255.0))
ax.set_ylabel("Count")
ax.set_yscale("log")
ax.set_xlabel("Exonic distance (bp)")
ax.set_title("Average variant distance, without exclusion or coverage filtering")
ax.axvline(np.mean(lengths),label="Variant distance mean = " + '%.2f' % np.mean(lengths), color='r', ls='dashed', lw=1)
ax.axvline(np.median(l95),label=">=95% CCR percentile median = " + '%.2f' % np.median(l95), color='b', ls='dashed', lw=1)
ax.axvline(np.median(l99),label=">=99% CCR percentile median = " + '%.2f' % np.median(l99), color='k', ls='dashed', lw=1)
ax.legend(loc='upper right')
sns.despine()
plt.savefig(plotout, bbox_inches='tight')

# zoomed plot
matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
ma = 100
step = 20 
p,p_edges=np.histogram(lengths, bins=np.linspace(0, ma, step), range=(0,ma)) #bins=10
width_p = [p_edges[i+1] - p_edges[i] for i in range(0, len(p_edges)-1)]
ax.bar(p_edges[:-1], p, width = width_p, color = (161/255.0,218/255.0,215/255.0), edgecolor = (96/255.0, 133/255.0, 131/255.0))
ax.set_ylabel("Count")
ax.set_yscale("log")
ax.set_xlabel("Exonic distance (bp)")
ax.set_title("Average variant distance, without exclusion or coverage filtering")
ax.axvline(np.mean(lengths),label="Variant distance mean = " + '%.2f' % np.mean(lengths), color='r', ls='dashed')
ax.axvline(np.median(l95),label=">=95% CCR percentile median = " + '%.2f' % np.median(l95), color='b', ls='dashed')
ax.axvline(np.median(l99),label=">=99% CCR percentile median = " + '%.2f' % np.median(l99), color='k', ls='dashed')
ax.legend(loc='upper right')
sns.despine()
outsplit=plotout.rpartition(".")
outfile="".join(outsplit[0:-2])+'zoomed'+"".join(outsplit[-2:])
plt.savefig(outfile,bbox_inches='tight')
