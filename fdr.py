import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
import seaborn as sns
sns.set_style('white')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

benign = open(sys.argv[1], 'r')
patho = open(sys.argv[2], 'r')
plotout = sys.argv[3] # pdf location

benign = np.array([float(line) for line in benign])
patho = np.array([float(line) for line in patho])

print benign.size

fdr, fpr, cutoffs = [], [], []
print "FPR"
for cutoff in (90, 95, 99):
    cutoffs.append(cutoff)
    fpr.append(benign[benign >= cutoff].size / float(len(benign)))
    print cutoff, (benign >= cutoff).sum(), benign[benign >= cutoff].size / float(len(benign))

print "FDR"
for cutoff in (90, 95, 99):
    tp = (patho >= cutoff).sum()
    fp = (benign >= cutoff).sum()
    fdr.append(fp / float(tp + fp))
    print cutoff, "tp:", tp, "fp:", fp, fp / float(tp + fp)


fdr, fpr, cutoffs = np.array(fdr), np.array(fpr), np.array(cutoffs)
bins = range(0,101,5)
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
width=1
rects1 = ax.bar(cutoffs, fpr, width, color = (161/255.0,218/255.0,215/255.0), edgecolor = (96/255.0, 133/255.0, 131/255.0), label = "FPR")
rects2 = ax.bar(cutoffs+width, fdr, width, color = (56/255.0,138/255.0,172/255.0), edgecolor = (96/255.0, 133/255.0, 131/255.0), label = "FDR")
ax.set_xticks(cutoffs+width/2)
ax.set_xticklabels(map(str,cutoffs))
ax.legend()
ax.set_xlabel('CCR Percentile Cutoff') # FDR = (benigns/benigns+pathogenics in bin)
ax.set_ylim(0,0.1)
plt.tight_layout()
plt.savefig(plotout, bbox_inches='tight')
plt.close()
