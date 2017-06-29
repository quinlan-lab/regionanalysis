import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import subprocess as sp
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genes", help = "gene bed file")
parser.add_argument("-p", "--pathogenic", help = "pathogenic variant file")
parser.add_argument("-b", "--benign", help = "benign variant file")
parser.add_argument("-n", "--name", help = "ad gene distribution, for example")
#parser.set_defaults(genes = '/uufs/chpc.utah.edu/common/home/u1021864/serial/analysis/ad.bed')
#parser.set_defaults(pathogenic = '/scratch/ucgd/lustre/u1021864/serial/clinvar-patho-exac.vcf.gz')
#parser.set_defaults(benign = '/scratch/ucgd/lustre/u1021864/serial/gnomad-benign-exac.vcf.gz')
#parser.set_defaults(name = 'AD gene distribution')
args = parser.parse_args()
genes = args.genes
patho = args.pathogenic
benign = args.benign
name = args.name

def intersect(genes, variants, wo = False):
    def killproc(p):
        try:
            p.kill()
        except OSError:
            pass
    l = ['bedtools', 'intersect', '-a', genes, '-b', variants, '-sorted']
    if wo:
        l.append('-wo')
    p1 = sp.Popen(l, stdout = sp.PIPE)
    output,error = p1.communicate()
    killproc(p1)
    return output.strip()

Y = []; Y2 = [] #Y = benign, Y2 = patho
prevgene = []; ct = 0
for line in intersect(genes, patho).split("\n"):
    currgene = line.split('\t')[3]
    if prevgene is None:
        pass
    elif prevgene!=currgene:
        Y2.append(ct)
        ct = 0
    prevgene = currgene; ct += 1

for line in intersect(genes, benign).split("\n"):
    currgene = line.split('\t')[3]
    if prevgene is None:
        pass
    elif prevgene!=currgene:
        Y.append(ct)
        ct = 0
    prevgene = currgene; ct += 1

matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure(figsize=(5,5)) # adjust figsize to change shape of plot and proportions
fig.subplots_adjust(hspace=.05, wspace=.05, left=.01, bottom=.01)
ax = fig.add_subplot(1,1,1)
import seaborn
mi = np.min(Y+Y2)
ma = np.max(Y+Y2)
rng = (mi,ma)
#rng = (0,10)
p,p_edges=np.histogram(Y2, bins=40, range=(10,ma)) #bins=10
#p.append(p[-1])
width_p = (p_edges[1]-p_edges[0])
b,b_edges=np.histogram(Y, bins=40, range=(10,ma)) #bins=10
#b.append(b[-1])
width_b = (b_edges[1]-b_edges[0])
ax.bar(p_edges[:-1], p, width = width_p, color = 'orange', label = 'pathogenic', alpha = 0.7)
ax.bar(b_edges[:-1], b, width = width_b, color = 'navy', label = 'benign', alpha = 0.7)
#ax.plot(p_edges, p, label = 'pathogenic', color = 'crimson', ls = 'steps', lw=2, alpha=0.7)#pathogenic
#ax.plot(b_edges, b, label = 'benign', color = 'deepskyblue', ls = 'steps', lw=2, alpha=0.7)#benign
#ax.set_xlim(0,10)
ax.legend(loc='best')
ax.set_ylabel('Variant count per gene')
ax.set_xlabel('Genes')
plt.title(name)
plt.savefig('addist.pdf',bbox_inches='tight')
