from cyvcf2 import VCF
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
import sys

control = sys.argv[1] # samocha benign set
plotout = sys.argv[2]

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

fig, axes = plt.subplots(3, sharey=True)

for row, cutoff in enumerate((99, 95, 90)):
    rel_posns = []
    n = 0

    for v in VCF(control):
        if v.INFO.get('CCR', 0) >= cutoff:
            n += 1
            start = v.start
            exons = v.INFO['ccr_ranges']
            exons = [(int(a), int(b)) for a, b in (x.split("-") for x in exons.split(","))]
            assert all(b > a for a, b in exons)

            k = 0
            while k < len(exons) - 1 and start > exons[k][1]:
                k += 1
            # sum all the previous exons, add the start, and subtract the current exon start
            into = sum(b - a for a, b in exons[:k]) + v.start - exons[k][0]
            total = float(sum(b - a for a, b in exons))

            rel_posns.append(into / total)

            #dists.append(d)
            # divide by 2 because we take min so the most we can be is in the center.
    print n

    axes[row].hist(rel_posns, 10)
    axes[row].set_ylabel("cutoff: %d n: %d" % (cutoff, n))

axes[len(axes)-1].set_xlabel("relative position of benign variant into CCR")
plt.tight_layout()
plt.savefig(plotout, bbox_inches='tight')
plt.close() 
