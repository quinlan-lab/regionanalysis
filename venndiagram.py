from quicksect import Interval, IntervalTree
import sys
from gzip import open
from collections import defaultdict
import matplotlib 
matplotlib.use('Agg')
import seaborn as sns
sns.set_style('white')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from itertools import groupby
from operator import itemgetter


# create list of tuple intervals for each ccr, pli, and missense z interval and a tree of those intervals
# search each tree with both lists

ccrs = open(sys.argv[1], "rb") #exacresiduals/gnomad10x.5syn-ccrs.bed.gz
mpcs = open(sys.argv[2], "rb") #essentials/mpc.regions.clean.sorted.bed.gz
plis = open(sys.argv[3], "rb") #$HOME/software/pathoscore/score-sets/GRCh37/pLI/pLI.bed.gz

# generate data: we want to search the trees with the lists to get the numbers in each venn diagram
ccrtree = defaultdict(lambda: IntervalTree())
ccrlist = defaultdict(list)
sorter = itemgetter(0,3,6)
grouper = itemgetter(0,3,6)
ccrtemp = []
ccrgenes = set(); pligenes = set(); mpcgenes = set()
for ccr in ccrs:
    ccr = ccr.strip().split("\t")
    if float(ccr[-1]) < 95: continue
    ccrtemp.append(ccr)
for key, grp in groupby(sorted(ccrtemp, key = sorter), grouper):
    grp = list(grp)
    chrom = grp[0][0]; gene = grp[0][3]; ranges = grp[0][6]
    r = ranges.split(",")
    ccrgenes.add(gene)
    ccrtree[chrom].add(int(r[0].split("-")[0]), int(r[-1].split("-")[-1]), gene)
    ccrlist[chrom].append((int(r[0].split("-")[0]), gene, int(r[-1].split("-")[-1])))
#sanity check: 
#print ccrtree['1'].search(2105429,2105455)

mpctree = defaultdict(lambda: IntervalTree())
mpclist = defaultdict(list)
for mpc in mpcs:
    mpc = mpc.strip().split("\t")
    if float(mpc[-1]) > 0.4: continue
    chrom = mpc[0]; start = int(mpc[1]); end = int(mpc[2]); gene = mpc[3]
    mpcgenes.add(gene)
    mpctree[chrom].add(start, end, gene)
    mpclist[chrom].append((start, gene, end))
#sanity check:
#print mpctree['1'].search(69090,70008)

plitree = defaultdict(lambda: IntervalTree())
plilist = defaultdict(list)
for pli in plis: # pli only has one score for each gene
    pli = pli.strip().split("\t")
    if float(pli[-1]) < 0.9: continue
    chrom = pli[0]; start = int(pli[1]); end = int(pli[2]); gene = pli[3]
    pligenes.add(gene)
    plitree[chrom].add(start, end, gene)
    plilist[chrom].append((start, gene, end))
#sanity check:
#print plitree['1'].search(9770513,9787104)

# now the search (autosomes only)
import multiprocessing as mp
p = mp.Pool(12)

# REGION BASED

pct, mct, cct, pmct, pcct, mcct, pmcct = 0, 0, 0, 0, 0, 0, 0
ct = 0
prevtree = defaultdict(lambda: IntervalTree())
for chrom in map(str,range(1, 23)):
    for ccr in ccrlist[chrom]:
        ct += 1
        plisearch = plitree[chrom].search(ccr[0], ccr[-1])
        if plisearch:
            if len(plisearch)>1:
                ct += len(plisearch)-1
            for intr in plisearch: # we deliberately double count this.
                prevtree[chrom].insert(intr)
                mpcsearch = mpctree[chrom].find(intr) # we deliberately double count this.
                if mpcsearch:
                    ct += len(mpcsearch)-1
                    for intr2 in mpcsearch:
                        pmcct += 1
                else:
                    pcct += 1
        else:
            mpcsearch = mpctree[chrom].search(ccr[0], ccr[-1])
            if mpcsearch:
                if len(mpcsearch)>1:
                    ct += len(mpcsearch)-1
                for intr in mpcsearch:
                    prevtree[chrom].insert(intr)
                    mcct += 1
            else:
                cct += 1
        prevtree[chrom].add(ccr[0], ccr[-1])

# now copy above, but make sure it is not in prevtree
for chrom in map(str,range(1, 23)):
    for mpc in mpclist[chrom]:
        plisearch = plitree[chrom].search(mpc[0], mpc[-1])
        if plisearch:
            for intr in plisearch:
                prevsearch = prevtree[chrom].find(intr)
                if prevsearch: continue
                prevtree[chrom].insert(intr)
                pmct += 1
        else:
            prevsearch = prevtree[chrom].find(intr)
            if prevsearch: continue
            mct += 1
        prevtree[chrom].add(mpc[0], mpc[-1])

for chrom in map(str,range(1, 23)):
    for pli in plilist[chrom]:
        prevsearch = prevtree[chrom].search(pli[0], pli[-1])
        if prevsearch: continue
        pct += 1

print "Region/Intersection based"
print pct, pmct, pcct, pmcct, mct, mcct, cct
print ct

v = venn3_unweighted(subsets = (cct, mct, mcct, pct, pcct, pmct, pmcct), set_labels = ('CCR >= 95', 'Missense depletion <= 0.4', 'pLI >= 0.9'))
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/venn.pdf', bbox_inches='tight')

# just to see if there are regions in genes that would be unlooked at
# GENE BASED

print "Gene based"
print len(ccrgenes)
plt.clf()
# basic set theory
v = venn3_unweighted(subsets = (len(ccrgenes - mpcgenes - pligenes), len(mpcgenes - ccrgenes - pligenes), len(mpcgenes & ccrgenes - pligenes), len(pligenes - ccrgenes - mpcgenes), len(pligenes & ccrgenes - mpcgenes), len(pligenes & mpcgenes - ccrgenes), len(pligenes & mpcgenes & ccrgenes)), set_labels = ('CCR >= 95', 'Missense depletion <= 0.4', 'pLI >= 0.9'))
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/venngene.pdf', bbox_inches='tight')
