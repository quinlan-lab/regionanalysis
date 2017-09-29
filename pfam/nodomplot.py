import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from bisect import bisect_left
import gzip
from collections import defaultdict
from operator import itemgetter
import numpy as np
import sys
import seaborn as sns
sns.set_style('white')

pdict = defaultdict(list)
cdict = defaultdict(list)
gdict = defaultdict(list)

ccrs=sys.argv[1] # topnonpfamccrs.bed
pfamlist=sys.argv[2] # pfam.exonic.bed.gz
flatexome=sys.argv[3] # /uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/flatexome.bed
plotout=sys.argv[4] # /uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/nodom.pdf

offset=defaultdict(int)

# creates gene-key dictionary of genomic starts and length offsets
with open(flatexome, 'r') as f:
    for line in f:
        fields = line.strip().split("\t")
        s = int(fields[1]); e = int(fields[2]); g = fields[3] # start, end, gene
        if gdict[g]:
            gdict[g].append((s, gdict[g][-1][1] + offset[g], e))
        else:
            gdict[g].append((s, 0, e))
        offset[g] = e - s

# custom bisection function for BED files
def bisect(lst, value, key=None): 
    if key is None:
        key = lambda x: x
    def bis(lo, hi=len(lst)):
        while lo < hi:
            mid = (lo + hi) // 2
            if key(lst[mid]) <= value:
                lo = mid + 1
            else:
                hi = mid
        return lo
    return bis(0)

# maps genome coordinates to gene-based flattened exon-only coordinates using gene-key dictionary of genomic starts and exon-based length offsets
def exomemap(s, e, g, gdict):
    #Es = s + gdict[g][bisect_left(gdict[g], (s, ))-1][1] - gdict[g][bisect_left(gdict[g], (s, ))-1][0] # exon-based start
    #Ee = e + gdict[g][bisect_left(gdict[g], (e, ))-1][1] - gdict[g][bisect_left(gdict[g], (e, ))-1][0] # exon-based end 
    Es = s + gdict[g][bisect(gdict[g], s, key=itemgetter(0))-1][1] - gdict[g][bisect(gdict[g], s, key=itemgetter(0))-1][0] # exon-based start
    Ee = e + gdict[g][bisect(gdict[g], e, key=itemgetter(0))-1][1] - gdict[g][bisect(gdict[g], e, key=itemgetter(0))-1][0] # exon-based end 
    return Es, Ee

with open(ccrs, 'r') as f:
    rangeprev = None
    for line in f:
        fields = line.strip().split("\t")
        r = fields[6]; g = fields[3] # ranges, gene
        if rangeprev and rangeprev + geneprev != r + g:
            se = [i.split("-") for i in rangeprev.split(",")]
            s = int(se[0][0]); e = int(se[-1][-1])
            Es, Ee = exomemap(s, e, geneprev, gdict)
            cdict[geneprev].append((Es, Ee, s, e))
        rangeprev = r; geneprev = g; lineprev = line

se = [i.split("-") for i in rangeprev.split(",")]
s = int(se[0][0]); e = int(se[-1][-1])
Es, Ee = exomemap(s, e, geneprev, gdict)
cdict[geneprev].append((Es, Ee, s, e))
 
with gzip.open(pfamlist, 'rb') as f:
    for line in f:
        fields = line.strip().split("\t")
        s = int(fields[1]); e = int(fields[2]); p = fields[3]; g = fields[4] # start, end, family/domain, gene
        Es, Ee = exomemap(s, e, g, gdict)
        pdict[g].append((Es, Ee)) 

#print "Exon Offset"
#for i in gdict:
#    print i, gdict[i]
#print "CCRs in Exon Space"
#for i in cdict:
#    print i, cdict[i]
#print "Pfams in Exon Space"
#for j in pdict:
#    print j, pdict[j]

distance, coords = [], []
for gene in cdict:
    for ccr in cdict[gene]:
        si = bisect(pdict[gene], ccr[0], key=itemgetter(0))
        ei = bisect_left(pdict[gene], (ccr[1],))
        if pdict[gene]:
            if si != 0:
                ds = ccr[0] - pdict[gene][si-1][1]
            else:
                ds = pdict[gene][0][0] - ccr[1]
        else:
            ds = ccr[0]
        if ei < len(pdict[gene]):
            de = pdict[gene][ei][0] - ccr[1]
        else:
            if pdict[gene]:
                de = ccr[0] - pdict[gene][-1][1]
                if de < 0:
                    print ccr, pdict[gene][-1], len(pdict[gene]), ei
            else:
                de = gdict[gene][-1][1] + (gdict[gene][-1][2] - gdict[gene][-1][0]) - ccr[1]
        distance.append(ds if ds < de else de)
        coords.append((ccr[2],ccr[3], gene))

# full-size plot
matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
mi, ma = min(distance), max(distance)
print "distances, ccrs"
for i, j in zip(distance, coords):
    print i, j 
step=200
p,p_edges=np.histogram(distance, bins=np.linspace(0, ma, ma/step), range=(0,ma)) #bins=10
#print p, p_edges
#width_p = (p_edges[-1]-p_edges[-2])-(p_edges[-1]-p_edges[-2])/5
width_p = [p_edges[i+1] - p_edges[i] for i in range(0, len(p_edges)-1)]
ax.bar(p_edges[:-1], p, width = width_p, color = 'orange', label = 'pathogenic', alpha = 0.7)
ax.set_ylabel("Count")
ax.set_yscale("log")
ax.set_xlabel("Exonic Distance in bp")
ax.set_title("Non-Pfam-intersecting CCR (>=95%) distance from Pfam domains")
sns.despine()
plt.savefig(plotout,bbox_inches='tight')

# zoomed plot
matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(10,10)) # adjust figsize to change shape of plot and proportions
ax = fig.add_subplot(1,1,1)
mi, ma = min(distance), 2000
#print "distances, ccrs"
for i, j in zip(distance, coords):
    print i, j
step=100
p,p_edges=np.histogram(distance, bins=np.linspace(0, ma, step), range=(0,ma)) #bins=10
#print p, p_edges
#width_p = (p_edges[-1]-p_edges[-2])-(p_edges[-1]-p_edges[-2])/5
width_p = [p_edges[i+1] - p_edges[i] for i in range(0, len(p_edges)-1)]
ax.bar(p_edges[:-1], p, width = width_p, color = 'orange', label = 'pathogenic', alpha = 0.7)
ax.set_ylabel("Count")
ax.set_yscale("log")
ax.set_xlabel("Exonic Distance in bp")
ax.set_title("Non-Pfam-intersecting CCR (>=95%) distance from Pfam domains")
sns.despine()
outsplit=plotout.rpartition(".")
outfile="".join(outsplit[0:-2])+'zoomed'+"".join(outsplit[-2:])
plt.savefig(outfile,bbox_inches='tight')
