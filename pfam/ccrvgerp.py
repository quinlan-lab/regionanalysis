from bw import BigWig
import numpy as np
import tabix
import cPickle as pickle
import sys
from itertools import groupby
from operator import itemgetter

def read_gerp(gerp, region):
    chrom = region[0]
    if "chr" not in chrom:
        chrom="chr"+chrom
    start=int(region[1]); end=int(region[2])
    gerps=np.frombuffer(gerp.values(chrom, start, end), dtype='f')
    return np.nanmean(gerps), len(gerps)

gerppath = sys.argv[1] #'/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
ccrpath = sys.argv[2] #'/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz'

def perchrom(ccr_gerp_chrom):
    ccrpath, gerppath, chrom = ccr_gerp_chrom
    ccr = tabix.open(ccrpath)
    gerp = BigWig(gerppath)

    gerpdict={}; gerps=[]
    for region in ccr.querys(chrom):
        gene=region[3]; ranges=region[6]; pctile=float(region[-1])
        gerpscore, overlap = read_gerp(gerp, region) # _ = pfam, redundant variable
        gerps.append((gerpscore, overlap, ranges, gene, pctile))
    sorter = itemgetter(2,3)
    grouper = itemgetter(2,3)
    for key, grp in groupby(sorted(gerps, key = sorter), grouper):
        lengths = []; scores = []
        grp = list(grp)
        ranges = grp[0][2]; gene = grp[0][3]; pctile = grp[0][-1]
        for i, elem in enumerate(grp):
            scores.append(grp[i][0])
            lengths.append(grp[i][1])
        gerpscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
        gerpdict[key]=(gerpscore,pctile,gene,sum(lengths),ranges,chrom)

    return gerpdict

import multiprocessing as mp
p = mp.Pool(12)

data=[]
for outs in p.imap_unordered(perchrom, ((ccrpath, gerppath, str(chrom)) for chrom in range(1, 23))):
   for d in outs:
       data.append(outs[d])
output = open('ccrgerp.pkl', 'wb')
pickle.dump(data,output)
output.close()
# for chrom in range(1,23):
    # for outs in perchrom((ccrpath, gerppath, str(chrom))):
        # print chrom
        # for d in outs:
            # print("\t".join(map(str, (d[k] for k in keys))))
