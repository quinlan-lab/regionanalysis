from bw import BigWig
import numpy as np
import tabix
import cPickle as pickle

def read_gerp(gerp, region):
    chrom = region[0]
    if "chr" not in chrom:
        chrom="chr"+chrom
    start=int(region[1]); end=int(region[2]); ranges=region[6]
    gerps=np.frombuffer(gerp.values(chrom, start, end), dtype='f')
    return ranges, np.nanmean(gerps), len(gerps)

gerppath='/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
ccrpath='/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz'

def perchrom(ccr_gerp_chrom):
    ccrpath, gerppath, chrom = ccr_gerp_chrom
    ccr = tabix.open(ccrpath)
    gerp = BigWig(gerppath)

    gerpdict={}; lengths=[]; scores=[]; rangeprev = None
    for region in ccr.querys(chrom):
        ranges, gerpscore, overlap = read_gerp(gerp, region) # _ = pfam, redundant variable
        if rangeprev is None:
            pass # will be dealt with at next pass through loop
        elif rangeprev == ranges:
            lengths.append(overprev)
            scores.append(gerpprev)
        elif rangeprev != ranges:
            lengths.append(overprev)
            scores.append(gerpprev)
            ccrprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
            pctile=float(region[-1])
            gerpdict[rangeprev]=(ccrprevscore,pctile)
            lengths=[]; scores=[]
        rangeprev = ranges; overprev = overlap; gerpprev = gerpscore; pctile=region[-1]

    lengths.append(overprev)
    scores.append(gerpprev)
    ccrprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
    gerpdict[rangeprev]=(ccrprevscore,pctile)
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
