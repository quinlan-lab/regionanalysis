from bw import BigWig
import numpy as np
import tabix 
from itertools import groupby
from operator import itemgetter
import subprocess as sp
import cPickle as pickle
import sys

def intersect(variants, regions, wo = True, sortedbool=False):
    def killproc(p):
        try:
            p.kill()
        except OSError:
            pass
    l = ['bedtools', 'intersect', '-a', variants, '-b', regions]
    if wo:
        l.append('-wo')
    if sortedbool:
        l.append('-sorted')
    p1 = sp.Popen(l, stdout = sp.PIPE)
    output,error = p1.communicate()
    killproc(p1)
    return output.strip()

def score_average(pfams,ccr):
    ccrdict={}; isxns = []
    for isxn in intersect(pfams, ccr).split("\n"):
        fields = isxn.strip().split("\t")
        fam = fields[3]
        overlap = float(fields[-1]) # need to add ccr*overlap/total overlap?`
        ccrscore = float(fields[-2])
        isxns.append((fam, overlap, ccrscore))
    sorter = itemgetter(0)
    grouper = itemgetter(0)
    for key, grp in groupby(sorted(isxns, key = sorter), grouper):
        grp = list(grp)
        fam = grp[0][0]
        lengths = []; scores = []
        for i, elem in enumerate(grp):
            lengths.append(grp[i][1])
            scores.append(grp[i][-1])
        famscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
        ccrdict[fam]=famscore

    return ccrdict

def read_ccrs(ccr, pfam):
    ccrs=[]
    chrom = pfam[0]
    if "chr" in chrom:
        chrom=chrom.split("chr")[1]
    start=str(pfam[1]); end=str(pfam[2]); fam=pfam[3]
    region=chrom+":"+start+"-"+end
    try:
        for r in ccr.querys(region):
            ccrs.append((r[0],int(r[1]),int(r[2]),r[-1]))
    except:
        pass  #just because chrom isn't present in a file, pytabix will throw a "TabixError"
    
    return ccrs

def read_gerp(gerp, pfam):
    chrom = pfam[0]
    if "chr" not in chrom:
        chrom="chr"+chrom
    start=pfam[1]; end=pfam[2]; fam=pfam[3]
    gerps=np.frombuffer(gerp.values(chrom, start, end), dtype='f')
    meangerp=np.nanmean(gerps); lengerps=len(gerps)
    if not np.isnan(meangerp): 
        return fam, meangerp, lengerps

def read_pfam(pfam):
    pfam=open(pfam,'r')
    pfams=[]
    for line in pfam:
        fields=line.strip().split()
        chrom=fields[0]; start=int(fields[1]); end=int(fields[2]); fam=fields[3]
        fam=fam.split("_exon")[0]
        pfams.append((chrom,start,end,fam)) 
    pfam.close()
    pfams.sort(key=itemgetter(3)) # much faster than lambda tup: tup[3]
    return pfams

def score_gerp(pfams, gerp):
    gerpdict={}; fams = []
    for pfam in pfams:
        result = read_gerp(gerp, pfam) # _ = pfam, redundant variableA
        if result:
            pfam, gerpscore, overlap = result
            fams.append((pfam, gerpscore, overlap))
    sorter = itemgetter(0)
    grouper = itemgetter(0)
    for key, grp in groupby(sorted(fams, key = sorter), grouper):
        lengths=[]; scores=[]
        grp = list(grp)
        family = grp[0][0]
        for i, elem in enumerate(grp):
            scores.append(grp[i][1])
            lengths.append(grp[i][-1])
        famscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
        gerpdict[family]=famscore
    
    return gerpdict

#pfampath='pfam.hg19.bed' # pfam doms incl. introns
pfampath = sys.argv[1] # from '/uufs/chpc.utah.edu/common/home/u1021864/analysis/pfam/pfam.genome.bed' # sorted by pfam name, in genome space
gerppath = sys.argv[2] # '/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
ccrpath = sys.argv[3] # '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz'

gerp = BigWig(gerppath)
pfams = read_pfam(pfampath)
#ccr = tabix.open(ccrpath)

ccrs=score_average(pfampath, ccrpath)
gerps=score_gerp(pfams, gerp)
    
#for i in ccrs:
#    print i, ccrs[i]
#for i in gerps:
#    print i, gerps[i]
cscores, gscores, labels = [], [], []
for pfam in ccrs:
    cscores.append(ccrs[pfam])
    gscores.append(gerps[pfam])
    labels.append(pfam)
data=[cscores,gscores,labels]
output = open('ccrgerppfam.pkl', 'wb')
pickle.dump(data,output)
output.close()
