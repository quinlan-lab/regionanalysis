from bw import BigWig
import numpy as np
import tabix 
from operator import itemgetter
import subprocess as sp

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
    famprev = None
    lengths=[]; scores=[]; ccrdict={}
    for isxn in intersect(pfams, ccr).split("\n"):
        fields = isxn.strip().split("\t")
        fam = fields[3]
        overlap = float(fields[-1]) # need to add ccr*overlap/total overlap?`
        ccrscore = float(fields[-2])
        if famprev is None:
            pass # will be dealt with at next pass through loop
        elif famprev == fam:  
            lengths.append(overprev)
            scores.append(ccrprev)
        elif famprev != fam:
            lengths.append(overprev)
            scores.append(ccrprev)
            famprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
            ccrdict[famprev]=famprevscore
            lengths=[]; scores=[]
        famprev = fam; overprev = overlap; ccrprev = ccrscore

    famprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
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
    return fam, np.nanmean(gerps), len(gerps)

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
    gerpdict={}; lengths=[]; scores=[]; famprev = None
    for pfam in pfams:
        pfam, gerpscore, overlap = read_gerp(gerp, pfam) # _ = pfam, redundant variable
        if famprev is None:
            pass # will be dealt with at next pass through loop
        elif famprev == pfam:  
            lengths.append(overprev)
            scores.append(gerpprev)
        elif famprev != pfam:
            lengths.append(overprev)
            scores.append(gerpprev)
            famprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
            gerpdict[famprev]=famprevscore
            lengths=[]; scores=[]
        famprev = pfam; overprev = overlap; gerpprev = gerpscore
    
    lengths.append(overprev) 
    scores.append(gerpprev)
    famprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
    gerpdict[famprev]=famprevscore
    return gerpdict

#pfampath='pfam.hg19.bed' # pfam doms incl. introns
pfampath='/uufs/chpc.utah.edu/common/home/u1021864/analysis/pfam/pfam.exonic.bed' # pfam doms by exon; sorted by pfam name, in pfamconvert.sh
gerppath='/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
ccrpath='/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz'

gerp = BigWig(gerppath)
pfams = read_pfam(pfampath)
#ccr = tabix.open(ccrpath)

ccrs=score_average(pfampath, ccrpath)
gerps=score_gerp(pfams, gerp)
    
for i in ccrs:
    print i, ccrs[i]
for i in gerps:
    print i, gerps[i]
#for pfam in ccrs:
