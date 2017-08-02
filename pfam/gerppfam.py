from bw import BigWig
import numpy as np
import tabix 
from operator import itemgetter

def score_average(pfams,ccrs):
    varprev = None; varprevscore = 0.0
    for variant in intersect(variants, ccr).split("\n"):
        fields = variant.strip().split("\t")
        var = [fields[x] for x in [0,1,3,4]]
        overlap = int(fields[-1]) + 1 # until aaron fixes bedtools off-by-one bug, need to add +1
        ccrscore = float(fields[-2])
        if varprev is None:
            pass # will be dealt with at next pass through loop
        elif varprev == var:  
            varprevscore+=overlap
        elif varprev != var:
            if varprevscore == 0.0:
                if average:
                    varprevscore += ccrprev * overprev / float(len(varprev[2]))
            print "%.3f" % varprevscore
            varprevscore = 0.0
        varprev = var; overprev = overlap; ccrprev = ccrscore

    if varprevscore == 0.0:
        if average:
            varprevscore += ccrprev * overprev / float(len(varprev[2]))
    print "%.3f" % varprevscore

def read_ccrs(ccr, pfam):
    chrom = pfam[0]
    if "chr" in chrom:
        chrom=chrom.split("chr")[1]
    start=str(pfam[1]); end=str(pfam[2]); fam=pfam[3]
    region=chrom+":"+start+"-"+end
    try:
        for r in ccr.querys(region):
            pfam
    except:
        pass  #just because chrom isn't present in a file, pytabix will throw a "TabixError"

    return 0

def read_gerp(gerp, pfam):
    chrom = pfam[0]
    if "chr" not in chrom:
        chrom="chr"+chrom
    start=pfam[1]; end=pfam[2]; fam=pfam[3]
    #print fam, np.nanmean(np.frombuffer(gerp.values(chrom, start, end), dtype='f'))

    return 0

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

#pfampath='pfam.hg19.bed' # pfam doms incl. introns
pfampath='pfam.browser.bed' # pfam doms by exon
gerppath='/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
ccrpath='/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz'

gerp = BigWig(gerppath)
pfams=read_pfam(pfampath)
ccr = tabix.open(ccrpath)
for pfam in pfams:
    read_ccrs(ccr, pfam)
    read_gerp(gerp, pfam)
