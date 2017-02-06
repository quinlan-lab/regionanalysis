import sys
import numpy as np
import toolshed as ts
import pyBigWig

def read_cadd(path='/scratch/ucgd/lustre/u1021864/serial/whole_genome_SNVs.tsv.gz'):
    var=None;vals=[];cadd=[]
    for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen("| tabix " + path + " {region}".format(region=region)) if x[0] != "#"):
        if var==None or var==toks[1]:
            vals.append(float(toks[5]))
        elif var!=toks[1] and var!=None:
            cadd.append(np.mean(vals))
            vals=[] 
        var=toks[1]
    return np.mean(cadd) 

def readccrs(path, gerp, phast, cadd):
    for i, d in enumerate(ts.reader(path, header="ordered")):
        d['gerp'] = ",".join(map(str,gerp.values("chr" + d['chrom'], int(d['start']), int(d['end']))))
        d['phast'] = ",".join(map(str,phast.values("chr" + d['chrom'], int(d['start']), int(d['end']))))
        region=d['chrom']+":"+d['start']+"-"+d['end']
        var=None;vals=[];caddvals=[]
        for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen("| tabix " + cadd + " {region}".format(region=region)) if x[0] != "#"):
            if var==None or var==toks[1]:
                vals.append(float(toks[5]))
            elif var!=toks[1] and var!=None:
                caddvals.append(np.mean(vals))
                vals=[]
            var=toks[1]
        d['cadd'] = ",".join(map(str,caddvals))
        if i == 0:
            print "\t".join(d.keys())
        print "\t".join(map(str, d.values()))

path = '/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'
gerp = pyBigWig.open(path)
path = '/scratch/ucgd/lustre/u1021864/serial/hg19.100way.phastCons.bw'
phast = pyBigWig.open(path)
path = '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/results/2016_12_10/weightedresiduals.txt'
cadd = '/scratch/ucgd/lustre/u1021864/serial/whole_genome_SNVs.tsv.gz'
readccrs(path,gerp,phast,cadd)
