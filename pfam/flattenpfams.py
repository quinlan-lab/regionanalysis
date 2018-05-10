import subprocess
import sys
from itertools import groupby
from operator import itemgetter

f=open(sys.argv[1],'r') # $DATA/pfam.bed already sorted pfam bed but for all transcripts domain definitions; cut -f -3, 10 are the important columns, split by spaces and semicolons
gbed=''
domains=[]
for line in f:
    fields=line.strip().split('\t')
    chrom=fields[0];start=int(fields[1]);end=int(fields[2]);info=fields[9];
    domain=[i.split('"') for i in info.split(";")][0][1]
    domains.append((line,chrom,start,end,domain))
sorter = itemgetter(4,1,2,3)
grouper = itemgetter(4)
for key, grp in groupby(sorted(domains, key = sorter), grouper):
    grp=list(grp)
    domain=grp[0][-1]
    for i, elem in enumerate(grp):
        chrom=grp[i][1]; start=str(grp[i][2]); end=str(grp[i][3])
        gbed+="\t".join([chrom,start,end])+"\n"
    p=subprocess.Popen(['bedtools merge'],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    output,err = p.communicate(gbed)
    for line in output.strip().split("\n"):
        print line + "\t" + domain
    gbed=''
