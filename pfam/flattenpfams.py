import subprocess
import sys

f=open(sys.argv[1],'r') # $DATA/pfam.bed already sorted pfam bed but for all transcripts domain definitions; cut -f -3, 10 are the important columns, split by spaces and semicolons
prevdomain=None
gbed=''
for line in f:
    fields=line.strip().split('\t')
    chrom=fields[0];start=fields[1];end=fields[2];info=fields[9];
    domain=[i.split('"') for i in info.split(";")][0][1]
    if domain==prevdomain or prevdomain is None:
        gbed+="\t".join([chrom,start,end])+"\n"
    else:
        p=subprocess.Popen(['bedtools merge'],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        output,err = p.communicate(gbed)
        gbed="\t".join([chrom,start,end])+"\n"
        if prevdomain:
            for line in output.strip().split("\n"):
                print line + "\t" + prevdomain
        else:
            for line in output.strip().split("\n"):
                print line + "\t" + domain
    prevdomain=domain
p=subprocess.Popen(['bedtools merge'],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
output,err = p.communicate(gbed)
for line in output.strip().split("\n"):
    print line + "\t" + prevdomain

