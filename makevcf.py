import sys
from pyfaidx import Fasta
import random

print "##fileformat=VCFv4.1"
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
fa=Fasta('/scratch/ucgd/lustre/u1021864/serial/hg19.fasta')
f=open('exacresiduals/fullresiduals.txt','r')
f.readline()
nucs=set(['A','C','G','T'])
for line in f:
        fields=line.strip().split("\t")
        chrom=fields[0];start=int(fields[1]);end=int(fields[2]);
        ref=fa['chr'+chrom][start].seq
        alt=random.sample(nucs-set([ref]),1)[0]
        print "\t".join(map(str,[chrom,start+1,'.',ref,alt,100,'PASS','AC=1']))
        ref=fa['chr'+chrom][end].seq
        alt=random.sample(nucs-set([ref]),1)[0]
        print "\t".join(map(str,[chrom,end+1,'.',ref,alt,100,'PASS','AC=1']))
