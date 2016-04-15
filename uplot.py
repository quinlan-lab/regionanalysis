import pyupset
import sys
import numpy as np

S = ['ad','ar','chromatin','clinvar','dnarep','dnase','doms','haplo','histone','meth','mirna','rna','tfbs']
f=sys.argv[1]
overlaps=open(f,"r")
outfile=open("upset.txt","w")
outfile.write('chrom\tstart\tend\tgene\tad\tar\tchromatin\tclinvar\tdnarep\tdnase\tdoms\thaplo\thistone\tmeth\tmirna\trna\ttfbs\n')

for line in overlaps:

    fields=line.strip().split("\t") #assumes filenames are already sorted, and they are thanks to me submitting them with * glob
    names=fields[-1] 
    mat="\t".join(["1" if i in names else "0" for i in S])
    outfile.write("\t".join(fields[:-1])+"\t"+mat+"\n")

overlaps.close()
outfile.close()

from collections import defaultdict

infile=open("upset.txt","r")
header=infile.readline().strip().strip("#").split("\t")
uset=defaultdict(list)
for line in infile:
    for i,j in zip(header,line):
        uset[i].append(j)
