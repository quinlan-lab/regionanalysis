import subprocess as sp
import sys
import numpy as np
import atexit

@atexit.register
def runcalls():
    a=sys.argv[1]
    b=sys.argv[2]
    tf=sys.argv[3]

    arg1=set(a.split('.')[:-1])
    arg2=set(b.split('.')[:-1])
    for i in arg1:
        if i in arg2:
            exit()

    S = ['ad','ar','chromatin','clinvar','dnarep','dnase','doms','haplo','histone','meth','mirna','rna','tfbs']
    union=sorted(arg1 | arg2)
    mat=",".join(["1" if i in union else "0" for i in S])

    p = sp.Popen("sed '1d' " + a + " | cut -f -10", shell=True, stdout=sp.PIPE) 
    if tf == '1': 
        p = sp.Popen(['bedtools', 'intersect', '-a', 'stdin', '-b', b, '-sorted', '-wa'],stdin=p.stdout,stdout=sp.PIPE)
    if tf == '0':
        p = sp.Popen(["awk",'NR==FNR{a[$1]}; {for (i in a) if (i == $4) print}', b, a], stdout=sp.PIPE)

    f1=open(".".join(union)+".txt","w")
    p = sp.Popen('sort -k1,1 -k2,2n | uniq', shell=True, stdin=p.stdout, stdout=sp.PIPE)
    f1.write('#chrom\tstart\tend\tgene\texon\tlength\tresid\tresidperc\tcovperc\tzperc\tad\tar\tchromatin\tclinvar\tdnarep\tdnase\tdoms\thaplo\thistone\tmeth\tmirna\trna\ttfbs\n')
    p = sp.Popen("awk '{print $0,"+mat+"}' OFS='\t' - ", shell=True, stdin=p.stdout, stdout=f1)
