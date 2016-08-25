import numpy as np
import toolshed as ts
from collections import OrderedDict, defaultdict

def rank(genes):
    kbm7,raji,jiyoye,k562=defaultdict(int),defaultdict(int),defaultdict(int),defaultdict(int)
    for gene in genes:
        if genes[gene][1] < 5e-2 and genes[gene][3] < 5e-2 and genes[gene][5] < 5e-2 and genes[gene][7] < 5e-2 and \
            genes[gene][0] < -1e-1 and genes[gene][2] < -1e-1 and genes[gene][4] < -1e-1 and genes[gene][6] < -1e-1:
            for agene in genes:
                if genes[gene][0]<genes[agene][0]:
                    kbm7[gene]+=1
                if genes[gene][2]<genes[agene][2]:    
                    raji[gene]+=1
                if genes[gene][4]<genes[agene][4]:
                    jiyoye[gene]+=1
                if genes[gene][6]<genes[agene][6]:
                    k562[gene]+=1
    genelist=defaultdict(int)
    for a,b,c,d in zip(kbm7,raji,jiyoye,k562):
        for e,f,g,h in zip(kbm7,raji,jiyoye,k562):
            if np.mean((kbm7[a],raji[a],jiyoye[a],k562[a]))<np.mean((kbm7[e],raji[e],jiyoye[e],k562[e])):
                genelist[a]+=1
    return genelist

it = ts.reader('ogfiles/wangetal2015.txt')
iterable = (i for i in it)
genes = {}
for gene in iterable:
    genes[gene['Gene']] = map(float,(gene['KBM7 CS'], gene['KBM7 adjusted p-value'], gene['Raji CS'], gene['Raji adjusted p-value'], gene['Jiyoye CS'], gene['Jiyoye adjusted p-value'], gene['K562 CS'],gene['K562 adjusted p-value']))

genelist=rank(genes)
for gene in genelist:
    print gene+"\t"+str(genelist[gene])

