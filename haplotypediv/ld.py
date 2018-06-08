import sys

CCR={}

for l in open(sys.argv[1], 'r'):
    A = l.rstrip().split()
    CCR[A[0]] = float(A[1])

LD={}
for l in open(sys.argv[2], 'r'):
    A = l.rstrip().split()
    LD[A[0]] = float(A[5])


for gene in LD:
    if gene in CCR:
        #print gene,LD[gene],CCR[gene]
        print CCR[gene]

