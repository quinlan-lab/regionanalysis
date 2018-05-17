import sys

CCR={}

for l in open(sys.argv[1], 'r'):
    A = l.rstrip().split()
    CCR[A[0]] = float(A[1])

ES={}
for l in open(sys.argv[2], 'r'):
    A = l.rstrip().split()
    ES[A[0]] = 1


for gene in ES:
    if gene in CCR:
        print CCR[gene]

