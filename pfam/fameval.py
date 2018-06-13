from argparse import ArgumentParser
import numpy as np
from itertools import groupby
from operator import itemgetter
from collections import defaultdict

parser=ArgumentParser()
#parser.add_argument("-p","--pfam", help="pfam families for histogram")
parser.add_argument("-i","--intersection", help="intersections between pfams and ccrs")
args=parser.parse_args()

tccr=defaultdict(lambda: [0,0]); tlen=0.0; tcnt=0

bins=range(0,101,5)
#pfams=open(args.pfam, 'r')
intersect=open(args.intersection, 'r')
intersections=[]
for line in intersect:
    fields=line.strip().split("\t")
    start=int(fields[1]); end=int(fields[2]); family=fields[3]; gene=fields[8]; ranges=fields[11]; ccr=float(fields[-1])
    lenh=end-start
    intersections.append((ccr,lenh,gene,ranges,family))
    sorter = itemgetter(-1)
    grouper = itemgetter(-1) 
for key, grp in groupby(sorted(intersections, key = sorter), grouper): #sort intersection by family name
    grp=list(grp)
    family=grp[0][-1]
    for key2, grp2 in groupby(sorted(grp, key = itemgetter(2,3)), itemgetter(2,3)):
        grp2 = list(grp2)
        for j in range(0,len(bins)):
            ccr = grp2[0][0]
            if ccr <= bins[j] and ccr > bins[j-1]: # not putting >= bins[j-1] eliminates 0th percentile CCRs, because bins[0-1] = bins[-1]
                tccr[bins[j]][1]+=1
                tcnt+=1
    for i, elem in enumerate(grp):  
        for j in range(0,len(bins)):
            ccr = grp[i][0]
            length = grp[i][1]
            if ccr <= bins[j] and ccr > bins[j-1]: # not putting >= bins[j-1] eliminates 0th percentile CCRs, because bins[0-1] = bins[-1]
                tccr[bins[j]][0]+=length
                tlen+=length
    print family + "\t" + str(dict(tccr)) + "\t" + str(tlen) + "\t" + str(tcnt)
    tccr=defaultdict(lambda: [0,0]); tlen=0.0; tcnt=0
