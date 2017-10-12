from argparse import ArgumentParser
import numpy as np

parser=ArgumentParser()
#parser.add_argument("-p","--pfam", help="pfam families for histogram")
parser.add_argument("-i","--intersection", help="intersections between pfams and ccrs")
args=parser.parse_args()

prevfam=None; prevlen=None; prevccr=None; tccr={}; tlen=0.0; 

bins=range(0,101,10)
#pfams=open(args.pfam, 'r')
intersect=open(args.intersection, 'r')
for line in intersect: #sort intersection by family name
    fields=line.strip().split("\t")
    start=int(fields[1]); end=int(fields[2]); family=fields[3]; ccr=float(fields[-1]);
    lenh=end-start
    if family==prevfam:
        for i in range(0,len(bins)):
            if prevccr <= bins[i] and prevccr > bins[i-1]: # not putting >= bins[i-1] eliminates 0th percentile CCRs, also bug is bins[0-1] = bins[-1]
                try:
                    tccr[bins[i]]+=prevlen
                except KeyError:
                    tccr[bins[i]]=prevlen
                tlen+=prevlen
            #elif prevccr == 0:
            #    try:
            #        tccr[bins[i]]+=prevlen
            #    except KeyError:
            #        tccr[bins[i]]=prevlen
            #    tlen+=prevlen
    elif prevfam!=family and prevfam is not None:
        print prevfam + "\t" + str(tccr) + "\t" + str(tlen)
        tlen=0.0
        tccr={}
    prevfam=family
    prevlen=lenh
    prevccr=ccr 

print prevfam + "\t" + str(tccr) + "\t" + str(tlen)
