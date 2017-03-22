import toolshed as ts
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vep", help="prepare vep format", action="store_true", default=False)
args=parser.parse_args()
vep=args.vep

var='VVPvariants/varsets/LQTS_variants.txt'
trans='VVPvariants/varsets/transcriptTable.txt'

if not vep:
    f=open('VVPvariants/varsets/shorthand.txt','w')
else:
    f=open('VVPvariants/varsets/vepweb.txt','w')
def readgenes(trans):
    genes=defaultdict(str)
    for fields in (x.rstrip('\r\n').split("\t") for x in ts.nopen(trans)):
        gene=fields[0]; transcript=fields[1]
        genes[gene]=transcript
    return genes
        
genes=readgenes(trans)
for i, d in enumerate(ts.reader(var)):
    variant=d['Coding']; gene=d['Gene']; transcript=genes[gene]
    line=transcript+":"+variant # so VEP web interface can translate to VCF
    if vep:
        f.write(line+"\n")
    else:
        if 'LQTS' in d['ICC.diseases']:
            f.write(line+"\t"+"pathogenic"+"\n")
        else:
            if d['Classification'] == "Benign":
                f.write(line+"\t"+"benign"+"\n")
