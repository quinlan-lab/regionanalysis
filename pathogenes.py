from cyvcf2 import VCF
import sys

#this generates file of genes with variants CLNSIG=5 and exac freq<0.01
patho=VCF(sys.argv[1]) #patho.vcf
kcsq=patho["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
f1=open('genescreens/pathogenes.txt','wb')
genes=set()
for variant in patho:
    info = variant.INFO
    try:
        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
    except KeyError:
        continue
    for csq in (c for c in csqs):
        genes.add(csq['SYMBOL'])

for i in genes:
    f1.write(i+"\n")
