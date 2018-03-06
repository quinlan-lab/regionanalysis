from cyvcf2 import VCF, Writer
import toolshed as ts
import sys

var='varsets/shorthand.txt'
vcf='varsets/LQTSvariants.vcf' # had to add contigs to header for file to work (length of chromosomes)
lqts=VCF(vcf)
lqts.add_info_to_header({'ID': 'PATHOGENIC', 'Description': '0 if benign, 1 if pathogenic',
    'Type':'Character', 'Number': '1'})
f='varsets/LQTSvariantsannotated.vcf'
w = Writer(f, lqts)

for v in lqts:
    for toks in (x.rstrip('\r\n').split('\t') for x in ts.nopen(var)):
        if v.ID==toks[0]:
            if toks[1]=='pathogenic': #pathogenic
                v.INFO['PATHOGENIC']='1'
            else: #benign
                 v.INFO['PATHOGENIC']='0'
            w.write_record(v)

w.close()
