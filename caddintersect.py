from cyvcf2 import VCF
import sys
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-p", "--pathogenic", help="pathogenic variant file")
parser.set_defaults(pathogenic = 'patho.vcf.gz')
parser.add_argument("-b", "--benign", help="benign variant file")
parser.set_defaults(benign = 'benign.vcf.gz')
parser.add_argument("-c", "--cadd", help="cadd file, then cadd indel file")
parser.set_defaults(cadd = '/scratch/ucgd/lustre/u1021864/serial/CADD.vcf.gz')
parser.add_argument("-d", "--caddindels", help="cadd indel file")
parser.set_defaults(caddindels = '/scratch/ucgd/lustre/u1021864/serial/CADDindels.vcf.gz')
parser.add_argument("-f", "--files", nargs=2, help="output file locations")
args=parser.parse_args()
patho=args.pathogenic
benign=args.benign
cadd=args.cadd
caddindels=args.caddindels
files=args.files

cadd = VCF(cadd)
caddindels = VCF(caddindels)
patho = VCF(patho)
benign = VCF(benign)

f1=open(files[0], 'w')
f2=open(files[1], 'w')

for var in patho:
    position=var.CHROM+":"+str(var.POS)+"-"+str(var.POS)
    for v in cadd(position):
        if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT, and var.ALT is already decomposed
            f1.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score
    for v in caddindels(position):
        if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT, and var.ALT is already decomposed
            f1.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score

for var in benign:
    position=var.CHROM+":"+str(var.POS)+"-"+str(var.POS)
    for v in cadd(position):
        if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT
            f2.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score
    for v in caddindels(position):
        if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT, and var.ALT is already decomposed
            f2.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score
