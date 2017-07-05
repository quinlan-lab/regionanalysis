from cyvcf2 import VCF
import sys
import math
import subprocess as sp
import numpy as np
import argparse
import tabix

parser=argparse.ArgumentParser()
parser.add_argument("-p", "--pathogenic", help="pathogenic variant file, pathogenic output file", nargs=2)
#parser.set_defaults(pathogenic = 'patho.vcf.gz')
parser.add_argument("-b", "--benign", help="benign variant file, benign output file", nargs=2)
#parser.set_defaults(benign = 'benign.vcf.gz')
parser.add_argument("-r", "--ccrs", help="ccrs file", nargs=2)
#parser.set_defaults(ccrs = 'exac-ccrs.bed.gz')
parser.add_argument("-c", "--cadd", help="cadd file or mpc file")
#parser.set_defaults(cadd = '/scratch/ucgd/lustre/u1021864/serial/CADD.vcf.gz')
parser.add_argument("-d", "--caddindels", help="cadd indel file")
#parser.set_defaults(caddindels = '/scratch/ucgd/lustre/u1021864/serial/CADDindels.vcf.gz')
parser.add_argument("-i", "--intersection", action="store_true", help="output intersected variants")
args=parser.parse_args()
pathogenic=args.pathogenic
benigns=args.benign
cadd=args.cadd
ccrs=args.ccrs
caddindels=args.caddindels
intersect=args.intersection

cadd = VCF(cadd)
if caddindels:
    caddindels = VCF(caddindels)

def intersect(genes, variants, wb = True):
    def killproc(p):
        try:
            p.kill()
        except OSError:
            pass
    l = ['bedtools', 'intersect', '-a', variants, '-b', regions, '-sorted']
    if wb:
        l.append('-wb')
    p1 = sp.Popen(l, stdout = sp.PIPE)
    output,error = p1.communicate()
    killproc(p1)
    return output.strip()


if pathogenic:
    patho = VCF(pathogenic[0])
    f1=open(pathogenic[1], 'w')
    if intersect:
        f1.write(patho.raw_header)
    for var in patho:
        position=var.CHROM+":"+str(var.POS)+"-"+str(var.POS)
        for v in cadd(position):
            if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT, and var.ALT is already decomposed
                if intersect:
                    f1.write(str(v)) # 7th column is always CADD PHRED score or MPC
                else:
                    f1.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score or MPC
        if caddindels:
            for v in caddindels(position):
                if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT, and var.ALT is already decomposed
                    if intersect:
                        f1.write(str(v)) # 7th column is always CADD PHRED score or MPC
                    else:
                        f1.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score or MPC
if benigns:
    benign = VCF(benigns[0])
    f2=open(benigns[1], 'w')
    if intersect:
        f2.write(benign.raw_header)
    for var in benign:
        position=var.CHROM+":"+str(var.POS)+"-"+str(var.POS)
        for v in cadd(position):
            if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT
                if intersect:
                    f2.write(str(v)) # 7th column is always CADD PHRED score or MPC
                else:
                    f2.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score or MPC
        if caddindels:
            for v in caddindels(position):
                if v.REF == var.REF and v.ALT[0] == var.ALT[0]: # using v.ID because there is no ID field so ID field is REF and REF is ALT, and var.ALT is already decomposed
                    if intersect:
                        f2.write(str(v)) # 7th column is always CADD PHRED score or MPC
                    else:
                        f2.write(str(v).strip().split('\t')[6]+'\n') # 7th column is always CADD PHRED score or MPC

if ccrs:
    regions=[]
    ccrcadd = tabix.open(ccrs[0])
    f3=open(ccrs[1], 'w')
    f3.write(cadd.raw_header)
    for v in cadd:
        position=v.CHROM+":"+str(v.POS)+"-"+str(v.POS)
        phredcadd=float(str(v).strip().split('\t')[6])
        f3.write("\t".join(str(v).strip().split('\t')[:6])+"\t")
        for region in ccrcadd.querys(position):
            if float(region[-1]) == 100:
                phredccr=100
            else:
                phredccr=-10*math.log10(1-float(region[-1])/100)
            regions.append(phredccr)
        phredccr=np.mean(regions)
        if phredccr == 0:
            f3.write(str(0)+"\n")
        elif phredccr > 10 or phredcadd > 30:
            if math.isnan(phredccr): # if no ccr score, just use cadd score
                f3.write(str(phredcadd)+"\n")
            else:
                f3.write(str(2*np.max([phredccr,phredcadd]))+"\n")
        else:
            if math.isnan(phredccr): # if no ccr score, just use cadd score
                f3.write(str(phredcadd)+"\n")
            else:
                f3.write(str(phredcadd+phredccr)+"\n")
        regions=[]
