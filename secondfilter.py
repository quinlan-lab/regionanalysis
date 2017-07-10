from cyvcf2 import VCF
import sys
import exacresiduals.utils as u
reload(sys)
sys.setdefaultencoding("ISO-8859-1")
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("-x", "--variant", help = "variant file")
parser.add_argument("-d", "--dominant", help = "dominant genes file, if necessary")
parser.add_argument("-i", "--haploinsufficient", help = "clingen dosage haplosufficiency 3 genes file, if necessary") #-h overlaps help
parser.add_argument("-f", "--filter", help = "filter on presence in exac/gnomad dataset")
parser.add_argument("-r", "--recessive", help = "recessive genes file, if necessary")
parser.add_argument("-c", "--clinvar", help = "working with clinvar data", action = "store_true")
#parser.set_defaults(variants = '/scratch/ucgd/lustre/u1021864/serial/variants-vep-anno-vt.vcf.gz')
#parser.set_defaults(dominant = 'genescreens/ad_genecards_clean.txt')
#parser.set_defaults(recessive = 'genescreens/ar_genecards_clean.txt')
#parser.set_defaults(haploinsufficient = 'genescreens/clingen_level3_genes_2015_02_27.tsv')
args = parser.parse_args()
dom = args.dominant
rec = args.recessive
haplo = args.haploinsufficient
variants = args.variant
filter = float(args.filter)
clinvar = args.clinvar

kcsq = "ALLELE|Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM|cDNA_position".split("|") # later can extract header in unix so it isn't hard coded

dom_genes = set() # since our model best represents dominant negative phenotypes, we are only interested in autosomal dominant genes here (from genecards)
haplo_genes = set() # since our model best represents dominant negative and haploinsufficient phenotypes, here we incorporate ClinGen 3 genes
rec_genes = set() # since our model may capture recessive phenotypes as well, here we use AD genes
if dom:
    for line in open(dom): #genescreens/ad_genecards_clean.txt
        dom_genes.add(line.strip())
if haplo:
    for line in open(haplo): #genescreens/clingen_level3_genes_2015_02_27.tsv
        haplo_genes.add(line.strip())
if rec:
    for line in open(rec): #genescreens/all_ar.tsv
        rec_genes.add(line.strip())

variants=VCF(variants)
f = open('/scratch/ucgd/lustre/u1021864/serial/'+'gnomad-benign-exac-filtered.vcf', 'wb')
header = variants.raw_header        
f.write(header)        

for variant in variants:
    info = variant.INFO
    try:        
        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]        
    except KeyError:      
        continue
    if filter > 1: # AF max is 1, so we can go by allele count
        if info['AC'] < filter:
            continue
    elif info['AF'] < filter:
        continue
    for csq in csqs:
        gene = csq['SYMBOL']
        if dom or haplo or rec:
            if gene in dom_genes or haplo_genes or rec_genes:
                f.write(str(variant))
                break
        else:
            f.write(str(variant))
            break
