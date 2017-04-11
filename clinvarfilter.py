from cyvcf2 import VCF
import sys
import exacresiduals.utils as u
reload(sys)
sys.setdefaultencoding("ISO-8859-1")
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-c", "--clinvar", help="clinvar file")
parser.add_argument("-i", "--haploinsufficient", help="clingen dosage haplosufficiency 3 genes file, if necessary") #-h overlaps help
parser.add_argument("-d", "--dominant", help="dominant genes file, if necessary")
#parser.set_defaults(clinvar = '/scratch/ucgd/lustre/u1021864/serial/clinvar-vt-anno-vep.vcf.gz')
#parser.set_defaults(dominant = 'genescreens/ad_genecards_clean.txt')
#parser.set_defaults(haploinsufficient = 'genescreens/clingen_level3_genes_2015_02_27.tsv')
args=parser.parse_args()
dom=args.dominant
clinvar=args.clinvar
haplo=args.haploinsufficient

#this script makes appropriate pathogenic and benign clinvar variant files

clinvar=VCF(clinvar) #clinvar file, vcfanno'd and vep'd: /scratch/ucgd/lustre/u1021864/serial/clinvar-vt-anno-vep.vcf.gz
if dom:
    dom_genes = {} # since our model best represents dominant negative phenotypes, we are only interested in autosomal dominant genes here (from genecards)
    for line in open(dom): #genescreens/ad_genecards_clean.txt
        dom_genes[line.strip()] = 1
if haplo:
    haplo_genes = {} # since our model best represents dominant negative and haploinsufficient phenotypes, here we incorporate ClinGen 3 genes
    for line in open(haplo): #genescreens/clingen_level3_genes_2015_02_27.tsv
        haplo_genes[line.strip()] = 1
kcsq = clinvar["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
f1=open('benign.vcf','wb')
f2=open('patho.vcf','wb')
header=clinvar.raw_header
f1.write(header);f2.write(header)
for variant in clinvar:
    info = variant.INFO
    try:
        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
    except KeyError:
        continue
    for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
        if csq['Feature'] == '' or csq['EXON'] == '': continue
        if not u.isfunctional(csq): continue
        if info['CLNSIG'] == '2' and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']:
            f1.write(str(variant))
            break
        if info['CLNSIG'] == '5' and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']:
            exac_af = variant.INFO.get('ac_exac_all')
            gene_raw = variant.INFO.get('GENEINFO')
            if gene_raw is not None:
                gene = gene_raw.split(':')[0] 
            if exac_af is None and variant.CHROM != 'X' and variant.CHROM != 'Y':
                if dom or haplo:
                    if gene in dom_genes: #or gene in haplo_genes:
            #try:
            #    if info['max_aaf_all']>0.01: continue
            #except:
            #    pass
                        f2.write(str(variant))
                else:
                    f2.write(str(variant))
            break
