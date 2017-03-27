from cyvcf2 import VCF
import sys
import exacresiduals.utils as u
reload(sys)
sys.setdefaultencoding("ISO-8859-1")
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--dominant", help="dominant genes only", action="store_true", default=False)
parser.add_argument("-f", "--files", help="clinvar file, then dominant genes file, if necessary", nargs="*")
parser.set_defaults(files = ['/scratch/ucgd/lustre/u1021864/serial/clinvar-vt-anno-vep.vcf.gz','genescreens/ad_genecards_clean.txt'])
args=parser.parse_args()
dom=args.dominant
files=args.files

#this script makes appropriate pathogenic and benign clinvar variant files

clinvar=VCF(args.files[0]) #clinvar file, vcfanno'd and vep'd: /scratch/ucgd/lustre/u1021864/serial/clinvar-vt-anno-vep.vcf.gz
if dom:
    dom_genes = {} # since our model best represents dominant negative phenotypes, we are only interested in autosomal dominant genes here (from genecards)
    for line in open(args.files[1]): #genescreens/ad_genecards_clean.txt
        dom_genes[line.strip()] = 1
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
                if dom:
                    if gene in dom_genes:
            #try:
            #    if info['max_aaf_all']>0.01: continue
            #except:
            #    pass
                        f2.write(str(variant))
                else:
                    f2.write(str(variant))
            break
