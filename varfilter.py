from cyvcf2 import VCF
import sys
import exacresiduals.utils as u
reload(sys)
sys.setdefaultencoding("ISO-8859-1")
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-x", "--variant", help="variant file")
parser.add_argument("-d", "--dominant", help="dominant genes file, if necessary")
parser.add_argument("-i", "--haploinsufficient", help="clingen dosage haplosufficiency 3 genes file, if necessary") #-h overlaps help
parser.add_argument("-e", "--exacver", help="specify 'exac' or 'gnomad'")
parser.add_argument("-f", "--filter", help="filter on presence in exac/gnomad dataset", action="store_true")
parser.add_argument("-n", "--name", help="clinvar, mcrae, etc.")
parser.add_argument("-r", "--recessive", help="recessive genes file, if necessary")
parser.add_argument("-c", "--clinvar", help="working with clinvar data", action="store_true")
parser.add_argument("-s", "--status", help="variant status: benign or patho, type one of the two exactly")
parser.add_argument("-g", "--gnomad", help="generate benigns from non-exac gnomad variants", action="store_true")
#parser.set_defaults(variants = '/scratch/ucgd/lustre/u1021864/serial/variants-vep-anno-vt.vcf.gz')
#parser.set_defaults(dominant = 'genescreens/ad_genecards_clean.txt')
#parser.set_defaults(recessive = 'genescreens/ar_genecards_clean.txt')
#parser.set_defaults(haploinsufficient = 'genescreens/clingen_level3_genes_2015_02_27.tsv')
args=parser.parse_args()
dom=args.dominant
rec=args.recessive
haplo=args.haploinsufficient
variants=args.variant
exacver=args.exacver
filter=args.filter
name=args.name
clinvar=args.clinvar
varstatus=args.status
gnomad=args.gnomad


def cfilter(info, varstatus):
    if varstatus=="benign":
        return info['CLNSIG'] == '2' and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']
    if varstatus=="patho":
        return (info['CLNSIG'] == '5' or info['CLNSIG'] == '4') and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']

#this script makes appropriate pathogenic and benign variants variant files

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
    for line in open(rec): #genescreens/clingen_level3_genes_2015_02_27.tsv
        rec_genes.add(line.strip())

variants=VCF(variants) #variants file, vcfanno'd and vep'd: /scratch/ucgd/lustre/u1021864/serial/variants-vt-anno-vep.vcf.gz

f=open(name+'-'+varstatus+'-'+exacver+'.vcf','wb')
header=variants.raw_header
kcsq = variants["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
f.write(header)

gene = ''

for variant in variants:
    info = variant.INFO
    gene_raw = variant.INFO.get('GENEINFO')
    if gene_raw is not None:
        gene = gene_raw.split(':')[0]
    if dom or haplo or rec:
        if gene not in dom_genes or haplo_genes or rec_genes:
            continue
    if clinvar:
        if not cfilter(info, varstatus):
            continue
#    if info['CLNSIG'] == '2' and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']:
#        f1.write(str(variant))
#        continue
    exac_af = variant.INFO.get('ac_exac_all')
    gene_raw = variant.INFO.get('GENEINFO')
    try:
        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
    except KeyError:
        continue
    for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
        if csq['Feature'] == '' or csq['EXON'] == '': continue
        if not u.isfunctional(csq): continue
        if varstatus == "benign":
            f.write(str(variant))
            break
        if gene_raw is not None:
            gene = gene_raw.split(':')[0] 
            if filter:
                if exac_af is None and variant.CHROM != 'X' and variant.CHROM != 'Y':
                #try:
                    #if info['max_aaf_all']>0.01: continue
                #except:
                #    pass
                    f.write(str(variant))
            else:
                f.write(str(variant))
        break
