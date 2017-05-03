from cyvcf2 import VCF
import sys
import exacresiduals.utils as u
reload(sys)
sys.setdefaultencoding("ISO-8859-1")
import argparse
import os

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

folder=os.path.dirname(variants)+'/'
if folder == '/':
    folder = ''
f=open(folder+name+'-'+varstatus+'-'+exacver+'.vcf','wb')
variants=VCF(variants) #variants file, vcfanno'd and vep'd: /scratch/ucgd/lustre/u1021864/serial/clinvar_20170104-vep-anno-vt.vcf.gz
header=variants.raw_header
f.write(header)

gene = ''
kcsq = variants["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

prevpos=-1; idx=0
for variant in variants:
    info = variant.INFO
    if prevpos == variant.POS:
        idx+=1
    else:
        idx=0
        prevpos = variant.POS
    if clinvar:
        if not cfilter(info, varstatus):
            continue
    gnomad_af = variant.INFO.get('ac_gnomad_all')
    gnomad_filter = variant.INFO.get('gnomad_filter')
    exac_af = variant.INFO.get('ac_exac_all')
    exac_filter = variant.INFO.get('exac_filter')
    alt = variant.ALT[0]
    try:
        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
    except KeyError:
        continue
    cct=0
    if "gnomad" in name:
        try:
            if info['AS_FilterStatus'].split(",")[0] not in ["PASS", "SEGDUP", "LCR"]:
                cct=-1
                continue
        except KeyError:
            pass
    for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
        if csq['Feature'] == '' or csq['EXON'] == '': continue
        if not u.isfunctional(csq): continue
        cct+=1
        gene = csq['SYMBOL']
        if dom or haplo or rec:
            if gene not in dom_genes or haplo_genes or rec_genes:
                cct-=1
                continue
        if "benign" in varstatus and "clinvar" in name:
            f.write(str(variant))
            break
        try:
            exac_csqs = [dict(zip(kcsq, c.split("|"))) for c in info['exac_csq'].split(",")]
            gnomad_csqs = [dict(zip(kcsq, c.split("|"))) for c in info['gnomad_csq'].split(",")]
        except KeyError:
            pass
        if gene is not None:
            if filter:
                if "exac" in exacver:
                    if exac_af is not None and (exac_filter is None or exac_filter in ["PASS", "SEGDUP", "LCR"]):
                        for csq2 in exac_csqs:
                            if csq2['Feature'] == csq['Feature'] and (csq2['Amino_acids'] == csq['Amino_acids'] or csq2['Codons'] == csq['Codons']):
                                cct=0
                                break
                    #try:
                        #if info['max_aaf_all']>0.01: continue
                    #except:
                    #    pass
                elif "gnomad" in exacver:
                    try:
                    # add fix here to pull in alts from gnomad, not oldmultiallelic
                        alts=[j for i in info['gnomad_oldmultiallelic'].split(",") for j in i.split("/")[1:]]
                        if alt in alts:
                            altind=alts.index(alt)
                            if info['gnomad_filterstatus'].split(",")[altind] in ["PASS", "SEGDUP", "LCR"]:
                                cct=-1
                                continue
                    except KeyError:
                        pass
                    if gnomad_af is not None and (gnomad_filter is None or gnomad_filter in ["PASS", "SEGDUP", "LCR"]):
                        for csq2 in gnomad_csqs:
                            if csq2['Feature'] == csq['Feature'] and (csq2['Amino_acids'] == csq['Amino_acids'] or csq2['Codons'] == csq['Codons']):
                                cct=0
                                break
    if variant.POS == 949395:
        print variant
    if variant.CHROM != 'X' and variant.CHROM != 'Y' and cct!=0:
        f.write(str(variant))
