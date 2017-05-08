import sys
import multiprocessing as mp
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

kcsq="ALLELE|Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM|cDNA_position".split("|") # later can extract header in unix so it isn't hard coded

p = mp.Pool(12)

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

infile=open(variants, "r")

def parseinfo(info):
    d={}
    for field in info.split(";"):
        vals=field.split("=")
        if len(vals)<2:
            d[vals[0]]=''
        else:
            d[vals[0]]=vals[1]
    return d

def cfilter(info, varstatus):
    if varstatus=="benign":
        return info['CLNSIG'] == '2' and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']
    if varstatus=="patho":
        return (info['CLNSIG'] == '5' or info['CLNSIG'] == '4') and info['CLNREVSTAT'] not in ['no_assertion', 'no_criteria', 'conf']

def pervariant(varianttuple):
    autopass=False; varpass=True
    filter,dom,haplo,rec,varstatus,clinvar,name,variant = varianttuple
    gene = ''; cct=0; outs=''
    fields=variant.strip().split("\t")
    original=fields[:8] #gnomad/clinvar
    filterby=fields[8:] #exac/gnomad, a set of variants to filter by
    oinfo = parseinfo(original[-1]); finfo = parseinfo(filterby[-1])
    ofilter = original[-2]; ffilter = filterby[-2]
    
    if "gnomad" in name:
        try:
            if oinfo['AS_FilterStatus'].split(",")[0] not in ["PASS", "SEGDUP", "LCR"]:
                return outs
        except KeyError:
            pass
    #if ofilter not in ["PASS", "SEGDUP", "LCR", "."]:
    #    return outs
    if clinvar:
        if not cfilter(oinfo, varstatus):
            return outs

    try:
        ocsqs = [dict(zip(kcsq, c.split("|"))) for c in oinfo['CSQ'].split(",")]
    except KeyError:
        return outs
    try:
        fcsqs = [dict(zip(kcsq, c.split("|"))) for c in finfo['CSQ'].split(",")]
    except KeyError:
        autopass=True

    for ocsq in (c for c in ocsqs if c['BIOTYPE'] == 'protein_coding'):
        if ocsq['Feature'] == '' or ocsq['EXON'] == '': continue
        if not u.isfunctional(ocsq): continue
        gene = ocsq['SYMBOL']
        if dom or haplo or rec:
            if gene not in dom_genes or haplo_genes or rec_genes:
                continue
        if "benign" in varstatus and "clinvar" in name or autopass:
            return "\t".join(original)
            break
        if filter:
            if ffilter is None or ffilter in ["PASS", "SEGDUP", "LCR"]:
                if exacver == "gnomad":
                    try:
                        if finfo['AS_FilterStatus'].split(",")[0] not in ["PASS", "SEGDUP", "LCR"]:
                            return "\t".join(original)
                    except KeyError:
                        pass
                for fcsq in (c for c in fcsqs if c['BIOTYPE'] == 'protein_coding'):
                    if not u.isfunctional(fcsq): continue
                    if fcsq['Feature'] == ocsq['Feature'] and (fcsq['Amino_acids'] == ocsq['Amino_acids'] or fcsq['Codons'] == ocsq['Codons']):
                        varpass=False

        if original[0] != 'X' and original[0] != 'Y' and varpass:
            return "\t".join(original)
            break

#for outs in p.imap_unordered(pervariant, ((filter,dom,haplo,rec,varstatus,clinvar,name,variant) for variant in infile)):
#    if outs:
#        print outs
for variant in infile:
    outs=pervariant((filter,dom,haplo,rec,varstatus,clinvar,name,variant))
    if outs:
        #print outs
        f.write(outs)
