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
parser.add_argument("-e", "--exacver", help = "specify 'exac' or 'gnomad'")
parser.add_argument("-f", "--filter", help = "filter on presence in exac/gnomad dataset", action = "store_true")
parser.add_argument("-n", "--name", help = "clinvar, mcrae, etc.")
parser.add_argument("-r", "--recessive", help = "recessive genes file, if necessary")
parser.add_argument("-c", "--clinvar", help = "working with clinvar data", action = "store_true")
parser.add_argument("-s", "--status", help = "variant status: benign or patho, type one of the two exactly")
#parser.set_defaults(variants = '/scratch/ucgd/lustre/u1021864/serial/variants-vep-anno-vt.vcf.gz')
#parser.set_defaults(dominant = 'genescreens/ad_genecards_clean.txt')
#parser.set_defaults(recessive = 'genescreens/ar_genecards_clean.txt')
#parser.set_defaults(haploinsufficient = 'genescreens/clingen_level3_genes_2015_02_27.tsv')
args = parser.parse_args()
dom = args.dominant
rec = args.recessive
haplo = args.haploinsufficient
variants = args.variant
exacver = args.exacver
filter = args.filter
name = args.name
clinvar = args.clinvar
varstatus = args.status

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
    for line in open(rec): #genescreens/clingen_level3_genes_2015_02_27.tsv
        rec_genes.add(line.strip())

folder = os.path.dirname(variants) + '/'
if folder == '/':
    folder = ''
f = open(folder + name + '-' + varstatus + '-' + exacver + '.txt' , 'wb')

infile = open(variants, "r")

def parseinfo(info):
    d={}
    for field in info.split(";"):
        vals = field.split("=")
        if len(vals)<2:
            d[vals[0]] = ''
        else:
            d[vals[0]] = vals[1]
    return d

def cfilter(info, varstatus):
    clnsig = info['CLNSIG']; clnrev=info['CLNREVSTAT']
    clnsig = re.split('\||,',clnsig); clnrev = re.split('\||,',clnrev)
    var = False
    if varstatus == "patho":
        for sig, rev in zip(clnsig,clnrev):
            if (sig == '5' or sig == '4'):
                if rev not in ['no_assertion', 'no_criteria']:
                    var = True
                elif rev in ['conf']:
                    return False
            else:
                return False
        return var
    if varstatus == "benign":
        for sig, rev in zip(clnsig,clnrev):
            if sig == '2':
                if rev not in ['no_assertion', 'no_criteria']:
                    var = True
                elif rev in ['conf']:
                    return False
            else:
                return False
        return var

def pervariant(varianttuple):
    autopass = False; varpass = True
    filter,dom,haplo,rec,varstatus,clinvar,name,original,filterby = varianttuple
    gene = ''
    oinfo = parseinfo(original[-1]); finfo = parseinfo(filterby[-1])
    ofilter = original[-2]; ffilter = filterby[-2]
    filterpos = filterby[1]
    
    if "gnomad" in name:
        try:
            if oinfo['AS_FilterStatus'].split(",")[0] not in ["PASS"]: # SEGDUP, LCR are debatable
                return False
        except KeyError:
            pass
    if ofilter not in ["PASS", "."]: # SEGDUP, LCR are debatable
        return False
    if clinvar:
        if not cfilter(oinfo, varstatus):
            return False
    try:
        ocsqs = [dict(zip(kcsq, c.split("|"))) for c in oinfo['CSQ'].split(",")]
    except KeyError:
        return False
    try:
        fcsqs = [dict(zip(kcsq, c.split("|"))) for c in finfo['CSQ'].split(",")]
    except KeyError:
        autopass = True
    if not any (c for c in ocsqs if c['BIOTYPE'] == 'protein_coding'): return False
    for ocsq in (c for c in ocsqs if c['BIOTYPE'] == 'protein_coding'):
        if ocsq['Feature'] == '' or ocsq['EXON'] == '':
            varpass = False
            continue
        if not u.isfunctional(ocsq):
            varpass = False
            continue
        gene = ocsq['SYMBOL']
        if dom or haplo or rec:
            if gene not in dom_genes or haplo_genes or rec_genes:
                return False
        if ("benign" in varstatus and "clinvar" in name) or autopass or filterby == "-1":
            return True
        if filter:
            if ffilter is None or ffilter in ["PASS", "."]: # SEGDUP, LCR are debatable
                if exacver == "gnomad":
                    try:
                        if finfo['AS_FilterStatus'].split(",")[0] not in ["PASS"]: # SEGDUP, LCR are debatable
                            return True
                    except KeyError:
                        pass
                for fcsq in (c for c in fcsqs if c['BIOTYPE'] == 'protein_coding'):
                    if not u.isfunctional(fcsq):
                        varpass = True
                        continue
                    if fcsq['Feature'] == ocsq['Feature'] and (fcsq['Amino_acids'] == ocsq['Amino_acids'] or fcsq['Codons'] == ocsq['Codons']):
                        return False
            varpass = True
        else:
            return True
    return varpass

#for outs in p.imap_unordered(pervariant, ((filter,dom,haplo,rec,varstatus,clinvar,name,variant) for variant in infile)):
#    if outs:
#        print outs

# change code to group variants by original POS, REF, ALT and filter all duplicates of a variant if even one is filtered; probably need to remove the generator unfortunately
# should be a large speed decrease
varprev = None; varpass = True; vars = []

for variant in infile:
    fields = variant.strip().split("\t")
    original = fields[:8] #gnomad/clinvar
    filterby = fields[8:-1] #exac/gnomad, a set of variants to filter by, now :-1 because there is a -wo field
    if varprev is None:
        pass
    elif varprev != original:
        for pvar in vars:
            varprev, filterprev = pvar
            varpass = pervariant((filter,dom,haplo,rec,varstatus,clinvar,name,varprev,filterprev))
            if not varpass:
                break
        if varpass and varprev[0] != 'X' and varprev[0] != 'Y':
            f.write("\t".join(varprev)+"\n")
        vars=[]
    varprev = original; filterprev = filterby
    vars.append((varprev, filterprev))

if varprev != None and varprev[0] != 'X' and varprev[0] != 'Y' and varpass:
    f.write("\t".join(varprev)+"\n")
