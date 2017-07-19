import exacresiduals.utils as u
from cyvcf2 import VCF
import sys

VCF_PATH=sys.argv[1]

gnomad = VCF(VCF_PATH)
kcsq = gnomad["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
#kcsq = "ALLELE|Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM|cDNA_position".split("|") # later can extract header in unix so it isn't hard coded

def perchrom(vcf_chrom):
    outs=[0,0]
    vcf, chrom = vcf_chrom
    viter = VCF(VCF_PATH)(chrom)
    for v in viter:
        outs[1]+=1
        info = v.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'): # getting duplicate rows because of this, wastes memory and potentially compute time, could remove and replace with just if isfunctional, add to rows then move on?
            # skipping intronic
            if csq['Feature'] == '' or csq['EXON'] == '' : continue #or csq['cDNA_position'] == '': continue
            if u.isfunctional(csq): 
                outs[0]+=1
                break
    return outs

missense=0; total=0; chroms = range(1,23); chroms.append("X"); chroms.append("Y")
import multiprocessing as mp
p = mp.Pool(12)
for outs in p.imap_unordered(perchrom, ((VCF, str(chrom)) for chrom in chroms)):
        missense+=outs[0]; total+=outs[1]

print missense, total
