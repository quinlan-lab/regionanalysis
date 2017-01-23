from cyvcf2 import VCF
import sys
import exacresiduals.utils as u
reload(sys)
sys.setdefaultencoding("ISO-8859-1")
#this script makes appropriate pathogenic and benign clinvar variant files

def func(string):
    return any([i in string for i in 'stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant', 'missense_variant', 'protein_altering_variant', 'frameshift_variant'])

clinvar=VCF(sys.argv[1]) #clinvar file, vcfanno'd and vep'd: /scratch/ucgd/lustre/u1021864/serial/clinvar-anno.vep.vcf.gz
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
            try:
                if info['max_aaf_all']>0.01: continue
            except:
                pass
            f2.write(str(variant))
            break
