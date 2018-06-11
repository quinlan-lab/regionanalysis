from cyvcf2 import VCF
import sys
reload(sys)
sys.setdefaultencoding('utf8')

def isfunctional(csqs):
    for csq in csqs.split(","):
        eff = csq.split("|", 2)[0]
        for c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon', 'rare_amino_acid',
                     'missense', 'protein_altering', 'frameshift', 'inframe_insertion', 'inframe_deletion'):
            if c in eff or (('splice_donor' in eff or 'splice_acceptor' in eff) and 'coding_sequence' in eff):
                return True
    return False

vcf = VCF(sys.argv[1])
print vcf.raw_header,
for v in vcf:
    if v.REF == v.ALT[0]:
        continue
    csq = v.INFO.get("BCSQ")
    if csq is None or not isfunctional(csq) or v.INFO.get('_exclude'):
        continue
    print str(v).strip()
