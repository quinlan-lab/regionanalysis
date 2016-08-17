from cyvcf2 import VCF
import itertools as it
import operator
import sys

exac = VCF(sys.argv[1])
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
def isfunctional(csq):
    return any(c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant',
                     'missense_variant', 'protein_altering_variant', 'frameshift_variant')
               for c in csq['Consequence'].split('&'))
print exac.raw_header,
for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    for v in viter:
        info = v.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        for csq in (c for c in csqs if c['CANONICAL'] == 'YES' and c['BIOTYPE'] == 'protein_coding' and isfunctional(c)):
            print str(v).strip()

