import toolshed as ts
import itertools as it
import operator
from cyvcf2 import VCF

s=set()
exac = VCF('/scratch/ucgd/lustre/u1021864/serial/ExAC.r1.sites.vt.vep.vcf.gz')
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    for v in viter:
        if not (v.FILTER is None or v.FILTER == "PASS"):
            continue
        info = v.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        for csq in (c for c in csqs if c['CANONICAL'] == 'YES' and c['BIOTYPE'] == 'protein_coding'):
            if csq['Feature'] not in s:
                s.add(csq['Feature'])
                print csq['Feature']
