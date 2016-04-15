import toolshed as ts
import itertools as it
import operator
from cyvcf2 import VCF

s=set()
exac = VCF('/uufs/chpc.utah.edu/common/home/u6000771/Projects/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz')
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
