from __future__ import print_function
from cyvcf2 import VCF, Writer
import sys
key = "gnomad10x5_ccr"

vcf = VCF(sys.argv[1])
vcf.add_info_to_header({'ID': 'ccr+cadd', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'ccr+polyphen', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'ccr+gerp', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'ccr+grantham', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'ccr+misbadness', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})


wtr = Writer("-", vcf)

for v in vcf:

    ccr = v.INFO.get(key)
    if ccr is None:
        wtr.write_record(v)
        continue

    gerp = v.INFO.get("GERP")
    if gerp:
        if gerp < 2.5 or ccr > 90: # lowest is -12.36, anything below 0 is not conserved really
            v.INFO['ccr+gerp'] = ccr
        else:
            v.INFO['ccr+gerp'] = gerp * 16.18 #max is 6.18, 6.18 * 16.18 ~= 100

    polyphen = v.INFO.get("pp2hdiv")
    if polyphen:
        polyphen = polyphen * 100 #(scales from 0 to 1)
        if polyphen < 0.1 or polyphen > 0.9:
            v.INFO['ccr+polyphen'] = polyphen
        else:
            v.INFO['ccr+polyphen'] = ccr
        
    cadd = v.INFO.get("CADD") # max cadd is 99.0, 0 is min
    if cadd:
        if ccr == 0:
            v.INFO['ccr+cadd'] = cadd * 2
        if ccr > 90:
             v.INFO['ccr+cadd'] = ccr * 2 
        else:
             v.INFO['ccr+cadd'] = ccr + cadd

    grantham = v.INFO.get('Grantham') # 215 is max grantham, 0 is min
    if grantham:
        if grantham < 200 or ccr > 90:
            v.INFO['ccr+grantham'] = ccr
        else:
            v.INFO['ccr+grantham'] = 100*grantham/215.0 # scale it from 0 to 100
    
    mis_badness = v.INFO.get('mis_badness')
    if mis_badness:
        mis_badness = mis_badness * 100 #(scaled from 0 to 1)
        if mis_badness < 90 or ccr > 90:
            v.INFO['ccr+misbadness'] = ccr
        else:
            v.INFO['ccr+misbadness'] = mis_badness
            
    wtr.write_record(v)

wtr.close()
