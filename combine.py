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
vcf.add_info_to_header({'ID': 'ccr+mcap', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'ccr+revel', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'ccr+mpc', 'Description': 'combined score',
                        'Type':'Float', 'Number': '1'})


wtr = Writer("-", vcf)

for v in vcf:
    ccr = v.INFO.get(key)
    gerp = v.INFO.get("GERP")
    polyphen = v.INFO.get("pp2hdiv")
    cadd = v.INFO.get("CADD") # max cadd is 99.0, 0 is min
    grantham = v.INFO.get('Grantham') # 215 is max grantham, 0 is min
    mis_badness = v.INFO.get('mis_badness')
    mcap = v.INFO.get('mcap')
    revel = v.INFO.get('revel')
    MPC = v.INFO.get('MPC')

    if gerp:
        gerp = gerp * 16.18 #max is 6.18, 6.18 * 16.18 ~= 100
        if gerp > 1.7 or ccr < 90 or ccr is None: # lowest is -12.36, anything below 1.7 is not conserved really
            v.INFO['ccr+gerp'] = gerp
        else:
            v.INFO['ccr+gerp'] = ccr
    elif ccr is not None:
        v.INFO['ccr+gerp'] = ccr

    if polyphen:
        polyphen = polyphen * 100 #(scales from 0 to 1)
        if polyphen < 0.1 or polyphen > 0.9 or ccr is None:
            v.INFO['ccr+polyphen'] = polyphen
        else:
            v.INFO['ccr+polyphen'] = ccr
    elif ccr is not None:
        v.INFO['ccr+polyphen'] = ccr
        
    if cadd:
        if ccr == 0 or ccr is None:
            v.INFO['ccr+cadd'] = cadd * 2
        elif ccr > 90:
             v.INFO['ccr+cadd'] = ccr * 2 
        else:
             v.INFO['ccr+cadd'] = ccr + cadd
    elif ccr is not None:
        v.INFO['ccr+cadd'] = ccr

    if grantham:
        grantham = 100*grantham/215.0 # scale from 0 to 100
        if grantham > 200 or ccr < 90 or ccr is None:
            v.INFO['ccr+grantham'] = grantham
        else:
            v.INFO['ccr+grantham'] = ccr
    elif ccr is not None:
        v.INFO['ccr+grantham'] = ccr
    
    if mis_badness:
        mis_badness = mis_badness * 100 #(scaled from 0 to 1)
        if mis_badness > 90 or ccr < 90 or ccr is None:
            v.INFO['ccr+misbadness'] = mis_badness
        else:
            v.INFO['ccr+misbadness'] = ccr
    elif ccr is not None:
        v.INFO['ccr+misbadness'] = ccr
            
    if mcap:
        mcap = mcap * 100 #(scaled from 0 to 1)
        if mcap > 90 or ccr < 90 or ccr is None:
            v.INFO['ccr+mcap'] = mcap
        else:
            v.INFO['ccr+mcap'] = ccr
    elif ccr is not None:
        v.INFO['ccr+mcap'] = ccr
     
    if revel:
        revel = revel * 100 #(scaled from 0 to 1)
        if revel < 20 or revel > 80 or ccr < 95 or ccr is None:
            v.INFO['ccr+revel'] = revel
        else:
            v.INFO['ccr+revel'] = ccr
    elif ccr is not None:
        v.INFO['ccr+revel'] = ccr
    
    if MPC:
        MPC = MPC * 25 #(scaled from 0 to ~4.3)
        if (MPC > 50 and MPC < 75) or MPC < 25 or ccr < 95 or ccr is None:
            v.INFO['ccr+mpc'] = MPC
        else:
            v.INFO['ccr+mpc'] = ccr
    elif ccr is not None:
        v.INFO['ccr+mpc'] = ccr
     
    wtr.write_record(v)

wtr.close()
