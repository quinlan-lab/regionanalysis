from __future__ import print_function
from cyvcf2 import VCF
import sys
reload(sys)
sys.setdefaultencoding('utf8')

variants=VCF(sys.argv[1])
f1=open("expert.clinvar.patho.vcf", "w")
f2=open("other.clinvar.patho.vcf", "w")

header = """##fileformat=VCFv4.1
##source=pathoscore
##reference=GRCh37
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""
print(header,file=f1)
print(header,file=f2)

def cfilter(info, varstatus): # refer to https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
    exclude = "criteria_provided,_conflicting_interpretations no_assertion_criteria_provided no_assertion_provided no_interpretation_for_the_single_variant".split()
    pathogenic = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']
    confident = ['reviewed_by_expert_panel', 'practice_guideline']
    semiconfident = ['criteria_provided,_multiple_submitters,_no_conflicts', 'criteria_provided,_single_submitter']
    try:
        clnsig = info['CLNSIG']
        clnrev = info['CLNREVSTAT']
    except KeyError:
        return False
    if "pathogenic" in varstatus and clnsig in pathogenic:
            if "confident" in varstatus and clnrev in confident:
                return True
            elif "confident" not in varstatus and clnrev in semiconfident:
                return True
            else:
                return False

for v in variants:
    info=v.INFO
    if cfilter(info, "confident pathogenic"):
        f1.write(str(v))
    if cfilter(info, "other pathogenic"):
        f2.write(str(v))
