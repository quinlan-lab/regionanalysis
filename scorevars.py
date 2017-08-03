import sys
import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("-x", "--variants", help = "variant file")
parser.add_argument("-c", "--ccr", help = "ccr file")
parser.add_argument("-a", "--average", help = "created weighted average based on total length of intersection of ccrs with variant", action = "store_true")
parser.add_argument("-n", "--name", help = "clinvar, mcrae, etc.")
#parser.set_defaults(variants = '/scratch/ucgd/lustre/u1021864/serial/clinvar-patho-exac.vcf.gz')
#parser.set_defaults(ccr = 'exacresiduals/results/exacv1newweight/weightedresiduals-cpg-novariant.txt')
args = parser.parse_args()
variants = args.variants
ccr = args.ccr
average = args.average
name = args.name

#bedtools intersect the two w/ -wo
def intersect(variants, regions, wo = True):
    def killproc(p):
        try:
            p.kill()
        except OSError:
            pass
    l = ['bedtools', 'intersect', '-a', variants, '-b', regions, '-sorted']
    if wo:
        l.append('-wo')
    p1 = sp.Popen(l, stdout = sp.PIPE)
    output,error = p1.communicate()
    killproc(p1)
    return output.strip()

# -wo bp overlap divided by length of variant (i.e., length of REF) and multiplied by CCR score to get new variant score
varprev = None; varprevscore = 0.0; scores=[]; lengths=[]
for variant in intersect(variants, ccr).split("\n"):
    fields = variant.strip().split("\t")
    var = [fields[x] for x in [0,1,3,4]]
    overlap = int(fields[-1]) + 1 # until aaron fixes bedtools off-by-one bug, need to add +1
    ccrscore = float(fields[-2])
    if varprev is None:
        pass # will be dealt with at next pass through loop
    elif varprev == var:  
        lengths.append(overprev)
        scores.append(ccrprev)
    elif varprev != var:
        lengths.append(overprev)
        scores.append(ccrprev)
        if average:
            varprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
        print "%.3f" % varprevscore
    varprev = var; overprev = overlap; ccrprev = ccrscore
lengths.append(overprev)
scores.append(ccrprev)
if average:
    varprevscore=sum([a*b for a,b in zip(scores,lengths)])/sum(lengths)
print "%.3f" % varprevscore
