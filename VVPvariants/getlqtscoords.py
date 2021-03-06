import toolshed as ts
from collections import defaultdict
import pyhgvs as hgvs
from pygr.seqdb import SequenceFileDB
from pyhgvs.utils import read_transcripts


var='varsets/LQTS_variants.txt'
trans='varsets/transcriptTable.txt'
f=open('varsets/lqtsvariants.vcf','w')

infile='../Homo_sapiens.GRCh37.75.genePred' #obtained by running gtfToGenePred UCSC binary on the GTF of GRCh37 from ENSEMBL. added 0 as placeholder undocumented id field to the front of each line using sed
with open(infile) as reffile:
    transcripts=read_transcripts(reffile)
def get_transcript(name):
    return transcripts.get(name)
def translate(variant,transcripts,get_transcript):
    genome = SequenceFileDB('../hg19.fa') #pip install bsddb3 is required; fa in analysis folder
    try:
        chrom, offset, ref, alt = hgvs.parse_hgvs_name(variant, genome, get_transcript=get_transcript)
    except:
        return 1
    return chrom, offset, ref, alt

def readgenes(trans):
    genes=defaultdict(str)
    for fields in (x.rstrip('\r\n').split("\t") for x in ts.nopen(trans)):
        gene=fields[0]; transcript=fields[1]
        genes[gene]=transcript
    return genes
        
genes=readgenes(trans)
 # had to add contigs to header for file to work (length of chromosomes)
f.write('##fileformat=VCFv4.0\n')
f.write('##reference=GRCh37\n')
f.write('''##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>\n''')
f.write('##INFO=<ID=PATHOGENIC,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
for i, d in enumerate(ts.reader(var)):
    variant=d['Coding']; gene=d['Gene']; transcript=genes[gene]
    ID=transcript+":"+variant # so VEP web interface can translate to VCF
    try:
        fields=list(translate(ID,transcripts,get_transcript))
        fields.insert(2,ID)
        line='\t'.join(map(str,fields))
    except TypeError:
        continue
    if 'LQTS' in d['ICC.diseases']:
        f.write(line+"\t"+"PATHOGENIC=1"+"\n")
    else:
        if d['Classification'] == "Benign":
            f.write(line+"\t"+"PATHOGENIC=0"+"\n")
