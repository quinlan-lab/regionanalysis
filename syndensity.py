# script for calculating density of C-T synonymous variants in CpGs only versus CpG density
from cyvcf2 import VCF
import utils as u
from gzip import open
from pyfaidx import Fasta
import toolshed as ts

gnomad=VCF('exacresiduals/data/gnomad-vep-vt.vcf.gz') 
kcsq = gnomad["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
rfile = 'exacresiduals/results/gnomAD10x.5syn/exac-regions-novariant.txt'
FASTA_PATH = 'exacresiduals/data/hg19.fa'

ys, genes = [], []

fasta = Fasta(FASTA_PATH, as_raw=True)

def syn_cpg_density(pairs, d, gnomad, kcsq):
    syn=0
    nbases=0
    prevvar=None
    for pair in pairs:
        r0=str(int(pair[0])+1); r1=str(int(pair[1]))
        for v in gnomad(d['chrom']+':'+r0+'-'+r1):
            if prevvar is not None and str(v.start)+str(v.end)+str(v.ALT[0])==prevvar: continue
            if not (v.FILTER is None or v.FILTER in ["PASS", "SEGDUP", "LCR"]):
                continue
            info = v.INFO
            try:
                as_filter=info['AS_FilterStatus'].split(",")[0]
                if as_filter not in ["PASS", "SEGDUP", "LCR"] :
                    continue
            except KeyError:
                pass
            info = v.INFO
            try:
                csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
            except KeyError:
                continue
            for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
                if csq['Feature'] == '' or csq['EXON'] == '' or csq['cDNA_position'] == '' or csq['SYMBOL']!=d['gene']: continue #in case non-exonic or not the same gene
                if u.issynonymous(csq) and str(v.REF) == 'C' and str(v.ALT[0]) == 'T':
                    fa = str(fasta[d['chrom']][v.start+2]) # we want the 1-based position after this
                    if fa == 'G':
                        syn+=1
                        
            prevvar=str(v.start)+str(v.end)+str(v.ALT[0])

    return syn

cpg = []; syn = []
for i, d in enumerate(ts.reader(rfile)):
    if d['chrom'] in ['X', 'Y']: continue
    pairs = [x.split("-") for x in d['ranges'].strip().split(",")]
    syn_ct=syn_cpg_density(pairs, d, gnomad, kcsq)
    if int(d['n_bases'])>20:
        num_c_in_cpg = float(d['cg_content'])*float(d['n_bases'])/2.0
        try:
            syn.append(syn_ct/num_c_in_cpg)
        except ZeroDivisionError:
            syn.append(0)
        cpg.append(float(d['cg_content']))

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import pearsonr


outfile = 'syncpg.pdf'

fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(cpg, syn, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(1 if round(10**x-1,-1) == 0 else round(10**x-1,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(cpg, syn)
print pr, pval
plt.tight_layout()
sns.despine()
print "Pearson's r: " + str(pr) + "\np-value: " + str(pval)
plt.xlabel('CpG density')
plt.ylabel('Synonymous density')
plt.savefig(outfile, format='pdf', bbox_inches='tight')
