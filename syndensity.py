from cyvcf2 import VCF
import utils as u
import gzip
from pyfaidx import Fasta
import toolshed as ts
import sys
from itertools import groupby
from operator import itemgetter

# script for calculating density of C-T synonymous variants in CpGs only versus CpG density

gnomad=VCF('exacresiduals/data/gnomad-vep-vt.vcf.gz') 
kcsq = gnomad["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
rfile = 'exacresiduals/results/gnomAD10x.5syn/exac-regions-novariant.txt'
#rfile = 'test.txt'
FASTA_PATH = 'exacresiduals/data/hg19.fa'

ys, genes = [], []
fasta = Fasta(FASTA_PATH, as_raw=True)

def isfunctional(csqs):
    for csq in csqs.split(","):
        eff = csq.split("|", 2)[0]
        for c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon', 'rare_amino_acid',
                     'missense', 'protein_altering', 'frameshift', 'inframe_insertion', 'inframe_deletion'):
            if c in eff or (('splice_donor' in eff or 'splice_acceptor' in eff) and 'coding_sequence' in eff):
                return True
    return False

def syn_cpg_density(pairs, d, gnomad, kcsq): #anything but functional, since by definition it would be synonymous
    ct_vars = 0
    all_vars = 0
    prevvar=None
    for pair in pairs:
        r0=str(int(pair[0])+1); r1=str(int(pair[1]))
        for v in gnomad(d['chrom']+':'+r0+'-'+r1):
            if prevvar is not None and str(v.start)+str(v.end)+str(v.ALT[0])==prevvar: continue
            if not (v.FILTER is None or v.FILTER in ["PASS", "SEGDUP", "LCR", "RF"]): # may remove RF
                continue
            info = v.INFO
            try:
                as_filter=info['AS_FilterStatus'].split(",")[0]
                if as_filter not in ["PASS", "SEGDUP", "LCR", "RF"]: # may remove RF
                    continue
            except KeyError:
                pass
            info = v.INFO
            try:
                csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
            except KeyError:
                continue
            if isfunctional(info['CSQ']): continue
            if str(v.REF) == 'C' and str(v.ALT[0]) == 'T':
                fa = str(fasta[d['chrom']][v.start+1]) # we want the 0-based position after this
                if fa == 'G': ct_vars += 1
            elif str(v.REF) == 'G' and str(v.ALT[0]) == 'A':
                fa = str(fasta[d['chrom']][v.start-1]) # we want the 0-based position before this
                if fa == 'C': ct_vars +=1
            all_vars += 1
            prevvar=str(v.start)+str(v.end)+str(v.ALT[0])
    return ct_vars, all_vars

rangeprev = None
sorter = itemgetter('chrom','gene','ranges')
grouper = itemgetter('chrom','gene','ranges')
ccrtemp = []
for ccr in ts.reader(rfile):
    if ccr['chrom'] in ['X', 'Y']: continue
    if int(ccr['n_bases'])<20 or float(ccr['cg_content'])==0.0: continue
    ccrtemp.append(ccr)
cpg = []
ctd = []
alld = []
for key, grp in groupby(sorted(ccrtemp, key = sorter), grouper):
    grp = list(grp)
    ranges = grp[0]['ranges']
    pairs = [x.split("-") for x in ranges.strip().split(",")]
    ct_cnt, all_cnt = syn_cpg_density(pairs, grp[0], gnomad, kcsq)
    try:
        ct_density = ct_cnt/float(grp[0]['n_bases'])
        all_density = all_cnt/float(grp[0]['n_bases'])
    except ZeroDivisionError:
        continue
    print '\t'.join([grp[0]['chrom'], grp[0]['start'], grp[0]['end'], grp[0]['cg_content'], str(int(float(grp[0]['cg_content'])*float(grp[0]['n_bases']))), str(ct_cnt), str(all_cnt), str(ct_density), str(all_density)])
    cpg.append(float(grp[0]['cg_content']))
    ctd.append(ct_density)
    alld.append(all_density)

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import pearsonr

outfile = '/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/syncpg'

fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(cpg, ctd, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(round(10**x) if round(10**x,-1) == 0 else round(10**x,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(cpg, ctd)
print "Pearson's r for C->T: " + str(pr) + "\np-value: " + str(pval)

plt.tight_layout()
sns.despine()
plt.xlabel('CpG density')
plt.ylabel('C->T synonymous variant density')
plt.savefig(outfile+'.pdf', format='pdf', bbox_inches='tight')

# all variants, not C-T/G-A
fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(cpg, alld, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(round(10**x) if round(10**x,-1) == 0 else round(10**x,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(cpg, alld)
print "Pearson's r for all: " + str(pr) + "\np-value: " + str(pval)

plt.tight_layout()
sns.despine()
plt.xlabel('CpG density')
plt.ylabel('All synonymous variant density')
plt.savefig(outfile+'all.pdf', format='pdf', bbox_inches='tight')
