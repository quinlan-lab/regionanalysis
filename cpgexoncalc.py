# script for calculating variant count in flattened exons (filtered on coverage and segdups like CCRs) vs CpG density
# can make it CpG specific mutations
from __future__ import print_function
import sys
import exacresiduals.utils as u
from cyvcf2 import VCF
from pyfaidx import Fasta
import itertools as it
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import pearsonr
import os.path

outfile=os.path.expandvars('$HOME/public_html/randomplots/exoncpg')
GTF_PATH='exacresiduals/data/Homo_sapiens.GRCh37.75.gtf.gz'
VCF_PATH='exacresiduals/data/gnomad-vep-vt.vcf.gz'
COVERAGE_PATH='exacresiduals/data/exacv2.chr{chrom}.cov.txt.gz'
FASTA_PATH='exacresiduals/data/hg19.fa'
SELF_CHAINS = "exacresiduals/data/self-chains.gt90.bed.gz"
SEGDUPS = "exacresiduals/data/segmental.bed.gz"
exclude = [SELF_CHAINS, SEGDUPS]
chroms = [str(x) for x in range(1, 23)] # + ["X", "Y"]
gnomad = VCF(VCF_PATH)
kcsq = gnomad["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

def syn_cpg_density(fasta, gnomadranges, kcsq): #anything but functional, since by definition it would be synonymous
    ct_vars = 0
    all_vars = 0
    prevvar=None
    for v in gnomadranges:
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
        if str(v.REF) == 'C' and str(v.ALT[0]) == 'T':
            fa = str(fasta[v.start+1]) # we want the 0-based position after this
            if fa == 'G': ct_vars += 1
        elif str(v.REF) == 'G' and str(v.ALT[0]) == 'A':
            fa = str(fasta[v.start-1]) # we want the 0-based position before this
            if fa == 'C': ct_vars +=1
        all_vars += 1
        prevvar=str(v.start)+str(v.end)+str(v.ALT[0])
    return ct_vars, all_vars

def perchrom(vcf_chrom):
    vcf, chrom = vcf_chrom

    chrom=str(chrom)
    rows = []
    outs = []
    print("reading chrom " + chrom, file=sys.stderr)

    fasta = Fasta(FASTA_PATH, as_raw=True)
    fa = str(fasta[chrom])
    coverage_array = u.read_coverage(chrom, length=len(fa), cov=10,
                        path=COVERAGE_PATH)

    gene_exon_starts, gene_exon_ends, splitters = u.read_exons("|tabix {gtf} {chrom}"
                                                            .format(chrom=chrom,
                                                                gtf=GTF_PATH),
                                                            chrom, 0.5, #cutoff is 0.5
                                                            coverage_array,
                                                            exclude)

    for chrom_gene in gene_exon_starts:
        exon_starts = gene_exon_starts[chrom_gene]
        exon_ends = gene_exon_ends[chrom_gene]
        last = exon_starts[0]
        splitter = splitters.get(chrom_gene, None)
        mranges = []; varflags = []
        for start, end in zip(exon_starts, exon_ends):
            mranges.append((start, end)); varflags.append('VARFALSE')
        mranges2, varflags2 = u.split_ranges(mranges, splitter, varflags)
        for ranges, vf in zip(mranges2, varflags2):
            seqs = [fa[s:e] for s, e in ranges]
            for (s, e), seq in zip(ranges, seqs):
                if e - s < 2: continue
                cg = u.cg_content(seq)
                length = float(len(seq))
                viter = VCF(VCF_PATH)(chrom+":"+str(s+1)+"-"+str(e)) # vcf is 1-based
                ct_vars, all_vars = syn_cpg_density(fa, viter, kcsq)
                if all_vars/length >=2: print (ranges, chrom_gene) # post-hoc sanity check
                outs.append((cg, ct_vars, all_vars, length))

    return outs

import multiprocessing as mp
p = mp.Pool(22)

output = []
for outs in p.imap_unordered(perchrom, ((VCF, str(chrom)) for chrom in chroms)):
    output.append(outs)
        

cpg = [i[0] for i in outs]
ct = [i[1] for i in outs]
all_ct = [i[2] for i in outs]
length = [i[3] for i in outs]

#C-T density of variants plot
fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
dens = [float(c)/float(l) for c, l in zip(ct,length)]
g=ax.hexbin(cpg, dens, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(round(10**x) if round(10**x,-1) <= 100 else round(10**x,-1)) # more round number binning
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(cpg, dens)
print (pr, pval)
plt.tight_layout()
sns.despine()
print ("Pearson's r: " + str(pr) + "\np-value: " + str(pval))
plt.xlabel('CpG density')
plt.ylabel('C->T variant density')
plt.savefig(outfile+'C-Tvariantdensity.pdf', format='pdf', bbox_inches='tight')

fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
dens = [float(c)/float(l) for c, l in zip(all_ct,length)]
g=ax.hexbin(cpg, dens, cmap=cm, bins='log', alpha=0.5, mincnt=1)
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(cpg, dens)
print (pr, pval)
plt.tight_layout()
sns.despine()
print ("Pearson's r: " + str(pr) + "\np-value: " + str(pval))
plt.xlabel('CpG density')
plt.ylabel('All variant density')
plt.savefig(outfile+'allvariantdensity.pdf', format='pdf', bbox_inches='tight')
