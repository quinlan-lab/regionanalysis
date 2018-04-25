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

outfile='exoncpgvars'
GTF_PATH='exacresiduals/data/Homo_sapiens.GRCh37.75.gtf.gz'
VCF_PATH='exacresiduals/data/gnomad-vep-vt.vcf.gz'
COVERAGE_PATH='exacresiduals/data/exacv2.chr{chrom}.cov.txt.gz'
FASTA_PATH='exacresiduals/data/hg19.fa'
SELF_CHAINS = "exacresiduals/data/self-chains.gt90.bed.gz"
SEGDUPS = "exacresiduals/data/segmental.bed.gz"
exclude = [SELF_CHAINS, SEGDUPS]
chroms = [str(x) for x in range(1, 23)] + ["X", "Y"]
exac = VCF(VCF_PATH)
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

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
            if end - start < 4: continue
            mranges.append((start, end)); varflags.append('VARFALSE')
        mranges2, varflags2 = u.split_ranges(mranges, splitter, varflags)
        for ranges, vf in zip(mranges2, varflags2):
            seqs = [fa[s:e] for s, e in ranges]
            for (s, e), seq in zip(ranges, seqs):
                cg = u.floatfmt(u.cg_content(seq))
                length = float(len(seq))
                viter = VCF(VCF_PATH)(chrom+":"+str(s)+"-"+str(e))
                ct = 0
                for v in viter:
                    if not (v.FILTER is None or v.FILTER in ["PASS", "SEGDUP", "LCR"]):
                        continue
                    info = v.INFO
                    try:
                        as_filter=info['AS_FilterStatus'].split(",")[0]
                        if as_filter not in ["PASS", "SEGDUP", "LCR"] :
                            continue
                    except KeyError:
                        pass
                    try:
                        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
                    except KeyError:
                        continue
                    # NOTE: using max here for alternates to be conservative
                    try: # gnomad doesn't have adj like exacv1
                        ac = info['AC_Adj']
                    except KeyError:
                        ac = info['AC']
                    if not isinstance(ac, (int, long)):
                        ac = max(ac)
                    try:
                        af = ac / float(info['AN_Adj'] or 1)
                    except KeyError:
                        af = ac / float(info['AN'] or 1)
                    for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'): # getting duplicate rows because of this, wastes memory and potentially compute time, could remove and replace with just if isfunctional, add to rows then move on?
                        # skipping intronic
                        if csq['Feature'] == '' or csq['EXON'] == '': continue #or csq['cDNA_position'] == '': continue
                        if not u.ismissense(csq) or not u.issynonymous: continue
                        ct += 1
                        break
                outs.append((float(cg), ct, length))

    return outs

import multiprocessing as mp
p = mp.Pool(12)

for outs in p.imap_unordered(perchrom, ((VCF, str(chrom)) for chrom in chroms)):
    for out in outs:
        print("\t".join(map(str,out)))

cpg = [i[0] for i in outs]
ct = [i[1] for i in outs]
length = [i[2] for i in outs]
fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
g=ax.hexbin(cpg, ct, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(1 if round(10**x-1,-1) == 0 else round(10**x-1,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
counts,edges=np.histogram(g.get_array(),bins=8)
cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
cbar.set_label('Number of Regions', rotation=270, labelpad=20)
#pearson correlation
pr, pval = pearsonr(cpg, ct)
print (pr, pval)
plt.tight_layout()
sns.despine()
print ("Pearson's r: " + str(pr) + "\np-value: " + str(pval))
plt.xlabel('CpG density')
plt.ylabel('Variant count')
plt.savefig(outfile+'.pdf', format='pdf', bbox_inches='tight')

#density of variant plot
fig, ax = plt.subplots()
colors = sns.color_palette("GnBu", 10)
cmap_name='mymap'
cm = LinearSegmentedColormap.from_list(cmap_name, colors)
dens = [float(c)/float(l) for c, l in zip(ct,length)]
g=ax.hexbin(cpg, dens, cmap=cm, bins='log', alpha=0.5, mincnt=1)
def y_format(x,y):
    return '{:.0f}'.format(1 if round(10**x-1,-1) == 0 else round(10**x-1,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
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
plt.ylabel('Variant density')
plt.savefig(outfile+'density.pdf', format='pdf', bbox_inches='tight')
