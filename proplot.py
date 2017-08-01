import sys
from cyvcf2 import VCF
import argparse
from geneimpacts import VEP
from collections import defaultdict
import tabix
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, transforms
import seaborn as sns
from matplotlib import collections as mc
from matplotlib.colors import ListedColormap, BoundaryNorm
sns.set_style('white')
import doctest

def rainbow_text(x, y, strings, colors, ax=None, **kw):
    """
    Take a list of ``strings`` and ``colors`` and place them next to each
    other, with text strings[i] being shown in colors[i].

    This example shows how to do both vertical and horizontal text, and will
    pass all keyword arguments to plt.text, so you can set the font size,
    family, etc.

    The text will get added to the ``ax`` axes, if provided, otherwise the
    currently active axes will be used.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    # horizontal version
    for s, c in zip(strings, colors):
        text = ax.text(x, y, s + " ", color=c, transform=t, **kw)
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(text._transform, x=ex.width, units='dots')

def overlaps(s1, e1, s2, e2):
    """
    >>> overlaps(2, 4, 3, 5)
    True
    >>> overlaps(2, 4, 4, 5)
    False
    >>> overlaps(2, 200, 3, 5)
    True
    >>> overlaps(3, 5, 2, 200)
    True
    >>> overlaps (3, 5, 5, 6)
    False
    >>> overlaps (5, 6, 3, 5)
    False
    """
    return not (e1 <= s2 or s1 >= e2)

def get_pfam(pfam, transcript, region):
    tb = tabix.open(pfam)
    pfams = []
    for r in tb.querys(region):
        s = int(r[1])
        e = int(r[2])
        details = r[9]
        if transcript not in details: continue
        fam=details.split(";")[0].split('"')[1]
        pfams.append([s, e, fam])
    return pfams

def geneplot(exons, pfams, patho_variants, population_variants=None, constraint=None,
        density=None, coverage=None, repeats=None, filename=None,
        opts={'constraint_color': (0.7, 0.7, 0.7),
              'patho_variant_color': '#ff0000',
              'exon_color': (0.8,0.8, 0.8),
              'pop_variant_color': '#4daf4a',
              'density_window': 20,
              }):
    """
    >>> _, _ = geneplot(
    ... [[1234, 1298],[22222, 22349]],
    ... [(1244, 0.001), (1247, 0.1), (1255, 0.2), (22233, 0.001), (22244, 0.2)],
    ... [(1235, 0.1), (1236, 0.2), (22340, 0.1), (22341, 0.1), (22344, 0.1)],
    ... [(1234, 1244, 81), (1244, 1255, 99), (22233, 22341, 80), (22341, 22349, 10)],
    ... density=range(1244, 1290) + range(22244, 22255) + range(22261, 22269))
    >>> plt.show()
    """
    widths = [float(e[1] - e[0]) for e in exons]
    fig = plt.figure(figsize=(20, 2))
    #height_ratios = (1, 1)
    sgs = gridspec.GridSpec(3, 1, height_ratios=[0.5, 1, 1], hspace=0.0)
    gs = gridspec.GridSpecFromSubplotSpec(1, len(exons), width_ratios=widths, subplot_spec=sgs[0]) # 2, len
            #height_ratios=height_ratios, hspace=0.0)
    gs2 = gridspec.GridSpecFromSubplotSpec(2, len(exons), width_ratios=widths, subplot_spec=sgs[1], hspace=0.3) # 1, len
    gs3 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=sgs[2])

    doms=[]
    fams=set()
    colorlist = ["tomato", "dodgerblue", "darkorchid", "mediumvioletred", "coral", "black", "skyblue", "sienna", "gold"]
    colors={}; ct=0
    from matplotlib.ticker import FormatStrFormatter
    for i, exon in enumerate(exons):
        for j, domain in enumerate(pfams):
            if overlaps(exon[0], exon[1], domain[0], domain[1]):
                dom1 = (domain[0] if exon[0] < domain[0] else exon[0])
                dom2 = (domain[1] if exon[1] > domain[1] else exon[1])
                dom3 = domain[2]
                doms.append((dom1,dom2,dom3))
            
        ax_exon = fig.add_subplot(gs2[0, i])# sharex=ax_cons)
        ax_exon.set_xticks([])
        ax_exon.set_yticks([])
        ax_exon.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax_exon.set_ylim(0,1)
        ax_exon.set_xlim(exon[0],exon[1])
        #if i == 0:
        #    ax_exon.set_ylabel('Exons')
        ax_exon.axhspan(.6, 1, xmin=0.001, xmax=0.999, edgecolor=opts['exon_color'], facecolor = 'none',
            lw=1, zorder=10) # zorder makes sure it's always on top
        
        vs = [v for v in patho_variants if exon[0] <= v[0] <= exon[1]]
        pop = [v for v in population_variants if exon[0] <= v[0] <= exon[1]]
        ctr = [v for v in constraint if exon[0] <= v[0]+1 <= exon[1]] # GTF format (Exons) are 1-based, regions are in 0-based half-open BED format
        cov = [v for v in coverage if exon[0] <= v[0] <= exon[1]]
        
        xs, ys = [], [] # line width controls height of heatmap

        ax_cons = fig.add_subplot(gs[0, i],sharex=ax_exon)
        ax_cons.set_yticks([80,90,100])
        ax_cons.set_ylim(80,100) #80 is our low bar for constraint
        if len(pop) > 0:
            afs=[x[1] for x in pop]
            alphas=map(lambda x: 1.3--np.log10(x)/max(-np.log10([k for k in afs])), afs)
            alphas=[1 if k > 1 else k for k in alphas]
            for index, v in enumerate(pop):
                color = 'blue'; alpha=0.1; lw=1
                if v[0] == 62062708:
                    color='red'; alpha=1; lw=3; print exon[0], exon[1]
                ax_cons.axvline(x=v[0], ymin=0, ymax=100, lw=lw, color=color, alpha=alpha) #alpha=alphas[index])
                #markers, stemlines, baseline = axd.stem([v[0]], [100], linefmt='-', markerfmt=' ', lw=0.05, color='g', alpha=alphas[index])
                #plt.setp(baseline, 'linewidth', 0)
                #plt.setp(stemlines, 'color', opts['pop_variant_color'], 'zorder', -1) # may have to plot each stem in a loop to change alphas
                #plt.setp(markers, 'color', opts['pop_variant_color'], 'zorder', 1,
                #        'markeredgecolor', '#666666', 'mew', 1,
                #        'markersize', 4)

        for s, e, height in ctr:
            if height < 80: continue #only show constraint above our cutoff
            color = ('k' if height >= 80 else 'b')
            ax_cons.plot((s,e), (height,height), color=color)
        if i == 0:
            ax_cons.set_ylabel('Constraint')
        else:
            plt.setp(ax_cons.get_yticklabels(), visible=False)
            #plt.setp(ax_patho.get_yticklabels(), visible=False)
            ax_cons.set_yticks([])
        # ax_patho.set_yticks([])
        # ax_patho.set_xticks([])

        
        # ax_patho = fig.add_subplot(gs[1, i], sharex=ax_cons)

        if len(vs) > 0:
            #markers, stemlines, baseline = ax_patho.stem([v[0] for v in vs], -np.log10([v[1] for v in vs]), linefmt='-', markerfmt=' ', lw=0.01)
            for index, v in enumerate(vs):
                ax_exon.axvline(x=v[0], ymin=.62, ymax=1, color='k', lw=1, alpha=1, zorder=8)
            #plt.setp(baseline, 'linewidth', 0)
            #plt.setp(stemlines, 'color', 'black', 'zorder', 0, 'alpha', 0.7)

            #plt.setp(markers, 'color', opts['patho_variant_color'], 'zorder', 2,
            #        'alpha', 0.7, 'markeredgecolor', '#666666', 'mew', 1,
            #        'markersize', 4)
            #ax_patho.set_ylim(0,max(-np.log10([v[1] for v in vs]))+.5)

        for s, e, fam in doms:
            if not overlaps(s,e,exon[0],exon[1]):continue
            if fam not in fams:
                #colors[fam] = colorlist[ct]
                colors[fam] = 'lightgrey'
                ct+=1
                fams.add(fam) 
            xmin=(s-exon[0])/float(exon[1]-exon[0])
            xmax=1-(exon[1]-e)/float(exon[1]-exon[0])
            ax_exon.axhspan(.6, 1, xmin=xmin, xmax=xmax, edgecolor=opts['exon_color'], facecolor = colors[fam],
                lw=1, zorder=9, alpha=0.5) # zorder makes sure it's always on top
        #print colors.keys()
        #ax_leg.table(cellText=colors.keys(),cellColours=colors.values())

        #if not density:
        #    continue
        #idensity = np.zeros(exon[1] - exon[0])
        #for d in density:
        #    if d >= exon[0] and d < exon[1]:
        #        idensity[d-exon[0]]+=1
        #for j, d in enumerate(idensity):
        #    if i+1%opts['density_window'] == 0:
        #        idensity[i-opts['density_window']:i] = 0
        #idensity = np.convolve(idensity, np.ones(opts['density_window']) / float(opts['density_window']), mode='same')
        #axd = fig.add_subplot(gs[0, i], sharex=ax) #1,i
        #axd.plot(np.arange(exon[0], exon[1]), idensity)
        #if i == 0:
        #    axd.set_ylabel('Variant density')
        #axd.set_yticks([])
        #axd.set_xticks([])
        
        ax_coverage = fig.add_subplot(gs2[1, i], sharex=ax_cons)
        ax_coverage.plot([c[0] for c in cov], [c[1] for c in cov], color='g')
        ax_coverage.set_yticks([])
        ax_coverage.set_ylim(0,1)
    
        for s, e in repeats:
            if not overlaps(s,e,exon[0],exon[1]): continue
            ax_coverage.axhline(y=.5, xmin=s, xmax=e, lw=1.5, color='r')
    ax_leg = fig.add_subplot(gs3[0, 0]) # leg = legend
    ax_leg.set_ylim(0,1)
    ax_leg.set_xlim(0,1)
    rainbow_text(0,1,colors.keys(),colors.values(),ax=ax_leg, weight="semibold")#ax3.text(0.5,0.5,)
    ax_leg.set_yticks([])
    ax_leg.set_xticks([])
    sns.despine(left=True, bottom=True)
    #gs.tight_layout(fig, h_pad=0)
    #plt.tight_layout()
    #return fig, gs
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/' + filename + '.pdf', bbox_inches='tight')


def get_control_genome_positions(control_vcf, region, query_transcript):
    control_positions = []
    densities = []
    for v in control_vcf(region):
        if v.var_type != 'snp': continue
        effects = [e.split('|') for e in v.INFO.get('CSQ').split(',')]
        if v.FILTER is not None: continue
        for e in effects:
            curr_transcript = e[6]
            impact = e[1]
            if 'missense_variant' not in impact:
                continue
            if v.INFO['AC'] < 1: continue
            if curr_transcript == query_transcript:
                control_positions.append((v.POS, v.INFO['AF'])) # can do AC
                densities.append(v.POS)
    return control_positions, densities


def get_patho_genome_positions(patho_vcf, region):
    path_positions = []
    for v in patho_vcf(region):
        if len(v.REF) > 1: continue
        #if v.var_type != 'snp': continue
        #clinsig = v.INFO.get('CLNSIG')
        #if clinsig != '5': continue
        # can't track actual AAF for pathos
        path_positions.append((v.POS, 10E-4))
    return path_positions


def get_exons(gff, transcript, region):
    tb = tabix.open(gff)
    exons = []
    for r in tb.querys(region):
        t = r[2]
        s = int(r[3])
        e = int(r[4])
        details = r[8]

        if transcript not in details: continue
        if t != 'CDS': continue
        exons.append([s, e])
    return exons


def get_ccrs(ccr, region):
    tb = tabix.open(ccr)
    ccrs = [] 
    for r in tb.querys(region):
        s = int(r[1])
        e = int(r[2])
        val = float(r[13])
        ccrs.append((s,e,val))
    return ccrs

def get_coverage(region):
    chrom = region.split(":")[0]
    coverage = '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/data/exacv2.chr' + chrom + '.cov.txt.gz'
    tb = tabix.open(coverage)
    cov = [] 
    for r in tb.querys(region):
        pos = int(r[1])
        val = float(r[6]) # coverage at 10x
        if val < 0.5: # our CCR cutoff
            cov.append((pos,val))
    return cov

def get_repeats(region):
    chrom = region.split(":")[0]
    self = '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/data/self-chains.gt90.bed.gz'
    seg = '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/data/hgsegmental.bed.gz' 
    tb = tabix.open(self)
    rep = [] 
    for r in tb.querys(region):
        s = int(r[1])
        e = int(r[2])
        rep.append((s,e))
    tb = tabix.open(seg)
    for r in tb.querys(region):
        s = int(r[1])
        e = int(r[2])
        rep.append((s,e))
    return rep

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-b", "--ccr", dest='ccr', help="CCR file")
    parser.add_argument("-c", "--control-vcf", dest='control_vcf', help="VCF file of control variants")
    parser.add_argument("-p", "--path-vcf", dest='patho_vcf', help="VCF file of pathogenic variants")
    parser.add_argument("-g", "--gff", dest='gff', help="GFF of gene models")
    parser.add_argument("-d", "--pfam", dest='pfam', help="bed of Pfam domains (d)")
    parser.add_argument("-t", "--transcript", dest='transcript', help="Transcript to plot")
    parser.add_argument("-r", "--region", dest='region', help="Region to query")
    parser.add_argument("-f", "--filename", dest='filename', help="Name of output file")

    args=parser.parse_args()

    control_vcf = VCF(args.control_vcf)
    patho_vcf = VCF(args.patho_vcf)

    population_variants, density   = get_control_genome_positions(control_vcf, args.region, args.transcript)
    patho_variants = get_patho_genome_positions(patho_vcf, args.region)
    exons = get_exons(args.gff, args.transcript, args.region)
    ccrs = get_ccrs(args.ccr, args.region)
    pfams = get_pfam(args.pfam, args.transcript, args.region)
    cov = get_coverage(args.region)
    rep = get_repeats(args.region)
    geneplot(exons, pfams, patho_variants, population_variants, constraint=ccrs, density=density, coverage=cov, repeats=rep, filename=args.filename)

if __name__ == "__main__":
    main()
