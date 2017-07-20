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
from matplotlib import gridspec
import seaborn as sns
from matplotlib import collections as mc
sns.set_style('white')


def geneplot(exons, patho_variants, population_variants=None, constraint=None,
        density=None,
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
    fig = plt.figure(figsize=(20, 6))
    height_ratios =  (1, 1) if density is None else (0.8, 0.8, 1)
    sgs = gridspec.GridSpec(2, 1, height_ratios=[6, 1], hspace=0.0)
    gs = gridspec.GridSpecFromSubplotSpec(2 if density is None else 3, len(exons), width_ratios=widths,
            height_ratios=height_ratios, hspace=0.2, subplot_spec=sgs[0])
    gs2 = gridspec.GridSpecFromSubplotSpec(1, len(exons), width_ratios=widths, subplot_spec=sgs[1])

    ncols = len(exons)
    ax0 = None


    for i, exon in enumerate(exons):
        if i == 0:
            ax = plt.subplot(gs[0, i])
            ax2 = plt.subplot(gs[0, i])
            ax0 = ax
            ax2_0 = ax2
        else:
            ax = fig.add_subplot(gs[0, i], sharey=ax0)
            ax2 = fig.add_subplot(gs2[0, i], sharey=ax2_0)

        vs = [v for v in patho_variants if exon[0] <= v[0] <= exon[1]]
        pop = [v for v in population_variants if exon[0] <= v[0] <= exon[1]]
        ctr = [v for v in constraint if exon[0] <= v[0] <= exon[1]]

        xs, ys = [], [] # line width controls height of heatmap

        lines = []
        colors = []
        for s, e, height in ctr:
            lines.append([(s, 0), (e, 0)])
            if height < 80:
                colors.append((0, 0, 1, height/100))
            elif height < 100:
                colors.append((1, 0, 0, height/100))
            #xs.extend([s, e])
            #ys.extend([height, height])
        lc = mc.LineCollection(lines, colors=colors, linewidths=300) # line width controls height of heatmap
        ax.set_ylim(0,20)
        ax.add_collection(lc)
        #ax.step(xs, ys, color=opts['constraint_color'])
        #ax.plot([s, e], [height, height], color=opts['constraint_color'])

        axv = fig.add_subplot(gs[1 if density is None else 2, i], sharex=ax)
        #axv.set_ylim(ymax=2)

        if len(vs) > 0:
            markers, stemlines, baseline = axv.stem([v[0] for v in vs], -np.log10([v[1]
                for v in vs]), linefmt='-', markerfmt='o')
            plt.setp(baseline, 'linewidth', 0)
            plt.setp(stemlines, 'color', opts['patho_variant_color'], 'zorder', 0, 'alpha', 0.7)

            plt.setp(markers, 'color', opts['patho_variant_color'], 'zorder', 2,
                    'alpha', 0.7, 'markeredgecolor', '#666666', 'mew', 1,
                    'markersize', 5)

        #axe = fig.add_subplot(gs[2, i], sharex=ax)
        if len(pop) > 0:
            markers, stemlines, baseline = axv.stem([v[0] for v in pop], -np.log10([v[1] for v in pop]), '--' )
            plt.setp(baseline, 'linewidth', 0)
            plt.setp(stemlines, 'color', opts['pop_variant_color'], 'zorder', -1, 'alpha', 0.7)
            plt.setp(markers, 'color', opts['pop_variant_color'], 'zorder', 1,
                    'alpha', 0.7, 'markeredgecolor', '#666666', 'mew', 1,
                    'markersize', 4)

        if i == 0:
            ax.set_ylabel('Constraint')
            axv.set_ylabel("-log10(AF)") # common and pathogenics
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.setp(axv.get_yticklabels(), visible=False)
        ax.set_yticks([])
        axv.set_yticks([])
        axv.set_xticks([])

        axe = fig.add_subplot(gs2[0, i], sharex=ax2)
        axe.set_yticks([])
        axe.set_xticks([])
        axe.set_ylim(0,1)
        #if i == 0:
        #    axe.set_ylabel('Exons')
        axe.axhspan(.5, 1, xmin=0.001, xmax=0.999, facecolor='orange', edgecolor=opts['exon_color'],
                lw=4, zorder=9) # zorder makes sure it's always on top

        if not density:
            continue
        idensity = np.zeros(exon[1] - exon[0])
        for d in density:
            if d >= exon[0] and d < exon[1]:
                idensity[d-exon[0]]+=1
        idensity = np.convolve(idensity, np.ones(opts['density_window']) / float(opts['density_window']), mode='same')
        axd = fig.add_subplot(gs[1, i], sharex=ax)
        axd.plot(np.arange(exon[0], exon[1]), idensity)
        if i == 0:
            axd.set_ylabel('Variant density')
        axd.set_yticks([])
        axd.set_xticks([])

    sns.despine(left=True, bottom=True)
    #gs.tight_layout(fig, h_pad=0)
    #plt.tight_layout()
    #return fig, gs
    plt.savefig('proplot.pdf', bbox_inches='tight')


def get_control_genome_positions(control_vcf, region, query_transcript):
    control_positions = []
    densities = []
    for v in control_vcf(region):
        if v.var_type != 'snp': continue
        effects = [e.split('|') for e in v.INFO.get('CSQ').split(',')]
        for e in effects:
            curr_transcript = e[6]
            impact = e[1]
            if impact not in ['missense_variant']:
                continue
            if v.INFO['AC'] < 1: continue
            if curr_transcript == query_transcript:
                # TODO: track actual AAF
                control_positions.append((v.POS, v.INFO['AF']))
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
        #if val != 0: ccrs.append((s,e,val))
        ccrs.append((s,e,val))
    return ccrs



def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-b", "--ccr", dest='ccr', help="CCR file")
    parser.add_argument("-c", "--control-vcf", dest='control_vcf', help="VCF file of control variants")
    parser.add_argument("-p", "--path-vcf", dest='patho_vcf', help="VCF file of pathogenic variants")
    parser.add_argument("-g", "--gff", dest='gff', help="GFF of gene models")
    parser.add_argument("-t", "--transcript", dest='transcript', help="Transcript to plot")
    parser.add_argument("-r", "--region", dest='region', help="Region to query")

    args=parser.parse_args()

    control_vcf = VCF(args.control_vcf)
    patho_vcf = VCF(args.patho_vcf)

    population_variants, density   = get_control_genome_positions(control_vcf, args.region, args.transcript)
    patho_variants = get_patho_genome_positions(patho_vcf, args.region)
    exons = get_exons(args.gff, args.transcript, args.region)
    ccrs = get_ccrs(args.ccr, args.region)
    geneplot(exons, patho_variants, population_variants, constraint=ccrs, density=density)

if __name__ == "__main__":
    main()
