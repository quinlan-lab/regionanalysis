import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import cPickle as pickle
import sys
import seaborn as sns
sns.set_style('white')
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches
import numpy as np

infile = sys.argv[1]
pk = open(infile, 'rb')
data = pickle.load(pk)
pk.close()
if infile == "ccrgerppfam.pkl":
    cscores, gscores, labels = data

    plt.scatter(cscores,gscores,marker='o',s=1)
    sns.despine()
    for label, c, g in zip(labels, cscores, gscores):
        if (c > 90 and g > 4) or (g < 0 and c < 20):
            plt.annotate(
                label, alpha=0.6,
                xy=(c, g), xytext=(-25, 25),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.3),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'),#arc3,rad=0'),
                fontsize=5)
    plt.xlabel("CCR")
    plt.ylabel("GERP")
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/gerpvccr_pfam.pdf', bbox_inches='tight')
if infile == "ccrgerp.pkl":
    ct=0; highct=0; total=0
    gerp, ccr, genes, lengths, ranges, chroms = zip(*data)
    constgenes=[]
    for i, j, k in zip(gerp, ccr, genes):
        if j >= 95:
            highct+=1
            if i < 1.7: 
                ct+=1
                constgenes.append(k)
        total+=1
    print ct, highct, total
    genefile=open('geneset.txt','w')
    for gene in constgenes:
        genefile.write(gene+"\n")
    fig, ax = plt.subplots()
    ax.add_patch(patches.Rectangle((95, min(gerp)), 5, abs(1.7-min(gerp)), fill=False, edgecolor='red', linewidth=2))
    colors = [(.9, .9, .9), (0.627451, 0.12549, 0.941176), (0.254902, 0.411765, 0.882353), (0.603922, 0.803922, 0.196078), (1, 0.647059, 0,), (1, 0, 0)]
    cmap_name='mymap'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors)
    g=ax.hexbin(ccr, gerp, cmap=cm, bins='log', alpha=0.5, mincnt=1)
    def y_format(x,y):
        return '{:.0f}'.format(1 if round(10**x-1,-1) == 0 else round(10**x-1,-1)) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
        #return '{:.0f}'.format(round(10**x-1,-1) if 10**x-1 != 1 else 1) # not 100% accurate binning, but the -1 is so we can label the bottom of the colorbar as 0, doesn't throw off calc by much
        
    counts,edges=np.histogram(g.get_array(),bins=8)
    cbar = fig.colorbar(g, ax=ax, orientation='vertical', extend='both', extendrect=True, drawedges=False, ticks=edges, format=FuncFormatter(y_format))
    cbar.set_label('Number of Regions', rotation=270, labelpad=20)
    plt.tight_layout()
    sns.despine()
    plt.xlabel("CCR\n"+str(ct)+" out of "+str(total)+" unique regions in red box")
    plt.ylabel("Mean GERP++ score")
    plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/gerpvccr.pdf', bbox_inches='tight')
  
    # making table for supplement
    f=open('purifyingselectionregions.tsv','w')
    f.write("Chrom\tStart\tEnd\tGene\tRanges\tLength\tCCR Percentile Score\tMean GERP++ score\n")
    for rangegroup, gerp, ccr, gene, length, chrom in zip(ranges, gerp, ccr, genes, lengths, chroms):
        if ccr >= 95 and gerp < 1.7:
            for range in rangegroup.split(","):
                f.write("\t".join([chrom] + range.split("-") + [gene, rangegroup, str(length), str(ccr), str(gerp)])+"\n")
    f.close()
