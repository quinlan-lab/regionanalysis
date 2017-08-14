import cPickle as pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("white")
import sys

filename=sys.argv[1]

scored = pickle.load(open("scoredict.pkl", "rb"))
patho=np.array(scored["gnomad10x5_ccr"][1])
benign=np.array(scored["gnomad10x5_ccr"][0])

def step_plot(patho, benign, ax, ax2, **kwargs):
    p, p_edges = np.histogram(patho, bins=kwargs.pop('bins', 50), range=[patho.min(), patho.max()], normed=True)
    p=map(lambda x: float(x)/sum(p), p)
    p = list(p)
    p.append(p[-1])
    ax.plot(p_edges, p, ls='steps', lw=1.9, label="pathogenic", **kwargs)
    ax2.plot(p_edges, p, ls='steps', lw=1.9, label="pathogenic", **kwargs)
    b, b_edges = np.histogram(benign, bins=kwargs.pop('bins', 50), range=[benign.min(), benign.max()], normed=True)
    b=map(lambda x: float(x)/sum(b), b)
    b = list(b)
    b.append(b[-1])
    ax.plot(b_edges, b, ls='steps', lw=1.9, label="benign", **kwargs)
    ax2.plot(b_edges, b, ls='steps', lw=1.9, label="benign", **kwargs)
    ax.set_ylim(max(b)-.2, 1.)  # outliers only
    ax2.set_ylim(0, max(p))  # most of the data
    sns.despine(bottom=True)
    ax2.spines['bottom'].set_visible(True)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()

# get the red and blue colors for path, benign
sns.set_palette(sns.color_palette("Set1",10))
fig, (ax, ax2) = plt.subplots(2,1, sharex=True)
step_plot(patho, benign, ax, ax2, alpha=0.85)

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

ax2.set_xlabel("CCR Percentile") #gnomad10x5_ccr
rng = patho.max() - patho.min()
ax.set_xlim(patho.min() - 0.01 * rng, patho.max() + 0.01 * rng)
ax2.set_ylabel("Frequency")
ax2.yaxis.set_label_coords(-0.1,1)
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/fig2'+filename+'.pdf', bbox_inches='tight')
