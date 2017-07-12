#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from optparse import OptionParser
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt

parser = OptionParser()

parser.add_option("-t",
                  "--title",
                  dest="title",
                  help="Title")

parser.add_option("-l",
                  "--labels",
                  dest="labels",
                  nargs=2,
                  help="data labels")

parser.add_option("-x",
                  "--xlabel",
                  dest="xlabel",
                  help="X axis label")

parser.add_option("-y",
                  "--ylabel",
                  dest="ylabel",
                  help="Y axis label")

parser.add_option("-o",
                  "--output_file",
                  dest="output_file",
                  help="Data file")

parser.add_option("--x_max",
                  dest="max_x",
                  type="float",
                  help="Max x value")

parser.add_option("--x_min",
                  dest="min_x",
                  type="float",
                  help="Min x value")

parser.add_option("--y_max",
                  dest="max_y",
                  type="float",
                  help="Max y value")

parser.add_option("--y_min",
                  dest="min_y",
                  type="float",
                  help="Min y value")

parser.add_option("--delim",
                  dest="delim",
                  default='\t',
                  help="delimiter")

(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')
Y,Y2=[],[]
for l in sys.stdin:
    a = l.rstrip().split(options.delim)
    if len(a) == 2:
        try:
            Y.append(float(a[0]))
        except:
            pass
        try:
            Y2.append(float(a[1]))
        except:
            pass
    if len(a) == 1:
        Y.append(float(a[0]))

if options.max_x:
    x_max = options.max_x
if options.min_x:
    x_min = options.min_x

if ((options.max_y) and (options.min_y)):
    ax.set_ylim(float(options.min_y),float(options.max_y))

matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure(figsize=(5,5)) # adjust figsize to change shape of plot and proportions 
fig.subplots_adjust(hspace=.05, wspace=.05, left=.01, bottom=.01)
ax = fig.add_subplot(1,1,1)

import seaborn # if imported earlier, fails to draw spines

if not Y2:
    mi = np.min(Y)
    ma = np.max(Y)
    rng = (mi,ma)
    p,p_edges=np.histogram(Y, bins=40, range=rng)
    p=map(lambda x: float(x)/sum(p), p)
    p.append(p[-1])
    width_p = (p_edges[1]-p_edges[0])
    print p, p_edges
    ax.plot(p_edges, p, color = 'crimson', ls = 'steps', lw=2, label = 'pathogenic')
    #ax.bar(p_edges[:-1], p, width = width_p, color = 'navy', label = 'pathogenic')
    skew=ss.skew(Y)
    ax.legend(loc='upper right')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('CCR Percentile (higher = more constrained)')#ax.set_xlabel('CpG')
    if options.title:
        title=options.title
        plt.title(title+"\nskewness: "+str(skew)+[" (negative)"," (positive)"][int(skew>=0)]) #negative skew implies more pathogenic
    ax.set_xlim(-(ma-mi)/100, ma)
if Y2:
    if options.labels:
        labels=options.labels
        label1=labels[0]
        label2=labels[1]
    else:
        label1="benign"
        label2="pathogenic"
    mi = np.min(Y+Y2)
    if options.min_x:
        mi = options.min_x
    ma = np.max(Y+Y2)
    rng = (mi,ma)
    p,p_edges=np.histogram(Y2, bins=40, range=rng)
    p=map(lambda x: float(x)/sum(p), p)
    p.append(p[-1])
    width_p = (p_edges[1]-p_edges[0])
    print "patho", p, p_edges
    #ax.bar(p_edges[:-1], p, width = width_p, color = 'g', label = 'KBM7: CS <- 2, p < 1E-3', alpha = 0.3)
    b,b_edges=np.histogram(Y, bins=40, range=rng)
    b=map(lambda x: float(x)/sum(b), b)
    b.append(b[-1])
    width_b = (b_edges[1]-b_edges[0])
    print "benign", b, b_edges
    #ax.bar(b_edges[:-1], b, width = width_b, color = 'r', label = 'random', alpha = 0.3)
    ax.spines['left'].set_color('k')
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_color('k')
    ax.spines['bottom'].set_visible(True)
    ax.plot(p_edges, p, label = label2, color = 'crimson', ls = 'steps', lw=2, alpha=0.7)#Essential Genes: CS < -1, p-val < 1E-3
    ax.plot(b_edges, b, label = label1, color = 'deepskyblue', ls = 'steps', lw=2, alpha=0.7)#random
    #ax.legend(loc='upper left')
    ax.legend(loc='upper right')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('CCR Percentile (higher = more constrained)')
    #ax.set_xlabel('MPC')
    ustat,upval=ss.mannwhitneyu(Y2,Y,use_continuity=False,alternative='greater') #pathogenic should be greater than benign
    kstat,kpval=ss.ks_2samp(np.array(Y2),np.array(Y))
    upval = '%.1e' % upval; ustat = '%.2e' % ustat; kpval = '%.1e' % kpval; kstat = str(round(kstat,3));
    if options.title:
        title=options.title
        plt.title(title+"\nU-test: "+ustat+" p-value: "+upval+"; ks-test: "+kstat+" p-value: "+kpval)
    ax.set_xlim(-(ma-mi)/100, ma)

if options.max_x:
    ax.set_xlim(xmax=options.max_x)
if options.min_x:
    ax.set_xlim(xmin=options.min_x)
if options.max_y:
    ax.set_ylim(ymax=options.max_y)
if options.min_y:
    ax.set_ylim(ymin=options.min_y)

#if options.logy:
    #ax.set_yscale('log')

if options.xlabel:
    ax.set_xlabel(options.xlabel)

if options.ylabel:
    ax.set_ylabel(options.ylabel)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#if options.black:
#    ax.spines['bottom'].set_color('white')
#    ax.spines['left'].set_color('white')
#    ax.title.set_color('white')
#    ax.yaxis.label.set_color('white')
#    ax.xaxis.label.set_color('white')
#    ax.tick_params(axis='x', colors='white')
#    ax.tick_params(axis='y', colors='white')

plt.savefig(options.output_file,bbox_inches='tight')
