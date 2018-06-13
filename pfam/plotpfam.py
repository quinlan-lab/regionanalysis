import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import numpy as np
import ast
import seaborn
seaborn.set_style(style='white')
import matplotlib.backends.backend_pdf
from scipy.stats import fisher_exact
import os.path

parser=ArgumentParser()
parser.add_argument("-p","--pfam", help="pfam families for histogram") # top200doms
parser.add_argument("-c","--clans", help="pfam clans file, has proper family names in it") # Pfam-A.clans.tsv.gz
parser.add_argument("-s","--stats", help="pfams bp covered in ccr %ile bins") #pfamshist.txt
parser.add_argument("-q","--count", help="occurrences of each pfam domain") #pfamcounts.txt
parser.add_argument("-o","--output", help="histogram pic") #pfam_hists.pdf
args=parser.parse_args()
table=open(os.path.dirname(os.path.expandvars(args.output))+'/pfam_stats(supp_table_3).tsv','w')
table.write('domain\tfull_name(if_applicable)\tnumber_of_occurrences\ttotal_bp\tbp_in_95_ccr_bin\tbinomial_p_val\n')

with open(args.pfam) as f:
    pfams=[_ for _ in f]
counts={}
with open(args.count) as f:
    for line in f:
        fields=line.strip().split()
        counts[fields[1]]=fields[0] # pfam id is key, true family name is value
families={}
with open(args.clans) as f:
    for line in f:
        fields=line.strip().split("\t")
        families[fields[3]]=fields[4] # pfam id is key, true family name is value

count=len(pfams)
stats=open(args.stats,"r")
bins=range(5,101,5)

def make_axes_visible(axes):
    for ax in axes.flatten():
        for tk in ax.get_yticklabels():
            tk.set_visible(True)
        for tk in ax.get_xticklabels():
            tk.set_visible(True)

ct=0; maxval=0
plt.rcParams["figure.figsize"]=(5,2.5)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
#f, axarr = plt.subplots(count, sharex=True)
#axarr = [plt.subplot(count, 1, i+1) for i in range(count)]

# to get counts for "c" and "d" for fisher's exact test.
# cells are the count of intervals
#
#                        |  in 95% CCR  | not in 95% CCR  |
#     ivls in domain     | a            | b               |
#     ivls not in domain | c            | d               |
all95 = 0.0; allnon95 = 0.0; domcnt = 0
for stat in stats:
    fields=stat.strip().split("\t")
    pfam2=fields[0].strip()
    domcnt += 1
    dbins=ast.literal_eval(fields[1]) #(dictionary of bins)
    if 100 in dbins:
        all95 += dbins[100][1]
    for i in dbins:
        if i != 100:
            allnon95 += dbins[i][1]
print all95, allnon95
print domcnt

sign = 0
for pfam in pfams:
    pfam=pfam.strip()
    for stat in stats:
        fields=stat.strip().split("\t")
        pfam2=fields[0].strip()
        dbins=ast.literal_eval(fields[1]) #(dictionary of bins)
        if pfam == pfam2:
            print pfam
            try:
                name = families[pfam]
                title = families[pfam] + "\n(" + pfam + ", N=" + counts[pfam] + ")"
            except KeyError:
                name = "N/A"
                title = "(" + pfam + ", N=" + counts[pfam] + ")"
            tccr=[[i,0] for i in bins]
            totlen=float(fields[2])
            totcnt=float(fields[3])
            if 100 in dbins:
                x = dbins[100][1]
            else:
                x = 0
            #pval = binom_test(x, totlen, 0.05, 'greater')
            a = x
            b = totcnt - x
            c = all95 - x
            d = allnon95-(totcnt-x)
            o_r, pval = fisher_exact([[a, b], [c, d]])
            pval = min(pval * domcnt, 1)
            print >> table, "\t".join(map(str,[pfam, name, counts[pfam], totlen, x, '{:.3g}'.format(pval)]))
            print o_r, pval
            if pval < 0.05:
                sign+=1
            for i in dbins:
                for j in range(0,len(tccr)):
                    if tccr[j][0]==i:
                        tccr[j][1]=(dbins[i][0]/totlen)
            vals = [i[1] for i in tccr]; binar = [i[0]-2.5 for i in tccr] # start at 0 and position bar in the middle of the bin
            #binar.insert(0,0)
            print vals, [i + 2.5 for i in binar]
            width = (binar[1]-binar[0])+.1
            maxval = max(maxval,max(vals))
            f, axarr = plt.subplots(1)#, sharex=True)
            axarr.bar(binar, vals, width = width, color = (161/255.0,218/255.0,215/255.0), alpha = 0.7, edgecolor = (96/255.0, 133/255.0, 131/255.0))
            axarr.xaxis.set_ticks(np.arange(0, 110, 5)) #0-100
            axarr.set_ylim(0,1)
            axarr.axhline(0.05, color = 'k', lw = 0.1, ls = '--', dashes = (60, 30))
            #axarr.plot(binar, vals)
            title += "\nFisher's OR: " + '{:.3g}'.format(o_r) + "; Bonferroni p-val: " + '{:.3g}'.format(pval)
            axarr.set_title(title)
            axarr.set_xlim(0,100) #10,100
            axarr.set_xlabel('CCR Percentile')
            axarr.set_ylabel('Frequency')
            plt.tight_layout()
            plt.subplots_adjust(hspace=0.7)
            ct+=1
    stats.seek(0)
#[make_axes_visible(i) for i in axarr]
#[i.set_xlabel('CCR Percentile') for i in axarr]
#[i.set_ylabel('Frequency') for i in axarr]
#plt.savefig(args.output, bbox_inches='tight')
pdf = matplotlib.backends.backend_pdf.PdfPages(args.output)
for fig in xrange(1, plt.gcf().number+1): ## will open an empty extra figure :(
    pdf.savefig(fig, bbox_inches='tight' )
pdf.close()
print maxval
table.close()
print sign
