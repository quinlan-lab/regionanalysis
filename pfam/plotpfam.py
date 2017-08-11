import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import numpy as np
import ast
import seaborn
seaborn.set_style(style='white')

parser=ArgumentParser()
parser.add_argument("-p","--pfam", help="pfam families for histogram") # top100doms
parser.add_argument("-c","--clans", help="pfam clans file, has proper family names in it") # Pfam-A.clans.tsv.gz
parser.add_argument("-s","--stats", help="pfams bp covered in ccr %ile bins") #pfamshist.txt
parser.add_argument("-o","--output", help="histogram pic") #pfam_hists.pdf
args=parser.parse_args()
with open(args.pfam) as f:
    pfams=[_ for _ in f]
families={}
with open(args.clans) as f:
    for line in f:
        fields=line.strip().split("\t")
        families[fields[3]]=fields[4] # pfam id is key, true family name is value

count=len(pfams)
stats=open(args.stats,"r")
bins=range(10,101,10)

def make_axes_visible(axes):
    for ax in axes.flatten():
        for tk in ax.get_yticklabels():
            tk.set_visible(True)
        for tk in ax.get_xticklabels():
            tk.set_visible(True)

ct=0
plt.rcParams["figure.figsize"]=(5,20*count/10)
#f, axarr = plt.subplots(count, sharex=True)
axarr = [plt.subplot(count, 1, i+1) for i in range(count)]
for pfam in pfams:
    pfam=pfam.strip()
    print pfam
    for stat in stats:
        fields=stat.strip().split("\t")
        pfam2=fields[0].strip()
        dbins=ast.literal_eval(fields[1]) #(dictionary of bins)
        if pfam == pfam2:
            try:
                title=families[pfam]
            except KeyError:
                title=pfam
            tccr=[[i,0] for i in bins]
            totlen=float(fields[2])
            for i in dbins:
                for j in range(0,len(tccr)):
                    if tccr[j][0]==i:
                        tccr[j][1]=(dbins[i]/totlen)
            vals=[i[1] for i in tccr]; binar=[i[0]-5 for i in tccr] # start at 0 and position bar in the middle of the bin
            #binar.insert(0,0)
            print vals, [i + 5 for i in binar]
            width = (binar[1]-binar[0])+.1
            axarr[ct].bar(binar, vals, width = width, color = 'b', alpha = 0.7)
            axarr[ct].xaxis.set_ticks(np.arange(0, 110, 10)) #0-100
            axarr[ct].axhline(0.1, color = 'k')
            #axarr[ct].plot(binar, vals)
            axarr[ct].set_title(title)
            axarr[ct].set_xlim(0,100) #10,100
            ct+=1
    stats.seek(0)
#[make_axes_visible(i) for i in axarr]
plt.tight_layout()
plt.subplots_adjust(hspace=0.7)
[i.set_xlabel('CCR Percentile') for i in axarr]
[i.set_ylabel('Frequency') for i in axarr]
plt.savefig(args.output, bbox_inches='tight')
