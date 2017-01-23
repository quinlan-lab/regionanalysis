import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

f = open(sys.argv[1],'r')
f.readline()
d={}
def transfrac(fraction):
    return eval('/'.join(map(str,map(float,fraction.split('/')))))

for line in f:
    fields=line.strip().split("\t")
    frac1=transfrac(fields[1]); frac2=transfrac(fields[2]); frac3=transfrac(fields[3])
    frac4=transfrac(fields[4]); frac5=transfrac(fields[5]); frac6=transfrac(fields[6])
    d[fields[0]]=frac1,frac2,frac3,frac4,frac5,frac6

N = len(d)
ind = np.arange(N)
width=0.30

sns.set_style('ticks')
fig, ax = plt.subplots()
sns.despine(right=True,top=True)
rects1 = ax.barh(ind,[i[0] for i in d.values()], width-.05, color = 'w', lw=2, alpha=0.7) # hatch='x')
rects2 = ax.barh(ind,[i[1] for i in d.values()], width-.05, color = 'k', lw=2, alpha=0.8)
rects3 = ax.barh(ind+width,[i[2] for i in d.values()], width-.05, color = 'dodgerblue', lw=2, alpha=0.7) # hatch='|')
rects4 = ax.barh(ind+width,[i[3] for i in d.values()], width-.05, color = 'k', lw=2, alpha=0.8)
rects5 = ax.barh(ind+width*2,[i[4] for i in d.values()], width-.05, color = 'darkorange', lw=2, alpha=0.7) # hatch='//')
rects6 = ax.barh(ind+width*2,[i[5] for i in d.values()], width-.05, color = 'k', lw=2, alpha=0.8)
ax.set_yticks(ind+width+.05)
ax.set_yticklabels([i for i in d.keys()])
ax.set_xlabel('Fraction of Genes in top CCRs')
ax.set_xlim(0,1)
plt.title('Enrichment for Functionally Important Genes in CCRs')

ax.legend((rects3[0], rects1[0], rects5[0], rects2[0]), ('top 5% CCRs', 'top 1% CCRs', 'pLI > 0.9', 'random'),loc='upper right')

def autolabel(rects1, rects2):
    # attach some text labels
    for rect1,rect2, in zip(rects1,rects2):
        width = rect1.get_width()
        val = round(width/float(rect2.get_width()),1)
        ax.text(1.05*width,rect1.get_y() + rect1.get_height()/2.,
                str(val)+"x", weight='bold')

autolabel(rects1,rects2); autolabel(rects3,rects4); autolabel(rects5,rects6)
plt.savefig('bars.png', bbox_inches='tight')
