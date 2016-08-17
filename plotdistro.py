import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import subprocess as sp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", nargs="*", help="the list of your files, in order: constrained regions, pathogenic denovos, control denovos (if applicable)", required=True)
parser.add_argument("-s", "--study-name", help="name of your study, for plot label", required=True)
parser.add_argument("-r", "--range", help="range of second plot",nargs=2, default=[90,100])
parser.add_argument("-n", "--no-zoom", help="use this argument if you don't want a second plot", action="store_true", default=False)
args=parser.parse_args()
files=args.files
study_name=args.study_name
zoom=args.range

#TODO: argparse, add more options?

#plotting the distribution of residuals for pathogenic vs benign

lengths=[]
f = open(files[0], 'r')
f.readline()
for line in f:
    fields = line.strip().split("\t")
    lengths.append(int(fields[2])-int(fields[1]))

'''
intersecting function for beds, particularly, benign/pathogenic bed/vcf files
'''
def intersect(hotspot):
    def killproc(p):
        try:
            p.kill()
        except OSError:
            pass
    p1 = sp.Popen("sed '1d' "+files[0], shell = True, stdout = sp.PIPE)
    p2 = sp.Popen(['bedtools', 'intersect', '-a', 'stdin', '-b', hotspot], stdin = p1.stdout, stdout = sp.PIPE)
    output,error = p2.communicate()
    killproc(p1); killproc(p2)
    return output.strip()

bentile,pathtile = [],[]
if len(files)>2:
    benign = intersect(files[2]) #sys.argv[1]='benign.bed'
    for line in benign.split("\n"):
        fields = line.strip().split("\t")
        bentile.append(float(fields[11])) # percentile from residuals.txt
patho = intersect(files[1]) #sys.argv[2]='patho.bed'
for line in patho.split("\n"):
    fields = line.strip().split("\t")
    pathtile.append(float(fields[11]))

'''
distribution of residuals - benign vs pathogenic
'''
fig = plt.figure()
ax = fig.add_subplot(111)
mi = np.min(bentile+pathtile)
ma = np.max(bentile+pathtile)
rng = (mi,ma)
p,p_edges=np.histogram(pathtile, bins=40, range=rng)
p=map(lambda x: float(x)/sum(p), p)
width_p = (p_edges[1]-p_edges[0])
ax.bar(p_edges[:-1], p, width = width_p, color = 'r', label = 'pathogenic', alpha = 0.3)
if len(files)>2:
    b,b_edges=np.histogram(bentile, bins=40, range=rng)
    b=map(lambda x: float(x)/sum(b), b)
    width_b = (b_edges[1]-b_edges[0])
    ax.bar(b_edges[:-1], b, width = width_b, color = 'b', label = 'benign', alpha = 0.3)
ax.set_xlabel("Percentile")
ax.set_ylabel("Frequency")
plt.title("Percentile Pathogenicity Comparison for "+study_name)
ax.legend(loc='upper left')
plt.savefig("plots/"+study_name+"distribution.png")
plt.close()

if not args.no_zoom:

    '''
    top 10% distribution, benign vs pathogenic
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #bentile = [i for i in bentile if i >= 90]
    #pathtile = [i for i in pathtile if i >= 90]
    mi = np.min(bentile+pathtile)
    ma = np.max(bentile+pathtile)
    rng = (mi,ma)
    p,p_edges=np.histogram(pathtile, bins=400, range=rng)
    p=map(lambda x: float(x)/sum(p), p)
    width_p = (p_edges[1]-p_edges[0])
    ax.bar(p_edges[:-1], p, width = width_p, color = 'r', label = 'pathogenic', alpha = 0.3)
    if len(files)>2:
        b,b_edges=np.histogram(bentile, bins=400, range=rng)
        b=map(lambda x: float(x)/sum(b), b)
        width_b = (b_edges[1]-b_edges[0])
        ax.bar(b_edges[:-1], b, width = width_b, color = 'b', label = 'benign', alpha = 0.3)
    ax.set_xlim(float(zoom[0]), float(zoom[1]))
    ax.set_xlabel("Percentile")
    ax.set_ylabel("Frequency")
    plt.title("Percentile Pathogenicity Comparison: Top 10% "+study_name)
    ax.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig("plots/"+study_name+str(zoom[0])+"-"+str(zoom[1])+"distribution.png")
    plt.close()
