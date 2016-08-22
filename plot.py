import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import json
import sys

folder=sys.argv[2]

figname = sys.argv[1].split("/")[-1].split(".")[0]
filename = "plots/"+folder+"/"+figname+".png"
f = open(sys.argv[1],'r')
j = f.readline()
j = json.loads(j)["'wc -l'"]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(j['sims'],alpha=0.5) # if you get an error here, uncomment the commented out met['sims'] = sims part of poverlap.py
plt.xlabel("Intersections")
plt.ylabel("Frequency")
plt.title(figname)
plt.axvline(j['observed'], color='r', linestyle='dashed', linewidth=2)
xmin = ax.get_xlim()[0]
ymin = ax.get_ylim()[0]
xmax = ax.get_xlim()[-1]
ymax = ax.get_ylim()[-1]
plt.text(j['observed']+.1,.5*ymax,'observed',rotation=90)
plt.axvline(j['simulated mean metric'], color='b', linestyle='dashed', linewidth=2)
plt.text(j['simulated mean metric']+.1,.5*ymax,'simulated mean',rotation=90)
plt.text(.04*(xmax-xmin)+xmin,.85*(ymax-ymin)+ymin,"pval:\n"+"%.3e" % j['simulated_p'],bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))
plt.savefig("plots/"+folder+"/"+figname+".png", bbox_inches = 'tight')
plt.close()
