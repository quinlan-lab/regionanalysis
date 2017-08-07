import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import cPickle as pickle
pk = open('gerp.pkl', 'rb')
data = pickle.load(pk)
pk.close()
cscores, gscores, labels = data

plt.scatter(cscores,gscores,marker='o',s=1)
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
plt.savefig('/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/gerpvccr.pdf', bbox_inches='tight')
