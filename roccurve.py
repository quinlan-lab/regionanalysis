import sys
import numpy as np
from sklearn import metrics
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

title=sys.argv[1]
p=open('ryan/patho','r')
b=open('ryan/benign','r')
y=[]; scores=[]
for line in p:
    y.append(1)
    scores.append(float(line))
for line in b:
    y.append(0)
    scores.append(float(line))
y=np.array(y); scores=np.array(scores)
fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
AUC=metrics.roc_auc_score(y, scores)
plt.plot(fpr,tpr,label='AUC = '+str(AUC))
plt.title(title)
plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
plt.legend(loc='upper left')
plt.savefig('roc.png',bbox_inches='tight')
