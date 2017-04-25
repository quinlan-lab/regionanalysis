import sys
import numpy as np
from sklearn import metrics
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-t", "--title", help="Title")
parser.add_argument("-c", "--ccr", help="ccr pathogenic percentile file, then benign percentile file", nargs="*")
parser.set_defaults(ccr = ['tmp/ccrpatho','tmp/ccrbenign'])
parser.add_argument("-p", "--pli", help="pli pathogenic score file, then benign pli score file", nargs="*")
parser.set_defaults(ccr = ['tmp/plipatho','tmp/plibenign'])
parser.add_argument("-d", "--cadd", help="cadd pathogenic score file, then benign cadd score file", nargs="*")
parser.set_defaults(ccr = ['tmp/caddpatho','tmp/caddbenign'])
parser.add_argument("-g", "--gnomad", help="ccr pathogenic percentile file, then benign percentile file", nargs="*")
parser.set_defaults(ccr = ['tmp/ccr2patho','tmp/ccr2benign'])
parser.add_argument("-r", "--rvis", help="pathogenic rvis score file, then benign rvis score file", nargs="*")
parser.set_defaults(ccr = ['tmp/rvispatho','tmp/rvisbenign'])
args=parser.parse_args()
ccr=args.ccr
pli=args.pli
cadd=args.cadd
gnomad=args.gnomad
rvis=args.rvis

title=args.title
p=open(ccr[0],'r') #pathogenic percentiles from CCR
b=open(ccr[1],'r') #benign percentiles from CCR
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
plt.plot(fpr,tpr,label='CCR = ' + "%.3f" % (AUC), color = 'b')
p=open(gnomad[0],'r') #pathogenic percentiles from CCR
b=open(gnomad[1],'r') #benign percentiles from CCR
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
plt.plot(fpr,tpr,label='CCR w/gnomAD = ' + "%.3f " % (AUC), color = 'm')
p=open(pli[0],'r') #pathogenic pLI
b=open(pli[1],'r') #benign pLI
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
plt.plot(fpr,tpr,label='pLI = ' + "%.3f " % (AUC), color = 'g')
p=open(cadd[0],'r') #pathogenic CADD
b=open(cadd[1],'r') #benign CADD
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
plt.plot(fpr,tpr,label='CADD = ' + "%.3f " % (AUC), color = 'k')
p=open(rvis[0],'r') #pathogenic rvis
b=open(rvis[1],'r') #benign rvis
y=[]; scores=[]
for line in p:
    y.append(1)
    scores.append(100-float(line))
for line in b:
    y.append(0)
    scores.append(100-float(line))
y=np.array(y); scores=np.array(scores)
fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
AUC=metrics.roc_auc_score(y, scores)
plt.plot(fpr,tpr,label='RVIS = ' + "%.3f " % (AUC), color = 'r')
plt.title(title)
plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
plt.legend(loc='best')
plt.savefig('roc.pdf',bbox_inches='tight')
