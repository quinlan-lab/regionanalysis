import sys
import numpy as np
from sklearn import metrics
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-t", "--title", help="Title")
parser.add_argument("-c", "--ccr", help="ccr based on exac; pathogenic percentile file, then benign percentile file", nargs="*")
#parser.set_defaults(ccr = ['tmp/ccrpatho','tmp/ccrbenign'])
parser.add_argument("-p", "--pli", help="pli pathogenic score file, then benign pli score file", nargs="*")
#parser.set_defaults(pli = ['tmp/plipatho','tmp/plibenign'])
parser.add_argument("-d", "--cadd", help="cadd pathogenic score file, then benign cadd score file", nargs="*")
#parser.set_defaults(cadd = ['tmp/caddpatho','tmp/caddbenign'])
parser.add_argument("-g", "--gnomad", help="ccr based on gnomad; pathogenic percentile file, then benign percentile file", nargs="*")
#parser.set_defaults(gnomad = ['tmp/ccr2patho','tmp/ccr2benign'])
parser.add_argument("-r", "--rvis", help="pathogenic rvis score file, then benign rvis score file", nargs="*")
#parser.set_defaults(rvis = ['tmp/rvispatho','tmp/rvisbenign'])
parser.add_argument("-b", "--combined", help="CADD+CCR pathogenic score file, then CADD+CCR benign score file", nargs="*")
#parser.set_defaults(combined = ['tmp/eccaddpatho','tmp/eccaddbenign'])
parser.add_argument("-m", "--mpc", help="MPC pathogenic score file, then MPC benign score file", nargs="*")
#parser.set_defaults(mpc = ['tmp/MPCpatho','tmp/MPCbenign'])
parser.add_argument("-o", "--outfile", help="output file, preferably something roc.pdf")
#parser.set_defaults(outfile = 'roc.pdf')
parser.add_argument("-n", "--numvars", help="number of pathogenic and then number of benign variants", nargs="*")
#parser.set_defaults(numvars = ['1000', '1000'])
args=parser.parse_args()
ccr=args.ccr
pli=args.pli
cadd=args.cadd
combined=args.combined
gnomad=args.gnomad
rvis=args.rvis
mpc=args.mpc
title=args.title
outfile=args.outfile
numvars=args.numvars
pathoct=numvars[0]
benignct=numvars[1]

variants=""

if ccr:
    ep=0; eb=0
    p=open(ccr[0],'r') #pathogenic percentiles from CCR
    b=open(ccr[1],'r') #benign percentiles from CCR
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(float(line))
        ep+=1
    for line in b:
        y.append(0)
        scores.append(float(line))
        eb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='CCR = ' + "%.3f" % (AUC), color = 'b')
    variants+="eCCRp: "+str(ep)+"/"+str(pathoct)+"; eCCRb: " +str(eb)+"/"+str(benignct)+"\n"
if gnomad:
    gp=0; gb=0
    p=open(gnomad[0],'r') #pathogenic percentiles from CCR
    b=open(gnomad[1],'r') #benign percentiles from CCR
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(float(line))
        gp+=1
    for line in b:
        y.append(0)
        scores.append(float(line))
        gb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='CCR w/gnomAD = ' + "%.3f " % (AUC), color = 'm')
    variants+="gCCRp: "+str(gp)+"/"+str(pathoct)+"; gCCRb: " +str(gb)+"/"+str(benignct)+"\n"
if pli:    
    pp=0; pb=0
    p=open(pli[0],'r') #pathogenic pLI
    b=open(pli[1],'r') #benign pLI
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(float(line))
        pp+=1
    for line in b:
        y.append(0)
        scores.append(float(line))
        pb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='pLI = ' + "%.3f " % (AUC), color = 'g')
    if len(pli)>2:
        pp=pli[2]; pb=pli[3]
    variants+="pLIp: "+str(pp)+"/"+str(pathoct)+"; pLIb: " +str(pb)+"/"+str(benignct)+"\n"
if cadd:
    cp=0; cb=0
    p=open(cadd[0],'r') #pathogenic CADD
    b=open(cadd[1],'r') #benign CADD
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(float(line))
        cp+=1
    for line in b:
        y.append(0)
        scores.append(float(line))
        cb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='CADD = ' + "%.3f " % (AUC), color = 'k')
    variants+="CADDp: "+str(cp)+"/"+str(pathoct)+"; CADDb: " +str(cb)+"/"+str(benignct)+"\n"
if rvis:
    rp=0; rb=0
    p=open(rvis[0],'r') #pathogenic rvis
    b=open(rvis[1],'r') #benign rvis
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(100-float(line))
        rp+=1
    for line in b:
        y.append(0)
        scores.append(100-float(line))
        rb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='RVIS = ' + "%.3f " % (AUC), color = 'r')
    if len(rvis)>2:
        rp=rvis[2]; rb=rvis[3]
    variants+="RVISp: "+str(rp)+"/"+str(pathoct)+"; RVISb: " +str(rb)+"/"+str(benignct)+"\n"
if combined:
    cp=0; cb=0
    p=open(combined[0],'r') #pathogenic CCR+CADD
    b=open(combined[1],'r') #benign CCR+CADD
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(float(line))
        cp+=1
    for line in b:
        y.append(0)
        scores.append(float(line))
        cb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='CADD+CCR = ' + "%.3f " % (AUC), color = 'c')
    variants+="CADD+CCRp: "+str(cp)+"/"+str(pathoct)+"; CADD+CCRb: " +str(cb)+"/"+str(benignct)+"\n"

if mpc:
    mp=0; mb=0
    p=open(mpc[0],'r') #pathogenic MPC
    b=open(mpc[1],'r') #benign MPC
    y=[]; scores=[]
    for line in p:
        y.append(1)
        scores.append(float(line))
        mp+=1
    for line in b:
        y.append(0)
        scores.append(float(line))
        mb+=1
    y=np.array(y); scores=np.array(scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    AUC=metrics.roc_auc_score(y, scores)
    plt.plot(fpr,tpr,label='MPC = ' + "%.3f " % (AUC), color = 'm')
    variants+="MPCp: "+str(mp)+"/"+str(pathoct)+"; MPCb: " +str(mb)+"/"+str(benignct)+"\n"

plt.title(title+"\n"+variants)
plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
plt.legend(loc='best')
plt.savefig(outfile,bbox_inches='tight')
