import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import toolshed as ts

lens=[]
it = ts.reader(sys.argv[1]) #$DATA/unfilteredregions.txt
iterable = (i for i in it)
for region in iterable:
    lens.append(int(region['end'])-int(region['start']))

print sum(lens)
totlen=0.0; ordered=[]; exomecov=[]; wgtdpct=[]
it = ts.reader(sys.argv[2]) #$HOME/analysis/exacresiduals/results/2016_12_10/weightedresiduals.txt (sorted by weighted_pct)
iterable = (i for i in it)
for region in iterable:
    length=int(region['end'])-int(region['start'])
    totlen+=length
    ordered.append((length,region['weighted_pct']))
ordered=sorted(ordered,key=lambda x:x[1])
length=0.0
it = ts.reader(sys.argv[2])
iterable = (i for i in it)
prevregion=None
for region in ordered:
    length+=(int(region['end'])-int(region['start']))
    if region['ranges']+region['gene']!=prevregion or prevregion!=None:
        exomecov.append(prevvalues[0])
        wgtdpct.append(prevvalues[1])
    prevregion=region['ranges']+region['gene']
    prevvalues=(length/totlen,region['weighted_pct'])
print exomecov, wgtdpct
