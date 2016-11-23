import sys
import numpy as np
f=open(sys.argv[1],"r")
d={}
for l in f:
    l=l.strip()
    if "random" in l:
        key=l
        if key not in d:
            d[key]=[]
    else:
        d[key].append(int(l))
for i in d:
    print i
    print int(np.mean(d[i]))
