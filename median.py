import numpy as np
import gzip

l95,l99=[],[]
with gzip.open('essentials/gnomadbased-ccrs.bed.gz', 'rb') as f:
    for line in f:
        fields = line.strip().split("\t")
        ranges=fields[6]; ccr=float(fields[-1])
        length=sum([float(i.split("-")[1]) - float(i.split("-")[0]) for i in ranges.split(",")])
        if ccr >= 95:
            l95.append(length)
        if ccr >= 99:
            l99.append(length)

print np.median(l95), np.median(l99)
