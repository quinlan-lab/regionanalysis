import sys

f = open('regions/topresid.txt','r')

totbp = 0
for line in f:
    fields = line.strip().split("\t")
    totbp += int(fields[2]) - int(fields[1])
f.close()
f = open('residualgist/residuals.txt', 'r')
list_ = []
f.readline()
for line in f:
    fields = line.strip().split("\t")
    list_.append(fields)
bp = 0 
list_ = sorted(list_,key = lambda x: abs(float(x[10])))
for line in list_:
    if bp > totbp:
        break
    bp += int(line[2]) - int(line[1])
    print "\t".join(line)
