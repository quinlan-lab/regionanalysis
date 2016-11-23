import sys
import subprocess

f3=open('genefractions.txt','w')

print 'study name	fraction of study covered by top 1%	fraction covered by random genes (top 1% length)	fraction of study covered by top 5%	fraction covered by random genes(top 5% length)	fraction covered by top 10%	fraction covered by random genes (top 10% length)'
tops=(1,5,10); d={}; totals={}
command='bash genecompare.sh '
for i in tops:
    subprocess.check_call(command+str(i), shell=True)
    f1=open('intersects.txt','r')
    for line in f1:
        strings=line.strip().split('\t')
        if "total" in strings:
            key = strings[0]
            if key not in d.keys():
                d[key]=[key]
            line=next(f1)
            totals[key]=line.strip().split()[0]
        if "intersections" in strings:
            key = strings[0]
            line=next(f1)
            d[key].append(line.strip().split()[0]+"/"+totals[key])
    f1.close()
    f2=open('randomintersects.txt','r')
    for line in f2:
        strings=line.strip().split('\t')
        if "randoms" in strings:
            key = strings[0]
            line=next(f2)
            d[key].append(line.strip().split()[0]+"/"+totals[key])
    f2.close()

for i in d.values():
    print '\t'.join(i)
