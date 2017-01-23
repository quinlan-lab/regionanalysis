import sys
import subprocess

f3=open('genefractions.txt','w')
f3.write('study name	fraction of study covered by top 1%	fraction covered by random genes (top 1% length)	fraction of study covered by top 5%	fraction covered by random genes(top 5% length)	fraction covered by pLI > 0.9	fraction covered by random genes (pLI > 0.9 length)'+"\n")

tops=(1,5); d={}; totals={}
command='bash genecompare.sh '
command2='bash randomcompare.sh '
for i in tops:
    subprocess.check_call(command+str(i), shell=True)
    subprocess.check_call(command2+'/tmp/resids', shell=True)
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

command3='bash plicompare.sh'
subprocess.check_call(command3, shell=True)
subprocess.check_call(command2+'/tmp/pli', shell=True)
f4=open('pli.txt','r')
for line in f4:
    strings=line.strip().split('\t')
    if "pLI" in strings:
        key = strings[0]
        line=next(f4)
        d[key].append(line.strip().split()[0]+"/"+totals[key])
f4.close()
f2=open('randomintersects.txt','r')
for line in f2:
    strings=line.strip().split('\t')
    if "randoms" in strings:
        key = strings[0]
        line=next(f2)
        d[key].append(line.strip().split()[0]+"/"+totals[key])
f2.close()
for i in d.values():
    f3.write('\t'.join(i)+'\n')
