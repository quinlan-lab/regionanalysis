import itertools as it
import sys
import subprocess

def read_exons(gtf):
    transcripts = defaultdict(pyinter.IntervalSet)

    for toks in (x.rstrip('\r\n').split("\t") for x in open(gtf) if x[0] != "#"):
        if toks[2] not in("UTR", "exon"): continue
        #if toks[0] != "1": break
        start, end = map(int, toks[3:5])
        assert start <= end, toks
        transcript = toks[8].split('transcript_id "')[1].split('"', 1)[0]
        transcripts[transcript].add(pyinter.closedopen(start-1, end))

    # sort by start so we can do binary search.
    # TODO: need to remove overlapping exons so we don't double-count
    transcripts = dict((k, sorted(v)) for k, v in transcripts.iteritems())
    #ends = dict((k, sorted(v)) for k, v in ends.iteritems())
    starts, ends = {}, {}
    for tr, ivset in transcripts.iteritems():
        sends = sorted(list(ivset))
        starts[tr] = [x.lower_value for x in sends]
        ends[tr] = [x.upper_value for x in sends]
    return starts, ends

E = {}
U = {}

f = open(sys.argv[1],"r")
for l in f:
    e = Record(l.rstrip().split("\t"))
    if not e.transid in E:
        E[e.transid] = []
    E[e.transid].append(e)
f.close()

f = open(sys.argv[2],"r")
for l in f:
    u = Record(l.rstrip().split("\t"))
    if not u.transid in U:
        U[u.transid] = []
    U[u.transid].append(u)
f.close()

otal_dropped = 0

for e_transid in E:
    if not e_transid in U:
        x=1
        for e in E[e_transid]:
            print e
    else:
        f = open('u.tmp', 'w')
        f.write('\n'.join( \
                [str(u) for u in sorted(U[e_transid], \
                                        key=lambda x:(x.chr,x.start,x.end))]))
        f.close()

        f = open('e.tmp', 'w')
        f.write('\n'.join( \
                [str(e) for e in sorted(E[e_transid], \
                                        key=lambda x:(x.chr,x.start,x.end))]))

        f.close()

        f = open('r.tmp', 'w')
        p = subprocess.Popen('bedtools subtract -sorted -a e.tmp -b u.tmp',
                             shell=True,
                             stdout=f)
        p.wait()
        f.close()

        f = open('r.tmp', 'r')
        num_lines = sum(1 for line in f)
        f.close()

        total_dropped += len(E[e_transid]) - num_lines

        f = open('r.tmp', 'r')
        for l in f:
            print l.rstrip()
        f.close()

sys.stderr.write('A total of ' + str(total_dropped) + 'exons removed' + '\n')

