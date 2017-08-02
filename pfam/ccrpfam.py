import sys

nodoms, doms = open("tmp/nodomccrs", "w"), open("tmp/domccrs", "w")
for line in sys.stdin:
    fields = line.strip().split("\t")
    ccr = fields[:14]
    pfam = fields[14:]
    #if float(ccr[-1]) == 0:
    #    continue
    if pfam[0] == ".":
        nodoms.write(ccr[-1]+"\n")
    else:
        doms.write(ccr[-1]+"\n")
