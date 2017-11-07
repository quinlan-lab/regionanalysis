import sys

prevkey, blockSizes, blockStarts = None, [], []
for line in sys.stdin: # exacresiduals/flatexome.bed sorted by gene
    fields = line.strip().split("\t")
    chrom = fields[0]; origstart = fields[1]; origend = fields[2]; gene = fields[3]
    d = {'chrom': chrom, 'start': origstart, 'end': origend, 'gene': gene}
    key = d['gene']
    if prevkey is None:
        start = d['start']
    elif key != prevkey:
        end = prevdict['end']
        print "\t".join([prevdict['chrom'], start, end, prevdict['gene']])
        start = d['start']
    prevdict = d; prevkey = key

end = prevdict['end']
print "\t".join([prevdict['chrom'], start, end, prevdict['gene']])
