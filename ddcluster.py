from __future__ import print_function
import sys
from itertools import groupby
from operator import itemgetter
import pyexcel as pe

book = pe.get_book(file_name=sys.argv[1])

fh = open('ddvarclusters.bed', 'w')
header = "chrom\tstart\tend\tgene"
#header = "chrom\tstart\tend\tgene\tref\talt\timpact\tgene_id\thgvs_gdna\thgvs_cdna\thgvs_protein\tsample_id\tpubmed_id"
print(header, file=fh)

sheet = book['Sheet1']
tempdicts=[]
for i, record in enumerate(sheet):
    if i == 0:
        keys = map(str, record)
        continue
    record = dict(zip(keys, record))
    tempdicts.append(record)
grouper = itemgetter("Gene name")
for key, grp in groupby(sorted(tempdicts, key = grouper), grouper):
    grp=list(grp)
    chrom = grp[0]["Chromosome"].split('chr')[-1]
    start = min(item["Start position"] for item in grp)
    end = max(item["End position"] for item in grp)

    print ("\t".join(map(str,[chrom,start,end,key])), file=fh)
    
   # print("{chrom}\t{start}\t{end}\t{gene}\t{ref}\t{alt}\t{impact}\t{gene_id}\t{hgvs_gdna}\t{hgvs_cdna}\t{hgvs_protein}\t{sample_id}\t{pubmed_id}".format(
   #     chrom=record['Chromosome'].split('chr')[-1], start=record['Start position'], end=record['End position'], gene=record['Gene name'],
   #     ref=record['Reference'], alt=record['Variant'], impact=record['Protein Effect'], gene_id=record['Gene id'],
   #     hgvs_gdna=record['HGVS gDNA'], hgvs_cdna=record['HGVS cDNA'], hgvs_protein=record['HGVS protein'],
   #     sample_id=record['SampleID'], pubmed_id=record['PubmedID']), file=fh)

fh.close()
