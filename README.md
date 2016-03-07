#miRNA (TS/TargetScan file, mirna.txt):
the higher the context score percentile, the lower the content score and the greater the fold repression of miRNA at the site
####try 85%+ score (top 35%)
#chromatin data (ChromHMM, chromatinseg.txt):
active promoter means active intergenic hypomethylated region
####ignore heterochromatin, CNVs, use 1,2,3,4,5,6,7,8 (top 58%)
#dna methylation (took RRBS, not 450k for quantity reasons, methylation.txt)
any site should be okay
####try percentMeth=60%+ (top 25%)
#DNAse I hypersens. sites (dnase.txt)
the higher the score, the more sensitive.
####try greater than 400 signal score, doesn't seem to be many greater than 900 (top 21%)
#TFBS sites (tfbs.txt)
the higher the score, the more signal seen.
####max is 1000, so try 600+ (top 18%)
#histone modifications (histonemod.txt)
just like TFBS and DNAse I, higher score = more signal
####max is 1000, try 700+ (top 33%)
#RNA binding sites (rnabinding.txt)
just like DNAse, the higher the score the more binding
####doesn't seem to pass 700, so 300+ (top 20%)
