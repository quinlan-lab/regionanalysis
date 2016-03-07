#miRNA (TS/TargetScan file, mirna.txt):
the higher the context score percentile, the lower the content score and the greater the fold repression of miRNA at the site
####try 80%+
#chromatin data (ChromHMM, chromatinseg.txt):
active promoter means active intergenic hypomethylated region
####ignore heterochromatin, CNVs, use 1,2,3,4,5,6,7,8
#dna methylation (took RRBS, not 450k for quantity reasons, methylation.txt)
any site should be okay
####try 50%+
#DNAse I hypersens. sites (dnase.txt)
the higher the score, the more sensitive.
####try greater than 400, doesn't seem to be many greater than 700
#TFBS sites (tfbs.txt)
the higher the score, the more signal seen.
####max is 1000, so try 700+
#histone modifications (histonemod.txt)
just like TFBS, higher score = more signal
####max is 1000, try 700+
#RNA binding sites (rnabinding.txt)
just like DNAse, the higher the score the more bindin
####doesn't seem to pass 700, so 400+
