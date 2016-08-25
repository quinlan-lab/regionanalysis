#To run analysis here use branch master:

Folder located at: /uufs/chpc.utah.edu/common/home/u1021864/analysis ("analysis" or as it is titled here regionanalysis)

To change how residuals are defined, and scripts for creating the top and middle regions for Monte Carlo comparisons:

Use https://github.com/quinlan-lab/exacresiduals/tree/test (branch test is the most up-to-date)

Folder located at: /uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals ("exacresiduals")

##Under the folder "analysis":

There is a folder titled "ogfiles" (original files) where I put the files as they originally were either converted in Excel to TXT form if that was the case, or gene tsvs from MacArthur's lab or TXTs from UCSC containing TFBS regions or the like.  You probably don't need to use this folder but there are scripts that create the files used in the analysis that run using the files in this folder.

There is a folder titled "genescreens" that contains "wang.bed" which is the bed file created from taking out the genes with an adjusted p-val of >=0.5 from ogfiles/CSbysgRNA.txt.  That file was converted from Excel format and was used to match sgRNAs to their corresponding genes' CS (CRISPR Score, the lower, the more "essential").  It contains 4 lines, chrom start end CS_score.  The script containing the code used to create these files from the folder "ogfiles" is in the analysis folder and is called "essentialgenes.sh".  I would not simply run the script as is, it could take a while -- look at the code first.

There is a folder titled "denovos" that contains the de novos used for the DNM (de novo mutation) analysis.  The script again pulls from the folder "ogfiles" and is called "generatedenovobeds.sh".  The resultant DNM bed files are stored in the folder "denovos".

The folder "intermeds" (intermediate files) is where I store intermediate files used for creating the end-user files so I can check the different steps of the analysis and data conversion if necessary.  I would not use files in here unless you know what you're looking for.

The folder "ogcopy" is an copy of the original "ogfiles" folder when I was trying to stratify by sites like TFBS and chromatin binding sites, it may be deleted in the future can be ignored.

The folder "regions" is where the script "regions.sh" in the folder "exacresiduals" stores the "topresid.txt" and "midresid.txt" files for Monte Carlo analysis.  At one point we were trying to create a gene based mean score for residuals and those are stored in there for now as well.

The folder "plots" is pretty self-explanatory, it's where the plots from "plotcrispr.py", "plotdistro.py", "plotresid.py" and "plot.py" are stored.

The folder "results" is where the Monte Carlo overlap files in JSON format are stored for use by "plot.py".



##old comparison data:
###miRNA (TS/TargetScan file, mirna.txt):
the higher the context score percentile, the lower the content score and the greater the fold repression of miRNA at the site
####try 85%+ score (top 35%)
###chromatin data (ChromHMM, chromatinseg.txt):
active promoter means active intergenic hypomethylated region
####ignore heterochromatin, CNVs, use 1,2,3,4,5,6,7,8 (top 58%)
###dna methylation (took RRBS, not 450k for quantity reasons, methylation.txt)
any site should be okay
####try percentMeth=60%+ (top 25%)
###DNAse I hypersens. sites (dnase.txt)
the higher the score, the more sensitive.
####try greater than 400 signal score, doesn't seem to be many greater than 900 (top 21%)
###TFBS sites (tfbs.txt)
the higher the score, the more signal seen.
####max is 1000, so try 600+ (top 18%)
###histone modifications (histonemod.txt)
just like TFBS and DNAse I, higher score = more signal
####max is 1000, try 700+ (top 33%)
###RNA binding sites (rnabinding.txt)
just like DNAse, the higher the score the more binding
####doesn't seem to pass 700, so 300+ (top 20%)
