HOME=/uufs/chpc.utah.edu/common/home/u1021864
DATA=/scratch/ucgd/lustre/u1021864/serial
# table from: https://www.nature.com/gim/journal/vaop/ncurrent/fig_tab/gim201731t2.html
# fasta combined from: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
# idea from: https://www.biostars.org/p/15992/
if [ ! -s $DATA/hg19.fasta.dict ]; then
    java -jar $HOME/software/picard/build/libs/picard.jar CreateSequenceDictionary R=$DATA/hg19.fasta O=$DATA/hg19.fasta.dict
fi
if [ ! -s hg19.fa ]; then
    sed 's/^>chr/^>/g' $DATA/hg19.fasta > hg19.fa # needed for fetalcoords.py
fi
python fetalcoords.py ogfiles/fetalgenomevariants.txt ogfiles/fetaltranscripts.txt fetalvariants.vcf
cat $HOME/analysis/ogfiles/fetalproteinvariants.txt | java -jar $HOME/software/jvarkit/dist/backlocate.jar -R $DATA/hg19.fasta | sed 's/\tchr/\t/g' > /tmp/backlocate.txt # used to make fetalgenomic.txt
# fetalgenomic.txt is made by pulling out wrong amino acid matching variant backlocations or ones that could not possibly change amina acid 1 to amino acid 2 by hand; http://bioinfo.bisr.res.in/project/crat/pictures/codon.jpg
# translatecoords.py /tmp/backlocate.txt fetalgenomic.txt
cat fetalgenomic.txt >> fetalvariants.vcf
