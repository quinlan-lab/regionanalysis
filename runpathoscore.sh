HOME=/uufs/chpc.utah.edu/common/home/u1021864
DATA=/scratch/ucgd/lustre/u1021864/serial
# make truth sets
if [ ! -s $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-benign.20170802.vcf.gz ] | [ ! -s $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170802.vcf.gz ]; then
    cd $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/
    bash make.sh
    cd -
fi
if [ ! -s $HOME/software/pathoscore/truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz ] | [ ! -s $HOME/software/pathoscore/truth-sets/GRCh37/samocha/samocha.benign.vcf.gz ]; then
    cd $HOME/software/pathoscore/truth-sets/GRCh37/samocha/
    bash make.sh
    cd -
fi
# make score sets
if [ ! -s $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz ]; then
    cd $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/
    bash make.sh
    cd -
fi
if [ ! -s $DATA/whole_genome_SNVs.tsv.gz ] && [ ! -s $HOME/software/pathoscore/score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz ]; then
    cd $HOME/software/pathoscore/score-sets/GRCh37/CADD/
    bash make.sh
    cd -
fi
if [ ! -s $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz ]; then
    cd $HOME/software/pathoscore/score-sets/GRCh37/GERP/
    bash make.sh
    cd -
fi
if [ ! -s $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz ]; then
    cd $HOME/software/pathoscore/score-sets/GRCh37/MPC/
    bash make.sh
    cd -
fi
###########################################
###### CLINVAR ############################
###########################################
# run annotate to add scores to truth set files
python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz:CCR:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MCAP/mcap.txt.gz:mcap:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores $HOME/software/pathoscore/score-sets/GRCh37/pLI/pLI.bed.gz:pLI:4:max $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-benign.20170802.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --exclude $DATA/gnomad-vep-vt.vcf.gz --prefix benign

python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz:CCR:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/CADD/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MCAP/mcap.txt.gz:mcap:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores $HOME/software/pathoscore/score-sets/GRCh37/pLI/pLI.bed.gz:pLI:4:max $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170802.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --exclude $DATA/gnomad-vep-vt.vcf.gz --prefix pathogenic # /uufs/chpc.utah.edu/common/home/u1021864/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170802.vcf.gz

# add combined scores
python combine.py pathogenic.vcf.gz | bgzip -c -@ 12 > pathogenic.combine.vcf.gz
python combine.py benign.vcf.gz | bgzip -c -@ 12 > benign.combine.vcf.gz

# run evaluate to generate plot data
python $HOME/software/pathoscore/pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscore/clinvar -s CCR -s CADD -s GERP -s MPC -s REVEL -i MTR -s pLI --suffix pdf pathogenic.combine.vcf.gz benign.combine.vcf.gz #-s gnomad5x1_ccr -s gnomad5x9_ccr -s gnomad50x1_ccr -s gnomad50x9_ccr -s gnomad30x5_ccr

# get intersection of pathogenic benign and ccrs for odds ratio
bedtools intersect -a pathogenic.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > patho-ccr.txt
bedtools intersect -a benign.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > benign-ccr.txt

# generate fig 2 plot
python fig2plot.py clinvar

# ad gene files
python $HOME/software/pathoscore/pathoscore.py annotate pathogenic.combine.vcf.gz --exclude $HOME/software/pathoscore/gene-sets/GRCh37/ad_genes/ad_gene_complement.bed.gz --prefix adgene.pathogenic.combine
python $HOME/software/pathoscore/pathoscore.py annotate benign.combine.vcf.gz --exclude $HOME/software/pathoscore/gene-sets/GRCh37/ad_genes/ad_gene_complement.bed.gz --prefix adgene.benign.combine

# get intersection of ad genes pathogenic benign and ccrs for odds ratio
bedtools intersect -a adgene.pathogenic.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > patho-adgene-ccr.txt
bedtools intersect -a adgene.benign.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > benign-adgene-ccr.txt

# create oddsratio plot
python oddsratio.py -f patho-ccr.txt benign-ccr.txt -a patho-adgene-ccr.txt benign-adgene-ccr.txt -o clinvar

# hi gene files
python $HOME/software/pathoscore/pathoscore.py annotate pathogenic.combine.vcf.gz --exclude $HOME/software/pathoscore/gene-sets/GRCh37/hi_genes/hi_gene_complement.bed.gz --prefix higene.pathogenic.combine
python $HOME/software/pathoscore/pathoscore.py annotate benign.combine.vcf.gz --exclude $HOME/software/pathoscore/gene-sets/GRCh37/hi_genes/hi_gene_complement.bed.gz --prefix higene.benign.combine

# get intersection of hi genes pathogenic benign and ccrs for odds ratio
bedtools intersect -a higene.pathogenic.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > patho-higene-ccr.txt
bedtools intersect -a higene.benign.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > benign-higene-ccr.txt

# create oddsratio plot
python oddsratio.py -f patho-ccr.txt benign-ccr.txt -a patho-higene-ccr.txt benign-higene-ccr.txt -o clinvar.hi

# hi gene files
python $HOME/software/pathoscore/pathoscore.py annotate pathogenic.combine.vcf.gz --exclude $HOME/software/pathoscore/gene-sets/GRCh37/hi_genes/clingen_hi_gene_complement.bed.gz --prefix clingen_higene.pathogenic.combine
python $HOME/software/pathoscore/pathoscore.py annotate benign.combine.vcf.gz --exclude $HOME/software/pathoscore/gene-sets/GRCh37/hi_genes/clingen_hi_gene_complement.bed.gz --prefix clingen_higene.benign.combine

# get intersection of hi genes pathogenic benign and ccrs for odds ratio
bedtools intersect -a clingen_higene.pathogenic.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > patho-clingen_higene-ccr.txt
bedtools intersect -a clingen_higene.benign.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > benign-clingen_higene-ccr.txt

# create oddsratio plot
python oddsratio.py -f patho-ccr.txt benign-ccr.txt -a patho-clingen_higene-ccr.txt benign-clingen_higene-ccr.txt -o clinvar.clingen_hi

# table of top bin (>95%) CCRs by number of intersections with pathogenic variants:
# this line grabs only non-gnomAD or ExAC pathogenics that are functional
python filtervars.py pathogenic.combine.vcf.gz > clinvarfunc.vcf
bedtools intersect -a <(zcat exacresiduals/gnomad10x.5syn-ccrs.bed.gz | awk '$14>=95') -b clinvarfunc.vcf -wao > all-ccr-95-patho.txt
bedtools intersect -a <(zcat exacresiduals/gnomad10x.5syn-ccrs.bed.gz | awk '$14>=99') -b clinvarfunc.vcf -wao > all-ccr-99-patho.txt
python pathotable.py all-ccr-95-patho.txt "95"
python pathotable.py all-ccr-99-patho.txt "99"
grep -v "^#" clinvarfunc.vcf | wc -l 

###########################################
###### SAMOCHA ############################
###########################################

# run annotate to add scores to truth set files
python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz:CCR:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $DATA/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $DATA/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MCAP/mcap.txt.gz:mcap:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores $HOME/software/pathoscore/score-sets/GRCh37/pLI/pLI.bed.gz:pLI:4:max $HOME/software/pathoscore/truth-sets/GRCh37/samocha/samocha.benign.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --exclude $DATA/gnomad-vep-vt.vcf.gz --prefix control

python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5syn-ccrs.bed.gz:CCR:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $DATA/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $DATA/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MCAP/mcap.txt.gz:mcap:5:max --scores $HOME/software/pathoscore/score-sets/GRCh37/REVEL/revel.txt.gz:REVEL:7:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MTR/mtrflatfile_1.0.txt.gz:MTR:11:max --scores $HOME/software/pathoscore/score-sets/GRCh37/pLI/pLI.bed.gz:pLI:4:max $HOME/software/pathoscore/truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --exclude $DATA/gnomad-vep-vt.vcf.gz --prefix neurodev

# add combined scores
python combine.py neurodev.vcf.gz | bgzip -c -@ 12 > neurodev.combine.vcf.gz
python combine.py control.vcf.gz | bgzip -c -@ 12 > control.combine.vcf.gz

# run evaluate to generate plot data
python $HOME/software/pathoscore/pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscore/samocha -s CCR -s CADD -s GERP -s MPC -s REVEL -i MTR -s pLI --suffix pdf neurodev.combine.vcf.gz control.combine.vcf.gz #-s gnomad5x1_ccr -s gnomad5x9_ccr -s gnomad50x1_ccr -s gnomad50x9_ccr -s gnomad30x5_ccr

# get intersection of pathogenic and ccrs for odds ratio
bedtools intersect -a neurodev.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > neurodev-ccr.txt
bedtools intersect -a control.combine.vcf.gz -b exacresiduals/gnomad10x.5syn-ccrs.bed.gz -wb > control-ccr.txt

# create oddsratio plot
python oddsratio.py -f neurodev-ccr.txt control-ccr.txt -o samocha
# generate fig 2 plot

python fig2plot.py samocha 
