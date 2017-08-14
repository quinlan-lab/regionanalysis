HOME=/uufs/chpc.utah.edu/common/home/u1021864
DATA=/scratch/ucgd/lustre/u1021864/serial
# make truth sets
if [ ! -s $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-benign.20170710.vcf.gz ] | [ ! -s $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz ]; then
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
if [ ! -s $DATA/whole_genome_SNVs.tsv.gz ] & [ ! -s $HOME/software/pathoscore/score-sets/GRCh37/CADD/whole_genome_SNVs.tsv.gz ]; then
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
python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz:gnomad10x5_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.1-ccrs.bed.gz:gnomad5x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.9-ccrs.bed.gz:gnomad5x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.1-ccrs.bed.gz:gnomad50x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.9-ccrs.bed.gz:gnomad50x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad30x.5-ccrs.bed.gz:gnomad30x5_ccr:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $DATA/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $DATA/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-benign.20170710.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --prefix benign

python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz:gnomad10x5_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.1-ccrs.bed.gz:gnomad5x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.9-ccrs.bed.gz:gnomad5x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.1-ccrs.bed.gz:gnomad50x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.9-ccrs.bed.gz:gnomad50x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad30x.5-ccrs.bed.gz:gnomad30x5_ccr:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $DATA/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $DATA/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max $HOME/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --prefix pathogenic # /uufs/chpc.utah.edu/common/home/u1021864/software/pathoscore/truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz

# add combined scores
python combine.py pathogenic.vcf.gz | bgzip -c -@ 12 > pathogenic.combine.vcf.gz
python combine.py benign.vcf.gz | bgzip -c -@ 12 > benign.combine.vcf.gz

# run evaluate to generate plot data
python $HOME/software/pathoscore/pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscore/clinvar -s gnomad10x5_ccr -s pp2hdiv -s CADD -s GERP -s MPC -s mpc_regions -s Grantham -s mis_badness -s ccr+cadd -s ccr+polyphen -s ccr+gerp -s ccr+grantham -s ccr+misbadness pathogenic.combine.vcf.gz benign.combine.vcf.gz #-s gnomad5x1_ccr -s gnomad5x9_ccr -s gnomad50x1_ccr -s gnomad50x9_ccr -s gnomad30x5_ccr

# generate fig 2 plot
python fig2plot.py clinvar

###########################################
###### SAMOCHA ############################
###########################################

# run annotate to add scores to truth set files
python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz:gnomad10x5_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.1-ccrs.bed.gz:gnomad5x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.9-ccrs.bed.gz:gnomad5x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.1-ccrs.bed.gz:gnomad50x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.9-ccrs.bed.gz:gnomad50x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad30x.5-ccrs.bed.gz:gnomad30x5_ccr:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $DATA/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $DATA/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max $HOME/software/pathoscore/truth-sets/GRCh37/samocha/samocha.benign.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --prefix control

python $HOME/software/pathoscore/pathoscore.py annotate --scores $HOME/analysis/exacresiduals/gnomad10x.5-ccrs.bed.gz:gnomad10x5_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.1-ccrs.bed.gz:gnomad5x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad5x.9-ccrs.bed.gz:gnomad5x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.1-ccrs.bed.gz:gnomad50x1_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad50x.9-ccrs.bed.gz:gnomad50x9_ccr:14:max --scores $HOME/analysis/exacresiduals/gnomad30x.5-ccrs.bed.gz:gnomad30x5_ccr:14:max --scores $HOME/software/pathoscore/score-sets/GRCh37/polyphen2/polyphen2/polyphen2.txt.gz:pp2hdiv:5:max --scores $DATA/whole_genome_SNVs.tsv.gz:CADD:6:max --scores $DATA/InDels.tsv.gz:CADD:6:max --scores $HOME/software/pathoscore/score-sets/GRCh37/GERP/gerp_rs.txt.gz:GERP:3:max --scores $HOME/software/pathoscore/score-sets/GRCh37/MPC/mpc.txt.gz:MPC:5:max --scores $HOME/analysis/essentials/mpc.regions.clean.sorted.bed.gz:mpc_regions:5:max $HOME/software/pathoscore/truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz --exclude $DATA/ExAC.r1.vt.vep.vcf.gz --prefix neurodev

# add combined scores
python combine.py neurodev.vcf.gz | bgzip -c -@ 12 > neurodev.combine.vcf.gz
python combine.py control.vcf.gz | bgzip -c -@ 12 > control.combine.vcf.gz

# run evaluate to generate plot data
python $HOME/software/pathoscore/pathoscore.py evaluate --functional --prefix $HOME/public_html/pathoscore/samocha -s gnomad10x5_ccr -s pp2hdiv -s CADD -s GERP -s MPC -s mpc_regions -s Grantham -s mis_badness -s ccr+cadd -s ccr+polyphen -s ccr+gerp -s ccr+grantham -s ccr+misbadness neurodev.combine.vcf.gz control.combine.vcf.gz #-s gnomad5x1_ccr -s gnomad5x9_ccr -s gnomad50x1_ccr -s gnomad50x9_ccr -s gnomad30x5_ccr

# generate fig 2 plot
python fig2plot.py samocha 
