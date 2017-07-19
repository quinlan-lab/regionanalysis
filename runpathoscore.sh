# make truth sets, copy score files
cd ~/software/pathoscore
cd truth-sets/GRCh37/clinvar/
bash make.sh
cd -
cp ~/analysis/gnomad-ccrs.bed.gz ~/analysis/gnomad-ccrs.bed.gz.tbi .
# run annotate to add scores to truth set files
python pathoscore.py annotate --scores gnomad-ccrs.bed.gz:gnomad_ccr:14:max truth-sets/GRCh37/clinvar/clinvar-benign.20170710.vcf.gz --prefix benign
python pathoscore.py annotate --scores gnomad-ccrs.bed.gz:gnomad_ccr:14:max truth-sets/GRCh37/clinvar/clinvar-pathogenic-likely_pathogenic.20170710.vcf.gz --prefix pathogenic
# run evaluate to generate plot data
python pathoscore.py evaluate -s exac_ccr pathogenic.vcf.gz benign.vcf.gz
