#LQTS_variants.txt is original file from APPRAISE paper set of variants; transcript table had to be made by me
python getlqtscoords.py
perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i varsets/lqtsvariants.vcf -o varsets/LQTSvariants.vcf --offline --fork 12 --force_overwrite --vcf
python lqtstranslate.py
python lqtsannotate.py
cd $HOME/software/VVP-pub
bash makeccr.sh
