#!/bin/bash
#SBATCH --account=ucgd-kp
#SBATCH --partition=ucgd-kp
#SBATCH -o %j-%N.out
#SBATCH -e %j-%N.err
perl software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i analysis/fake.vcf -o analysis/VEPfake.vcf --vcf --offline --plugin LoF --cache
