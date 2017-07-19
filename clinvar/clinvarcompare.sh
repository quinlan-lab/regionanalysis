date=20170710
set -e
if [[ ! -e clinvar_${date}.vcf.gz ]]; then
wget ftp://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/vcf_2.0/clinvar_${date}.vcf.gz
tabix -f clinvar_${date}.vcf.gz
fi
#makes expert pathogenic file and less confident "other" pathogenic files
python clinvarparse.py clinvar_${date}.vcf.gz
cat <(grep "^#" expert.clinvar.patho.vcf) <(grep -v "^#" expert.clinvar.patho.vcf | sort -k1,1 -k2,2n) | bgzip -c > expert.clinvar.patho.vcf.gz
cat <(grep "^#" other.clinvar.patho.vcf) <(grep -v "^#" other.clinvar.patho.vcf | sort -k1,1 -k2,2n) | bgzip -c > other.clinvar.patho.vcf.gz
tabix -f expert.clinvar.patho.vcf.gz; tabix -f other.clinvar.patho.vcf.gz
#tests them against CCR and creates weighted average scores for each variant in case of INDELs
python ../scorevars.py -x expert.clinvar.patho.vcf.gz -c ../gnomad-ccrs.bed.gz -a > ccrexp
python ../scorevars.py -x other.clinvar.patho.vcf.gz -c ../gnomad-ccrs.bed.gz -a > ccrother
echo "number of expert pathogenics overlapped SOLELY by ExAC v2 (gnomAD) variants"
grep "0.000" ccrexp | wc -l
echo "number of expert pathogenics NOT overlapped SOLELY by ExAC v2 (gnomAD) variants"
grep -v "0.000" ccrexp | wc -l
echo "number of single or multi submitter pathogenics overlapped SOLELY by ExAC v2 (gnomAD) variants"
grep "0.000" ccrother | wc -l
echo "number of single or multi submitter pathogenics NOT overlapped SOLELY by ExAC v2 (gnomAD) variants"
grep -v "0.000" ccrother | wc -l
echo "expert pathos covered by ExAC v2 INDELs"
CPI=$(bedtools intersect -a expert.clinvar.patho.vcf.gz -b <(zgrep "VARTRUE" ../gnomad-ccrs.bed.gz | awk '$3-$2>1') -sorted | wc -l)
echo $CPI
echo "expert pathos covered by ExAC v2 SNPs"
CPS=$(bedtools intersect -a expert.clinvar.patho.vcf.gz -b <(zgrep "VARTRUE" ../gnomad-ccrs.bed.gz | awk '$3-$2==1') -sorted | wc -l)
echo $CPS
echo "less confident pathos covered by ExAC v2 INDELs"
PI=$(bedtools intersect -a other.clinvar.patho.vcf.gz -b <(zgrep "VARTRUE" ../gnomad-ccrs.bed.gz | awk '$3-$2>1') -sorted | wc -l)
echo $PI
echo "less confident pathos covered by ExAC v2 SNPs"
PS=$(bedtools intersect -a other.clinvar.patho.vcf.gz -b <(zgrep "VARTRUE" ../gnomad-ccrs.bed.gz | awk '$3-$2==1') -sorted | wc -l)
echo $PS
echo "total proportion of pathogenics covered by ExAC v2 INDELs"
A=$(zgrep -v "^#" expert.clinvar.patho.vcf.gz | wc -l)
B=$(zgrep -v "^#" other.clinvar.patho.vcf.gz | wc -l)
bc <<< "scale=4; ($CPI+$PI)/($A+$B)"
echo "total proportion of pathogenics covered by ExAC v2 SNPs"
bc <<< "scale=4; ($CPS+$PS)/($A+$B)"
# adding in samocha for comparison temporarily
echo "samocha total proportion of pathogenics covered by ExAC v2 INDELs"
SPI=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz -b <(zgrep "VARTRUE" ../gnomad-ccrs.bed.gz | awk '$3-$2>1') -sorted | wc -l)
S=$(zgrep -v "^#"  ~/software/pathoscore/truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz | wc -l)
bc <<< "scale=4; $SPI/$S"
echo "samocha total proportion of pathogenics covered by ExAC v2 SNPs"
SPS=$(bedtools intersect -a ~/software/pathoscore/truth-sets/GRCh37/samocha/samocha.pathogenic.vcf.gz -b <(zgrep "VARTRUE" ../gnomad-ccrs.bed.gz | awk '$3-$2==1') -sorted | wc -l)
bc <<< "scale=4; $SPS/$S"
