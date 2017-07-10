# converting Homsy files, filtering out non "D-missense variants"
tr -s "" "\n" < ogfiles/homsy_database_S02.txt | grep -v '\t\t\t' | cut -f 2- | sed '1,2d' | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | perl -pe 's/Start\s*lost/Start lost/g' | perl -pe 's/Splice\s*site/Splice site/g' | perl -pe 's/Stop\s*lost/Stop lost/g' | grep -v -Pw 'Splice site|Synonymous' | grep -w 'D-Missense' | perl -pe 's/\*.*//g' > ~/analysis/denovos/homsy.bed

tr -s "" "\n" < ogfiles/homsy_database_S03.txt | grep -v '\t\t\t' | cut -f 2- | sed '1,2d' | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | perl -pe 's/Start\s*lost/Start lost/g' | perl -pe 's/Splice\s*site/Splice site/g' | perl -pe 's/Stop\s*lost/Stop lost/g' | grep -v -Pw 'Splice site|Synonymous' | grep -w 'D-Missense' | perl -pe 's/\*.*//g' > ~/analysis/denovos/homsycontrol.bed

# converting deLigt files, filtering out synonymous variants
tr -s "" "\n" < ogfiles/deligtetal.txt | grep -v "p.(=)" | sed '1d' | grep -v '\t\t\t\t\t\t\t\t' | cut -f 4- | sed 's/^chr//g' | grep -w -v 'NO' > ~/analysis/denovos/deligt.bed

tr -s "" "\n" < ogfiles/deligtetal.txt | grep -v "p.(=)" | sed '1d' | grep -v '\t\t\t\t\t\t\t\t' | cut -f 4- | sed 's/^chr//g' | grep -w 'NO' > ~/analysis/denovos/deligtcontrol.bed

# converting McRae files, filtering out synonymous variants
tr -s "" "\n" < ogfiles/mcraeetal.txt | sed '1d' | cut -f 3- | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | awk '$8<.5' | grep -Pw 'stop_gained|stop_lost|start_lost|initiator_codon_variant|rare_amino_acid_variant|missense_variant|protein_altering_variant|frameshift_variant' > denovos/mcraecontrol.bed

tr -s "" "\n" < ogfiles/mcraeetal.txt | sed '1d' | cut -f 3- | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | awk '$8>.5' | grep -Pw 'stop_gained|stop_lost|start_lost|initiator_codon_variant|rare_amino_acid_variant|missense_variant|protein_altering_variant|frameshift_variant' > denovos/mcrae.bed

tr -s "" "\n" < ogfiles/mpcdenovos.txt | sed '1d' | cut -f -5,17 | awk '{print $1,$2,".",$3,$4,$5,$6}' OFS='\t' | cat mpcheader - > denovos/mpcdenovos.vcf 
