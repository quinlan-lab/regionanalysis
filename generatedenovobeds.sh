# converting Homsy files, filtering out non "D-missense variants"
tr -s "" "\n" < ogfiles/homsy_database_S02.txt | grep -v '\t\t\t' | cut -f 2- | sed '1d' | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | perl -pe 's/Start\s*lost/Start lost/g' | perl -pe 's/Splice\s*site/Splice site/g' | perl -pe 's/Stop\s*lost/Stop lost/g' | grep -v 'Splice site' | perl -pe 's/\*.*//g' | grep 'D-Missense' > ~/analysis/denovos/homsy.bed

tr -s "" "\n" < ogfiles/homsy_database_S03.txt | grep -v '\t\t\t' | cut -f 2- | sed '1d' | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | perl -pe 's/Start\s*lost/Start lost/g' | perl -pe 's/Splice\s*site/Splice site/g' | perl -pe 's/Stop\s*lost/Stop lost/g' | grep -v 'Splice site' |  perl -pe 's/\*.*//g' | grep 'D-Missense' > ~/analysis/denovos/homsycontrol.bed

# converting deLigt files, filtering out synonymous variants
tr -s "" "\n" < ogfiles/deligtetal.txt | grep -v "p.(=)" | sed '1d' | grep -v '\t\t\t\t\t\t\t\t' | cut -f 4- | sed 's/^chr//g' | grep -w -v 'NO' > ~/analysis/denovos/deligt.bed

tr -s "" "\n" < ogfiles/deligtetal.txt | grep -v "p.(=)" | sed '1d' | grep -v '\t\t\t\t\t\t\t\t' | cut -f 4- | sed 's/^chr//g' | grep -w 'NO' > ~/analysis/denovos/deligtcontrol.bed

# converting McRae files, filtering out synonymous variants
tr -s "" "\n" < ogfiles/mcraeetal.txt | sed '1d' | cut -f 3- | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | awk '$8<.5' | grep -v -Pw "synonymous_variant|splice_donor_variant|splice_acceptor_variant|splice_region_variant" > denovos/mcraecontrol.bed

tr -s "" "\n" < ogfiles/mcraeetal.txt | sed '1d' | cut -f 3- | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | awk '$8>.5' | grep -v -Pw "synonymous_variant|splice_donor_variant|splice_acceptor_variant|splice_region_variant" > denovos/mcrae.bed
