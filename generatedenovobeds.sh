tr -s "" "\n" < ogfiles/homsy_database_S02.txt | grep -v '\t\t\t' | cut -f 2- | sed '1d' | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | perl -pe 's/Start\s*lost/Start lost/g' | perl -pe 's/Splice\s*site/Splice site/g' | perl -pe 's/Stop\s*lost/Stop lost/g' | grep -v 'Splice site' | perl -pe 's/\*.*//g' | grep 'D-Missense' > ~/analysis/denovos/homsy.bed

tr -s "" "\n" < ogfiles/homsy_database_S03.txt | grep -v '\t\t\t' | cut -f 2- | sed '1d' | awk '{t=$2-1} {printf $1"\t"t"\t"$2} {for (i=3; i<=NF; i++) printf "\t"$i} {print ""}' | perl -pe 's/Start\s*lost/Start lost/g' | perl -pe 's/Splice\s*site/Splice site/g' | perl -pe 's/Stop\s*lost/Stop lost/g' | grep -v 'Splice site' |  perl -pe 's/\*.*//g' | grep 'D-Missense' > ~/analysis/denovos/homsycontrol.bed


