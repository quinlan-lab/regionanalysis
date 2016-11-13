tr -s "" "\n" < ogfiles/VVP_KCNH2.txt | sed '1d' | awk '{printf $1"\t"$2-1"\t"} {for (i=2;i<=NF;i++) printf $i"\t"; printf "\n"}' > VVPvariants/KCNH2variants.bed
