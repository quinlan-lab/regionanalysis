rm results/*/*; rm plots/*/*
#python plotdistro.py
#python regression.py
bash analyze.sh 1 regions/topresid.txt denovos/homsy.bed
bash analyze.sh 1 regions/midresid.txt denovos/homsy.bed
bash analyze.sh 1 regions/topresid.txt denovos/homsycontrol.bed
bash analyze.sh 1 regions/midresid.txt denovos/homsycontrol.bed
bash plotloop.sh top
bash plotloop.sh mid
