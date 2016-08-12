rm results/*/*; rm plots/*/*
python plotdistro.py
python regression.py
bash analyze.sh 1 regions/topresid.txt patho.bed
bash analyze.sh 1 regions/midresid.txt patho.bed
bash analyze.sh 1 regions/topresid.txt benign.bed
bash analyze.sh 1 regions/midresid.txt benign.bed
bash plotloop.sh top
bash plotloop.sh mid
