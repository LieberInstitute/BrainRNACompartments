# For whatever reason, this only works on server 4 (and not on the biostat cluster or my LIBD workstation)
# in the directory /media/DATA/Amanda/RBPMap, run 'sh parser.sh'

cd /media/DATA/Amanda/RBPMap

filename=/media/DATA/Amanda/RBPMap/RBPMap_results_hg19.txt

awk '/^*/{i++};{print > (i".txt")}' $filename
