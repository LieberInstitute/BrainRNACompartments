cd /dcl01/lieber/ajaffe/Amanda/NucVsCyt/RBPMap

filename=/dcl01/lieber/ajaffe/Amanda/NucVsCyt/RBPMap/RBPMap_results_hg19.txt

awk '/^*/{i++};{print > (i".txt")}' $filename
