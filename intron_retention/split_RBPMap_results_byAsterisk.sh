filename=/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/RBPMap/RBPMap_results_hg19.txt

awk -vRS='\n\\*\\*+' '{print > "file" NR ".txt"}' $filename
