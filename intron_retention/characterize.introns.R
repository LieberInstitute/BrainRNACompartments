
# Check length


## Check evolutionary conservation via GERP score
2.	Conservation analysis: Base-wise Genomic Evolutionary Rate Profiling (GERP) scores on hg19 were downloaded from the GERP
website (http://mendel.stanford.edu/SidowLab/downloads/gerp/) (15, 16). For each OCR, we extracted the base-wise GERP scores 
in the region as well as its flanking regions of 1000bp (250bp and 2000bp in the stratified analysis). We then overlaid the 
extracted GERP score strings based on the center of each OCR, so that positions from different OCRs could be aligned to each 
other from centers to their proximities. We calculated the average GERP score for each aligned position and then smoothed the 
curve by applying a sliding window of 50bp. Means and 95% confidence intervals were calculated for each sliding window. 
Check UCSC genome browser up against down and matched non-regulated introns

http://compgen.cshl.edu/phast/faq.php
http://mendel.stanford.edu/SidowLab/downloads/gerp/

# Check RNA binding protein motifs

http://rbpmap.technion.ac.il/

# Repeat masker look for repeats

  http://www.repeatmasker.org/webrepeatmaskerhelp.html
  
# Intron position in transcript

ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("chromosome_name", "exon_chrom_start","exon_chrom_end",
                           "strand","ensembl_exon_id", "ensembl_transcript_id", "ensembl_gene_id"),
            mart=ensembl)
save(sym, file="/dcl01/lieber/ajaffe/Amanda/NucVsCyt/IRfinder_redo/biomaRt.transcript.annotation.rda")

geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]

