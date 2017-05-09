library(GenomicRanges)
library(ggplot2)

###############################################################
####################### intron position #######################
###############################################################

### Check intron position in the transcript

ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("chromosome_name", "exon_chrom_start","exon_chrom_end",
                           "strand","ensembl_exon_id", "ensembl_transcript_id", "ensembl_gene_id"),
            mart=ensembl)
save(sym, file="/dcl01/lieber/ajaffe/Amanda/NucVsCyt/IRfinder_redo/biomaRt.transcript.annotation.rda")

geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]