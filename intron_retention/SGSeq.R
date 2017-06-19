library("GenomicFeatures")
library("GenomicRanges")
library("SGSeq")

# http://www.bioconductor.org/packages/release/bioc/vignettes/SGSeq/inst/doc/SGSeq.html#overview

###  Create the TranscriptDb object from gtf file
# http://chitka-kalyan.blogspot.com/2014/02/creating-gencode-transcript-database-in.html

# Get the chromosome info as a dataframe
# One can use the script fethChromsomeSize from UCSC to the get info, filter it to remove information from non-std chromosomes

chrom.info=read.table(file="/users/aprice26/biotools/hg19.chrom.sizes", header=F)
colnames(chrom.info) = c("chrom", "length")
chrom.info$chrom = as.character(chrom.info$chrom)
chrom.info[which(chrom.info$chrom=="chr1_gl000192_random"),1] = "GL000192.1"
chrom.info[which(chrom.info$chrom=="chr4_gl000193_random"),1] = "GL000193.1"
chrom.info[which(chrom.info$chrom=="chr7_gl000195_random"),1] = "GL000195.1"
chrom.info[which(chrom.info$chrom=="chr9_gl000199_random"),1] = "GL000199.1"
chrom.info[which(chrom.info$chrom=="chr11_gl000202_random"),1] = "GL000202.1"
chrom.info[which(chrom.info$chrom=="chr17_gl000204_random"),1] = "GL000204.1"
chrom.info[which(chrom.info$chrom=="chr17_gl000205_random"),1] = "GL000205.1"
chrom.info[which(chrom.info$chrom=="chrUn_gl000212"),1] = "GL000212.1"
chrom.info[which(chrom.info$chrom=="chrUn_gl000220"),1] = "GL000220.1"
chrom.info[which(chrom.info$chrom=="chrUn_gl000228"),1] = "GL000228.1"
chrom.info[which(chrom.info$chrom=="chrUn_gl000241"),1] = "GL000241.1"

# Create the transcriptdb from the gencode gtf file
# Download the latest gencode comprehensive gtf file from gencode website

gencode <- makeTxDbFromGFF("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf", 
                           chrominfo=chrom.info, 
                           dataSource=paste("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz"),
                           organism="Homo sapiens")

# Save the transcriptDb object as a sql database object
saveDb(gencode, file="/dcl01/lieber/ajaffe/Amanda/NucVsCyt/MISO_out/annotation/hg19/gencode.v25lift37.annotation.sqlite")
loadDb("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/MISO_out/annotation/hg19/gencode.v25lift37.annotation.sqlite")


