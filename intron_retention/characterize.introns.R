library(GenomicRanges)
library(data.table)
library(ggplot2)

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames
lapply(IRres, head)
# Filter introns
allIntrons = IRres[[1]]
allIntrons = allIntrons[grep("clean", allIntrons$GeneIntronDetails, fixed=T),]
string = unlist(strsplit(as.character(allIntrons$GeneIntronDetails), "/", fixed = TRUE), recursive = FALSE)
allIntrons = data.frame(allIntrons[,c(1:4,6)], genes = string[grep("ENSG", string)],
                        intronID = paste0(allIntrons$Chr,":",allIntrons$Start,"-",allIntrons$End))

## Differential IR by group
# read in results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
nonconst = list()
for (i in 1:length(comps)){
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
}
names(nonconst) = comps
elementNROWS(nonconst)
string = lapply(nonconst, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE))
IR.diff = lapply(nonconst, function(y) y$A.IRratio - y$B.IRratio)
dIR = Map(cbind, nonconst, ensID = lapply(string, function(y) y[grep("ENSG", y)]), 
          comments = lapply(string, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])), 
          IR.diff = IR.diff, Sign = lapply(IR.diff, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")), 
          intronID = lapply(nonconst, function(x) paste0(x$Chr,":",x$Start,"-",x$End)))
dIRclean = lapply(dIR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" &
                                           y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                                           y$comments=="clean"),])
sigdIR = lapply(dIRclean, function(x) x[which(x$p.diff<=0.05),])
sigdIR = lapply(sigdIR, function(x) split(x, x$Sign))
sigdIR = unlist(sigdIR, recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                  "Cytosol-Increased", "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")
lapply(sigdIR, head)
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)])) 
introns = c("All Introns" = list(allIntrons), sigdIR)

# Check length distribution of introns
introns = Map(cbind, introns, length = lapply(introns, function(x) abs(x$End - x$Start)))
introns = Map(cbind, introns[c(1:9, 11:13)], 
              Comparison = list("All Introns", "Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased",
                                "Prenatal:Nucleus-Increased","Cytosol:Adult-Increased","Cytosol:Prenatal-Increased",
                                "Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                                "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased"))
lapply(introns, head)





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

# Gene Ontology of genes containing differentially retained introns

