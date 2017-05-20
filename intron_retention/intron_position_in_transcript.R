library(GenomicRanges)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

### read in Differential IR results files
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
          intronID = lapply(nonconst, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End)))
dIRclean = lapply(dIR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" &
                                             y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                                             y$comments=="clean"),])
sigdIR = lapply(dIRclean, function(x) x[which(x$p.diff<=0.05),])
sigdIR = lapply(sigdIR, function(x) split(x, x$Sign))
sigdIR = unlist(sigdIR, recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                  "Cytosol-Increased", "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)]))
introns = c("Introns (Fraction)" = list(do.call(rbind, lapply(dIRclean[1:2], function(x) x[,c(1:4,6,8,30,32,33,34)]))), 
            "Introns (Age)" = list(do.call(rbind, lapply(dIRclean[3:4], function(x) x[,c(1:4,6,8,30,32,33,34)]))), sigdIR)
introns[["Introns (Fraction)"]] = introns[["Introns (Fraction)"]][unique(introns[["Introns (Fraction)"]][,"intronID"]),]
introns[["Introns (Age)"]] = introns[["Introns (Age)"]][unique(introns[["Introns (Age)"]][,"intronID"]),]



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