library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)
library(jaffelab)

# create ENCODE pd table
table = read.table("./Downloads/GSE30567_run_info.txt", header=T)
sc = scan("./Desktop/Untitled.txt", what="character")
table = table[which(table$SRA_Run %in% sc),]
dim(table)
nuc = grep("Nucleus", table$File_Uploaded_to_SRA)
cyt = grep("Cytosol", table$File_Uploaded_to_SRA)
table$CellType = table$Fraction = "NA"
table$Fraction[nuc] ="Nucleus"
table$Fraction[cyt] ="Cytosol"
table$CellType = gsub("hg19_wgEncodeCshlLongRnaSeq", "", table$File_Uploaded_to_SRA)
table$CellType = gsub("wgEncodeCshlLongRnaSeq", "", table$CellType)
table$CellType = gsub("\\Cytosol.*","", table$CellType)
table$CellType = gsub("\\Nucleus.*","", table$CellType)
write.table(table, file="./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/ENCODEpd.txt", quote=F)

enc = read.table("./Dropbox/sorted_figures/new/ENCODEpd.txt")
dim(enc)
head(enc)
length(unique(enc$SRA_Run))


### gene count files
setwd("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/Dougherty_data")
do = read.table("dougherty.pd.txt", header=T)
geneFn = paste0("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/Dougherty_data/Counts/Genes/",
                do$SRA_ID, "_Ensembl_v75_Genes.counts")
names(geneFn) = do$SRA_ID
all(file.exists(geneFn))

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
rownames(geneMap) = geneMap$Geneid
geneMap$Chr = paste0("chr", ss(geneMap$Chr, ";"))
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
geneMap$Geneid = NULL

## counts
geneCountList = mclapply(geneFn, function(x) {
  cat(".")
  read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,do$SRA_ID] # put in order

# make RPKM
bg = matrix(rep(pd$totalMapped), nc = nrow(pd), 
            nr = nrow(geneCounts),    byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
              nc = nrow(pd),    byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

### gene count files (ENCODE)
setwd("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/ENCODE_data/")
enc = read.table("ENCODEpd.txt")
names = unique(enc$SRA_Run)
geneFn = paste0("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/ENCODE_data/Counts/Genes/",
                names, "_Ensembl_v75_Genes.counts")
names(geneFn) = names
all(file.exists(geneFn))

### gene count files (Dougherty)
setwd("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/Dougherty_data")
do = read.table("dougherty.pd.txt", header=T)
names = unique(do$SRA_ID)
geneFn = paste0("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/Dougherty_data/Counts/Genes/",
                names, "_Ensembl_v75_Genes.counts")
names(geneFn) = names
all(file.exists(geneFn))


### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
rownames(geneMap) = geneMap$Geneid
geneMap$Chr = paste0("chr", ss(geneMap$Chr, ";"))
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
geneMap$Geneid = NULL

### biomart 
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
            values=rownames(geneMap), mart=ensembl)
geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]

## counts (ENCODE)
geneCountList = mclapply(geneFn, function(x) {
  cat(".")
  read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,enc$SRA_Run] # put in order

## counts (Dougherty)
geneCountList = mclapply(geneFn, function(x) {
  cat(".")
  read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,do$SRA_ID] # put in order


# make RPKM
bg = matrix(rep(pd$totalMapped), nc = nrow(pd), 
            nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
              nc = nrow(pd),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

# number of reads assigned
geneStatList = lapply(paste0(geneFn, ".summary"), 
                      read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = pd$RNum
pd$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

###############
### exon counts

exonFn = paste0("/dcl01/lieber/RNAseq/Datasets/DLPFC_Ribozero/Counts/Exons/DLPFC_Ribozero_",
                pd$RNum, "_Ensembl_v75_Exons.counts")
names(exonFn) = pd$RNum
all(file.exists(exonFn))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$Symbol = sym$hgnc_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Chr = paste0("chr", exonMap$Chr)

## counts
exonCountList = mclapply(exonFn, function(x) {
  cat(".")
  read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,pd$RNum] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

## make RPKM
bgE = matrix(rep(pd$totalMapped), nc = nrow(pd), 
             nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
              nc = nrow(pd),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)

#############
##### junctions

## via primary alignments only
junctionFiles = paste0("/dcl01/lieber/RNAseq/Datasets/DLPFC_Ribozero/Junctions/DLPFC_Ribozero_", 
                       pd$RNum, "_junctions_primaryOnly_regtools.count")
all(file.exists(junctionFiles)) #  TRUE

juncCounts = junctionCount(junctionFiles, pd$RNum,
                           strandSpecific = TRUE, maxCores=12)

## annotate junctions
load("/home/epi/ajaffe/Lieber/Projects/RNAseq/ensembl_hg19_v75_junction_annotation.rda")

anno = juncCounts$anno
seqlevels(anno, force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))

## add additional annotation
anno$inEnsembl = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inEnsemblStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inEnsemblEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$ensemblGeneID = NA
anno$ensemblGeneID[queryHits(oo)] = as.character(theJunctions$ensemblID[subjectHits(oo)])
anno$ensemblSymbol = NA
anno$ensemblSymbol[queryHits(oo)] = theJunctions$symbol[subjectHits(oo)]
anno$ensemblStrand = NA
anno$ensemblStrand[queryHits(oo)] = as.character(strand(theJunctions)[subjectHits(oo)])
anno$ensemblTx = CharacterList(vector("list", length(anno)))
anno$ensemblTx[queryHits(oo)] = theJunctions$tx[subjectHits(oo)]
anno$numTx = elementLengths(anno$ensemblTx)

# clean up
anno$ensemblSymbol = geneMap$Symbol[match(anno$ensemblGeneID, rownames(geneMap))]

## junction code
anno$code = ifelse(anno$inEnsembl, "InEns", 
                   ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
                          ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")))

## b/w exons and junctions
exonGR = GRanges( exonMap$Chr,	IRanges(exonMap$Start, exonMap$End))
anno$startExon = match(paste0(seqnames(anno),":",start(anno)-1), 
                       paste0(seqnames(exonGR), ":", end(exonGR)))
anno$endExon = match(paste0(seqnames(anno),":",end(anno)+1),
                     paste0(seqnames(exonGR), ":", start(exonGR)))
g = data.frame(leftGene = exonMap$Geneid[anno$startExon],
               rightGene = exonMap$Geneid[anno$endExon],
               leftGeneSym = exonMap$Symbol[anno$startExon],
               rightGeneSym = exonMap$Symbol[anno$endExon],
               stringsAsFactors=FALSE)
g$newGene = NA
g$newGeneSym = NA
g$newGene[which(g$leftGene==g$rightGene)] = 
  g$leftGene[which(g$leftGene==g$rightGene)] 
g$newGeneSym[which(g$leftGene==g$rightGene)] = 
  g$leftGeneSym[which(g$leftGene==g$rightGene)] 
g$newGene[which(g$leftGene!=g$rightGene)] = 
  paste0(g$leftGene,"-",g$rightGene)[which(g$leftGene!=g$rightGene)] 
g$newGeneSym[which(g$leftGene!=g$rightGene)] = 
  paste0(g$leftGeneSym,"-",g$rightGeneSym)[which(g$leftGene!=g$rightGene)] 
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] = 
  g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))] 
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] = 
  g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] = 
  g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] = 
  g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] 
g$newGeneSym[g$newGeneSym==""] = NA
g$newGeneSym[g$newGeneSym=="-"] = NA
anno$newGeneID = g$newGene
anno$newGeneSymbol = g$newGeneSym
anno$isFusion = grepl("-", anno$newGeneID)

## extract out
jMap = anno
jCounts = juncCounts$countDF
jCounts = jCounts[names(jMap),pd$RNum]

mappedPer80M = pd$totalMapped/80e6
countsM = DataFrame(mapply(function(x,d) x/d, jCounts , mappedPer80M))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm = as.data.frame(countsM)

## sequence of acceptor/donor sites
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
left = right = jMap
end(left) = start(left) +1
start(right) = end(right) -1

jMap$leftSeq  = getSeq(Hsapiens, left)
jMap$rightSeq = getSeq(Hsapiens, right)

### save counts
save(pd, jMap, jCounts, geneCounts, geneMap, exonCounts, exonMap, compress=TRUE,
     file="/dcl01/lieber/RNAseq/Datasets/DLPFC_Ribozero/rawCounts_DLPFC_Ribozero.rda")
save(pd, jMap, jRpkm, geneRpkm,	geneMap, exonRpkm, exonMap, compress=TRUE,
     file="/dcl01/lieber/RNAseq/Datasets/DLPFC_Ribozero/rpkmCounts_DLPFC_Ribozero.rda")

## write out for coverage
write.table(pd[,c("RNum", "bamFile")], 
            "samples_with_bams.txt",
            row.names=FALSE, quote=FALSE, sep="\t")