library(org.Hs.eg.db)
library(biomaRt)

load("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/rpkmCounts_nucleusVsCytosol_n24.rda")
load("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/rawCounts_nucleusVsCytosol_n24.rda")

keepPeople = which(pd$totalMapped > 1e6) # removes low sequenced adult nuclear ribozero sample
keepGenes = which(rowMeans(geneCounts) > 0)
x = which(rowMeans(geneRpkm) > 0)
geneRpkm = geneRpkm[x,keepPeople]
jRpkm = jRpkm[,keepPeople]
geneCounts = geneCounts[keepGenes,keepPeople]
exonCounts = exonCounts[,keepPeople]
pd = pd[keepPeople,]
geneMap = geneMap[keepGenes,]
JCOUNTS = data.frame(jCounts)
JCOUNTS = JCOUNTS[,keepPeople]
jCounts = jCounts[,keepPeople]
pd$Label = as.factor(paste(pd$Fetal, pd$Zone,pd$Library, sep="\n"))
pd$Label = factor(pd$Label, 
                  levels = c("Adult\nCytosol\npolyA", "Fetal\nCytosol\npolyA", "Adult\nNucleus\npolyA",
                             "Fetal\nNucleus\npolyA", "Adult\nCytosol\nRiboZero", "Fetal\nCytosol\nRiboZero",
                             "Adult\nNucleus\nRiboZero", "Fetal\nNucleus\nRiboZero"))
pd$WorkingID = c("Adult1_Cytosol_polyA", "Adult1_Nucleus_polyA", "Adult2_Cytosol_polyA", "Adult2_Nucleus_polyA",
                         "Adult3_Cytosol_polyA", "Adult3_Nucleus_polyA", "Fetal1_Cytosol_polyA", "Fetal1_Nucleus_polyA",
                         "Fetal2_Cytosol_polyA", "Fetal2_Nucleus_polyA", "Fetal3_Cytosol_polyA", "Fetal3_Nucleus_polyA",
                         "Adult1_Cytosol_RiboZero", "Adult2_Cytosol_RiboZero", "Adult2_Nucleus_RiboZero",
                         "Adult3_Cytosol_RiboZero", "Adult3_Nucleus_RiboZero", "Fetal1_Cytosol_RiboZero", "Fetal1_Nucleus_RiboZero",
                         "Fetal2_Cytosol_RiboZero", "Fetal2_Nucleus_RiboZero", "Fetal3_Cytosol_RiboZero", "Fetal3_Nucleus_RiboZero")

head(geneMap)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id", "hgnc_symbol",
                           "entrezgene", "description", "gene_biotype", "transcript_count"), 
            values=rownames(geneMap), mart=ensembl)
geneMap = data.frame(geneMap, EnsID = rownames(geneMap))
geneMap$Type = as.factor(sym$gene_biotype[match(geneMap$EnsID, sym$ensembl_gene_id)])
geneMap$TranscriptCount = sym$transcript_count[match(geneMap$EnsID, sym$ensembl_gene_id)]

#Parse pd into groups
Cytosol <- pd[which(pd$Zone=="Cytosol"),]
Nucleus <- pd[which(pd$Zone=="Nucleus"),]
Adult <- pd[which(pd$Fetal=="Adult"),]
Fetal <- pd[which(pd$Fetal=="Fetal"),]
PolyA <- pd[which(pd$Library=="polyA"),]
Ribozero <- pd[which(pd$Library=="RiboZero"),]

Adult.Ribo <-pd[which(pd$Fetal=="Adult" & pd$Library=="RiboZero"),]
Fetal.Ribo <-pd[which(pd$Fetal=="Fetal" & pd$Library=="RiboZero"),]
Adult.polyA <-pd[which(pd$Fetal=="Adult" & pd$Library=="polyA"),]
Fetal.polyA <-pd[which(pd$Fetal=="Fetal" & pd$Library=="polyA"),]
Cyt.Ribo <-  pd[which(pd$Zone=="Cytosol" & pd$Library=="RiboZero"),]
Nuc.Ribo <-pd[which(pd$Zone=="Nucleus" & pd$Library=="RiboZero"),]
Nuc.polyA <-pd[which(pd$Zone=="Nucleus"  & pd$Library=="polyA"),]
Cyt.polyA <-pd[which(pd$Zone=="Cytosol"  & pd$Library=="polyA"),]
Adult.cyt <-pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Adult"),]
Adult.nuc <-pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Adult"),]
Fetal.cyt <-pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Fetal"),]
Fetal.nuc <-pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Fetal"),]

ACP <- pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Adult" & pd$Library=="polyA"),]
ACR <- pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Adult" & pd$Library=="RiboZero"),]
ANP <- pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Adult" & pd$Library=="polyA"),]
ANR <- pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Adult" & pd$Library=="RiboZero"),]
FCP <- pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Fetal" & pd$Library=="polyA"),]
FCR <- pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Fetal" & pd$Library=="RiboZero"),]
FNP <- pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Fetal" & pd$Library=="polyA"),]
FNR <- pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Fetal" & pd$Library=="RiboZero"),]

# Parse count matrices
polya.counts <- geneCounts[,which(colnames(geneCounts)%in%PolyA$SampleID)]
ribozero.counts <- geneCounts[,which(colnames(geneCounts)%in%Ribozero$SampleID)]
Adult.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult$SampleID)]
Fetal.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal$SampleID)]
Cytosol.counts <- geneCounts[,which(colnames(geneCounts)%in%Cytosol$SampleID)]
Nucleus.counts <- geneCounts[,which(colnames(geneCounts)%in%Nucleus$SampleID)]

Adult.Ribo.counts <-geneCounts[,which(colnames(geneCounts)%in%Adult.Ribo$SampleID)]
Fetal.Ribo.counts <-geneCounts[,which(colnames(geneCounts)%in%Fetal.Ribo$SampleID)]
Adult.polyA.counts<-geneCounts[,which(colnames(geneCounts)%in%Adult.polyA$SampleID)]
Fetal.polyA.counts<-geneCounts[,which(colnames(geneCounts)%in%Fetal.polyA$SampleID)]
Cyt.Ribo.counts <-  geneCounts[,which(colnames(geneCounts)%in%Cyt.Ribo$SampleID)]
Nuc.Ribo.counts <-  geneCounts[,which(colnames(geneCounts)%in%Nuc.Ribo$SampleID)]
Nuc.polyA.counts <- geneCounts[,which(colnames(geneCounts)%in%Nuc.polyA$SampleID)]
Cyt.polyA.counts <- geneCounts[,which(colnames(geneCounts)%in%Cyt.polyA$SampleID)]
Adult.cyt.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult.cyt$SampleID)]
Adult.nuc.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult.nuc$SampleID)]
Fetal.cyt.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal.cyt$SampleID)]
Fetal.nuc.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal.nuc$SampleID)]