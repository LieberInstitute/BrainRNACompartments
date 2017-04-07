library(ggplot2)
library(limma)

# Load IRFinder Results
IRres = list("1113C1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br1113C1/IRFinder-IR-nondir.txt", header = TRUE),
             "1113N1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br1113N1/IRFinder-IR-nondir.txt", header = TRUE),
             "2046C.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br2046C/IRFinder-IR-nondir.txt", header = TRUE),
             "2046N.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir.txt", header = TRUE),
             "2074C.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br2074C/IRFinder-IR-nondir.txt", header = TRUE),
             "2074N.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br2074N/IRFinder-IR-nondir.txt", header = TRUE),
             "5339C1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir.txt", header = TRUE),
             "5339N1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br5339N1/IRFinder-IR-nondir.txt", header = TRUE),
             "5340C1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br5340C1/IRFinder-IR-nondir.txt", header = TRUE),
             "5340N1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br5340N1/IRFinder-IR-nondir.txt", header = TRUE),
             "5341C1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br5341C1/IRFinder-IR-nondir.txt", header = TRUE),
             "5341N1.P" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/PolyA/Br5341N1/IRFinder-IR-nondir.txt", header = TRUE),
             "1113C1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br1113C1/IRFinder-IR-nondir.txt", header = TRUE),
             "2046C.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br2046C/IRFinder-IR-nondir.txt", header = TRUE),
             "2046N.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br2046N/IRFinder-IR-nondir.txt", header = TRUE),
             "2074C.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br2074C/IRFinder-IR-nondir.txt", header = TRUE),
             "2074N.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br2074N/IRFinder-IR-nondir.txt", header = TRUE),
             "5339C1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br5339C1/IRFinder-IR-nondir.txt", header = TRUE),
             "5339N1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br5339N1/IRFinder-IR-nondir.txt", header = TRUE),
             "5340C1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br5340C1/IRFinder-IR-nondir.txt", header = TRUE),
             "5340N1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br5340N1/IRFinder-IR-nondir.txt", header = TRUE),
             "5341C1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br5341C1/IRFinder-IR-nondir.txt", header = TRUE),
             "5341N1.R" = read.table("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/IRfinder/Ribozero/Br5341N1/IRFinder-IR-nondir.txt", header = TRUE))  
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings=="-"),]) 
# Associate the gene name with the intron
string = lapply(IRres, function(x) as.character(x$GeneIntronDetails))
string = lapply(string, function(x) strsplit(x, "/", fixed = TRUE))
x = lapply(string, function(x) unlist(x, recursive = FALSE))
a <- lapply(x, function (x) 1:length(x))
ge <- lapply(a, function(x) x[seq(2, length(x), 3)])
rownames = lapply(IRres, function(x) paste0("i",rownames(x)))
Gene=list()
for (i in 1:length(x)){
  tmp=x[[i]]
  Gene[[i]]=tmp[ge[[i]]]
  rownames(IRres[[i]]) = rownames[[i]]
}
names(Gene) = names(IRres)
IRres = Map(cbind, IRres, GeneID = Gene)
x = lapply(IRres, function(x) x$IRratio)
w = lapply(IRres, function(x) x$Warnings)

IRratios = data.frame(x[[1]], x[[2]], x[[3]], x[[4]], x[[5]], x[[6]], x[[7]], x[[8]], x[[9]], x[[10]], x[[11]], x[[12]], x[[13]], x[[14]], x[[15]], x[[16]], 
                      x[[17]], x[[18]], x[[19]], x[[20]], x[[21]], x[[22]], x[[23]], row.names = rownames(IRres[[1]]))
warnings = data.frame(w[[1]], w[[2]], w[[3]], w[[4]], w[[5]], w[[6]], w[[7]], w[[8]], w[[9]], w[[10]], w[[11]], w[[12]], w[[13]], w[[14]], w[[15]], w[[16]], 
                      w[[17]], w[[18]], w[[19]], w[[20]], w[[21]], w[[22]], w[[23]], row.names = rownames(IRres[[1]]))
colnames(IRratios) = colnames(warnings) = names(IRres)
IRratiosP = IRratios[,1:12]
IRratiosR = IRratios[,13:23]
wP = warnings[,1:12]
wR = warnings[,13:23]
wP = wP[which(wP$"1113C1.P"=="-" & wP$"1113N1.P"=="-" & wP$"2046C.P"=="-" & wP$"2046N.P"=="-" & wP$"2074C.P"=="-" & wP$"2074N.P"=="-" & 
               wP$"5339C1.P"=="-" & wP$"5339N1.P"=="-" & wP$"5340C1.P"=="-" & wP$"5340N1.P"=="-" & wP$"5341C1.P"=="-" & wP$"5341N1.P"=="-"),]
wR = wR[which(wR$"1113C1.R"=="-" & wR$"2046C.R"=="-" & wR$"2046N.R"=="-" & wR$"2074C.R"=="-" & wR$"2074N.R"=="-" & wR$"5339C1.R"=="-" & 
                wR$"5339N1.R"=="-" & wR$"5340C1.R"=="-" & wR$"5340N1.R"=="-" & wR$"5341C1.R"=="-" & wR$"5341N1.R"=="-"),]
IRratiosP = IRratiosP[which(rownames(IRratiosP) %in% rownames(wP)),] # 3264 introns remaining
IRratiosR = IRratiosR[which(rownames(IRratiosR) %in% rownames(wR)),] # 383 introns remaining
colnames(IRratiosP) = PolyA$SampleID
colnames(IRratiosR) = Ribozero$SampleID
IRratiosP = IRratiosP[which(rowSums(IRratiosP)!=0),] # 3111 introns remaining
IRratiosR = IRratiosR[which(rowSums(IRratiosR)!=0),] # 382 introns remaining
IR.AP = IRratiosP[,which(colnames(IRratiosP) %in% Adult.polyA$SampleID)]
IR.FP = IRratiosP[,which(colnames(IRratiosP) %in% Fetal.polyA$SampleID)]
IR.AP = IR.AP[which(rowSums(IR.AP)!=0),] # 2725 introns remaining
IR.FP = IR.FP[which(rowSums(IR.FP)!=0),] # 2942 introns remaining

## Differential expression
  #Adult PolyA
modAP = model.matrix(~Zone, data=Adult.polyA)
fitAP = lmFit(log2(IR.AP+1), modAP)
ebAP = ebayes(fitAP)
APstats = data.frame(log2FC = fitAP$coef[,2], pvalue = ebAP$p[,2])
APstats = APstats[order(APstats$pvalue),]
APstats = APstats[which(APstats$pvalue <=0.05),] # 331 introns

#Fetal PolyA
modFP = model.matrix(~Zone, data=Fetal.polyA)
fitFP = lmFit(log2(IR.FP+1), modFP)
ebFP = ebayes(fitFP)
FPstats = data.frame(log2FC = fitFP$coef[,2], pvalue = ebFP$p[,2])
FPstats = FPstats[order(FPstats$pvalue),]
FPstats = FPstats[which(FPstats$pvalue <=0.05),] # 216 introns

## run ANOVA
  # PolyA
# Percent Variance Explained Per Gene According to ANOVA 
explainListP = apply(IRratiosP, 1, function(y) {
  f = anova(lm(log2(y+1) ~ Fetal + Zone, data=PolyA))
  f[2]/sum(f[2])
})
percVarExplainMatP = matrix(unlist(explainListP), nc = 3, 
                           nrow = nrow(IRratiosP), byrow=TRUE, 
                           dimnames = list(rownames(IRratiosP), c("Fetal","Zone","Residuals")))

write.csv(percVarExplainMatP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/percVarExplained.IR.csv")
percVarExplainMatP = as.data.frame(percVarExplainMatP)
percVarExplainMatP = percVarExplainMatP[order(percVarExplainMatP$Zone, decreasing = TRUE),]

p = rbind(data.frame(Proportion = percVarExplainMatP[,1], Factor = "Age", Intron = rownames(percVarExplainMatP)), 
          data.frame(Proportion = percVarExplainMatP[,2], Factor = "Fraction", Intron = rownames(percVarExplainMatP)))

ggplot(p, aes(x=Proportion, fill=Factor, colour=Factor)) + geom_density(alpha=.3) +
  ggtitle("Proportion Individual Intron Retention Variance\nExplained By Each Factor") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))




