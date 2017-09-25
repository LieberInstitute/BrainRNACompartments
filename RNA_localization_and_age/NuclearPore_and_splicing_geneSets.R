library(readxl)
library(DESeq2)
library("ggplot2")

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Nucleoporin genes from http://amigo.geneontology.org/grebe
nups = read.table("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/nuclear_pore_genes.txt", sep = "\t", header = F)
nupGenes = c(as.character(unique(nups$V3)), "NUPL1", "POM121C")
geneMap[which(geneMap$ensemblID=="ENSG00000272391"),"Symbol"] = "POM121C"
nups = geneMap[which(geneMap$Symbol %in% nupGenes),]
nupCounts = geneCounts.down[which(rownames(geneCounts.down) %in% rownames(nups)),grep("polyA", colnames(geneCounts.down))]

# Splicing genes
splicepos = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splicingPositive.xlsx", col_names = FALSE))
splicing = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splicing.xlsx", col_names = FALSE))
spliceneg = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splicingNegative.xlsx", col_names = FALSE))
spl = rbind(splicepos, splicing, spliceneg)
geneMap[which(geneMap$ensemblID=="ENSG00000255663"),"Symbol"] = "I3L521"
geneMap[which(geneMap$ensemblID=="ENSG00000265241"),"Symbol"] = "RBM8A"
splGenes = c(as.character(unique(spl$X__3)), "C11orf35", "ZNF259")
spl = geneMap[which(geneMap$Symbol %in% splGenes),]
splCounts = geneCounts.down[which(rownames(geneCounts.down) %in% rownames(spl)),grep("polyA", colnames(geneCounts.down))]

# Differential expression of 2 gene sets
Nupdds.int <- DESeqDataSetFromMatrix(countData = nupCounts, colData = pd[which(pd$Library=="polyA"),], design = ~ Zone + Fetal + Fetal:Zone)
Nupdds.int <- DESeq(Nupdds.int)
Nupres.int = results(Nupdds.int)
Nupdds.age <- DESeqDataSetFromMatrix(countData = nupCounts, colData = pd[which(pd$Library=="polyA"),], design = ~ Zone + Fetal)
Nupdds.age <- DESeq(Nupdds.age)
Nupres.age = results(Nupdds.age)
Nupdds.frac <- DESeqDataSetFromMatrix(countData = nupCounts, colData = pd[which(pd$Library=="polyA"),], design = ~ Fetal + Zone)
Nupdds.frac <- DESeq(Nupdds.frac)
Nupres.frac = results(Nupdds.frac)

Nupres.int = cbind(Nupres.int, geneMap[match(rownames(Nupres.int), geneMap$gencodeID),])
write.csv(Nupres.int, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/NuclearPoreGene_DE_Age_Frac_interaction.csv")
Nupres.age = cbind(Nupres.age, geneMap[match(rownames(Nupres.age), geneMap$gencodeID),])
write.csv(Nupres.age, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/NuclearPoreGene_DE_Age.csv")
Nupres.frac = cbind(Nupres.frac, geneMap[match(rownames(Nupres.frac), geneMap$gencodeID),])
write.csv(Nupres.frac, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/NuclearPoreGene_DE_Frac.csv")

nupres = list(Age = Nupres.age[which(Nupres.age$padj<=0.05),], Fraction = Nupres.frac[which(Nupres.frac$padj<=0.05),], 
              Interaction = Nupres.int[which(Nupres.int$padj<=0.05),])
nupres = lapply(nupres, function(x) x[order(abs(x$log2FoldChange), decreasing =T),])

# plot expression of nuclear pore genes significantly DE by age, fraction, or an interaction between the two
elementNROWS(nupres)
rpkm = lapply(nupres, function(x) t(geneRpkm.down[which(rownames(geneRpkm.down) %in% rownames(x)), grep("polyA", colnames(geneRpkm.down))]))
rpkm = lapply(rpkm, as.data.frame)
rpkm = Map(cbind, rpkm, Age = list(c(rep.int("Adult", 6), rep.int("Prenatal",6)),c(rep.int("Adult", 6), rep.int("Prenatal",6)),c(rep.int("Adult", 6), rep.int("Prenatal",6))),
           Fraction = list(rep.int(c("Cytosol", "Nucleus"), 6), rep.int(c("Cytosol", "Nucleus"), 6), rep.int(c("Cytosol", "Nucleus"), 6)))
geneIDs = lapply(rpkm, colnames)

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/NuclearPoreGene_DE_polyAonly.pdf", height = 6, width = 6)
frac = rpkm[["Fraction"]]
for (j in 1:(ncol(frac)-2)){
  tmp = frac[,c(colnames(frac)[j],"Age","Fraction")]
  g = ggplot(tmp, aes(x=Age, y=tmp[,1], fill=Fraction)) + geom_boxplot() +
      ggtitle(paste0("Fraction: ", geneMap[which(geneMap$gencodeID==colnames(frac)[j]),"Symbol"],
                     "\nFDR= ", round(nupres[["Fraction"]][which(rownames(nupres[["Fraction"]])==colnames(frac)[j]),"padj"],3))) +
      ylab("RPKM") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)
}
age = rpkm[["Age"]]
for (j in 1:(ncol(age)-2)){
  tmp = age[,c(colnames(age)[j],"Age","Fraction")]
  g = ggplot(tmp, aes(x=Age, y=tmp[,1], fill=Fraction)) + geom_boxplot() +
    ggtitle(paste0("Age: ", geneMap[which(geneMap$gencodeID==colnames(age)[j]),"Symbol"],
                   "\nFDR= ", round(nupres[["Age"]][which(rownames(nupres[["Age"]])==colnames(age)[j]),"padj"],3))) +
    ylab("RPKM") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)
}
int = rpkm[["Interaction"]]
for (j in 1:(ncol(int)-2)){
  tmp = int[,c(colnames(int)[j],"Age","Fraction")]
  g = ggplot(tmp, aes(x=Age, y=tmp[,1], fill=Fraction)) + geom_boxplot() +
    ggtitle(paste0("Interaction: ", geneMap[which(geneMap$gencodeID==colnames(int)[j]),"Symbol"],
                   "\nFDR= ", round(nupres[["Interaction"]][which(rownames(nupres[["Interaction"]])==colnames(int)[j]),"padj"],3))) +
    ylab("RPKM") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)
}
dev.off()


# Splicing genes
spldds.int <- DESeqDataSetFromMatrix(countData = splCounts, colData = pd[which(pd$Library=="polyA"),], design = ~ Zone + Fetal + Fetal:Zone)
spldds.int <- DESeq(spldds.int)
splres.int = results(spldds.int)
spldds.age <- DESeqDataSetFromMatrix(countData = splCounts, colData = pd[which(pd$Library=="polyA"),], design = ~ Zone + Fetal)
spldds.age <- DESeq(spldds.age)
splres.age = results(spldds.age)
spldds.frac <- DESeqDataSetFromMatrix(countData = splCounts, colData = pd[which(pd$Library=="polyA"),], design = ~ Fetal + Zone)
spldds.frac <- DESeq(spldds.frac)
splres.frac = results(spldds.frac)

Nupres.int = cbind(Nupres.int, geneMap[match(rownames(Nupres.int), geneMap$gencodeID),])
write.csv(Nupres.int, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/NuclearPoreGene_DE_Age_Frac_interaction.csv")
Nupres.age = cbind(Nupres.age, geneMap[match(rownames(Nupres.age), geneMap$gencodeID),])
write.csv(Nupres.age, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/NuclearPoreGene_DE_Age.csv")
Nupres.frac = cbind(Nupres.frac, geneMap[match(rownames(Nupres.frac), geneMap$gencodeID),])
write.csv(Nupres.frac, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/NuclearPoreGene_DE_Frac.csv")

splres.int <- cbind(splres.int, geneMap[match(rownames(splres.int), geneMap$gencodeID),])
write.csv(splres.int, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/SplicingGenes_DE_Age_Frac_interaction.csv")
splres.age <- cbind(splres.age, geneMap[match(rownames(splres.age), geneMap$gencodeID),])
write.csv(splres.age, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/SplicingGenes_DE_Age.csv")
splres.frac <- cbind(splres.frac, geneMap[match(rownames(splres.frac), geneMap$gencodeID),])
write.csv(splres.frac, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/SplicingGenes_DE_Age.csv")

splres = list(Age = splres.age[which(splres.age$padj<=0.05),], Fraction = splres.frac[which(splres.frac$padj<=0.05),], 
              Interaction = splres.int[which(splres.int$padj<=0.05),])
splres = lapply(splres, function(x) x[order(abs(x$log2FoldChange), decreasing =T),])

# plot expression of splicing genes significantly DE by age, fraction, or an interaction between the two
elementNROWS(splres)
rpkm = lapply(splres, function(x) t(geneRpkm.down[which(rownames(geneRpkm.down) %in% rownames(x)), grep("polyA", colnames(geneRpkm.down))]))
rpkm = lapply(rpkm, as.data.frame)
rpkm = Map(cbind, rpkm, Age = list(c(rep.int("Adult", 6), rep.int("Prenatal",6)),c(rep.int("Adult", 6), rep.int("Prenatal",6)),c(rep.int("Adult", 6), rep.int("Prenatal",6))),
           Fraction = list(rep.int(c("Cytosol", "Nucleus"), 6), rep.int(c("Cytosol", "Nucleus"), 6), rep.int(c("Cytosol", "Nucleus"), 6)))
geneIDs = lapply(rpkm, colnames)

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/SplicingGenes_DE_polyAonly.pdf", height = 6, width = 6)
frac = rpkm[["Fraction"]]
for (j in 1:(ncol(frac)-2)){
  tmp = frac[,c(colnames(frac)[j],"Age","Fraction")]
  g = ggplot(tmp, aes(x=Age, y=tmp[,1], fill=Fraction)) + geom_boxplot() +
    ggtitle(paste0("Fraction: ", geneMap[which(geneMap$gencodeID==colnames(frac)[j]),"Symbol"],
                   "\nFDR= ", round(splres[["Fraction"]][which(rownames(splres[["Fraction"]])==colnames(frac)[j]),"padj"],3))) +
    ylab("RPKM") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)
}
age = rpkm[["Age"]]
for (j in 1:(ncol(age)-2)){
  tmp = age[,c(colnames(age)[j],"Age","Fraction")]
  g = ggplot(tmp, aes(x=Age, y=tmp[,1], fill=Fraction)) + geom_boxplot() +
    ggtitle(paste0("Age: ", geneMap[which(geneMap$gencodeID==colnames(age)[j]),"Symbol"],
                   "\nFDR= ", round(splres[["Age"]][which(rownames(splres[["Age"]])==colnames(age)[j]),"padj"],3))) +
    ylab("RPKM") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)
}
int = rpkm[["Interaction"]]
for (j in 1:(ncol(int)-2)){
  tmp = int[,c(colnames(int)[j],"Age","Fraction")]
  g = ggplot(tmp, aes(x=Age, y=tmp[,1], fill=Fraction)) + geom_boxplot() +
    ggtitle(paste0("Interaction: ", geneMap[which(geneMap$gencodeID==colnames(int)[j]),"Symbol"],
                   "\nFDR= ", round(splres[["Interaction"]][which(rownames(splres[["Interaction"]])==colnames(int)[j]),"padj"],3))) +
    ylab("RPKM") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)
}
dev.off()


### look at these gene sets in the global results: Nuclear pore genes

## Are nuclear pore genes more likely to be significantly DE than other genes?
# By fraction in adult
fisher.test(data.frame(c(nrow(Apres[which(rownames(Apres) %in% nups$gencodeID & Apres$padj<=0.05),]),
                         nrow(Apres[which(rownames(Apres) %in% nups$gencodeID & Apres$padj>0.05),])),
                       c(nrow(Apres[c(!which(rownames(Apres) %in% nups$gencodeID),which(Apres$padj<=0.05)),]),
                         nrow(Apres[c(!which(rownames(Apres) %in% nups$gencodeID),which(Apres$padj>0.05)),]))))
#p-value = 0.7909
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6064772 1.8155422
#sample estimates:
#  odds ratio 
#1.069531
            
# By fraction in prenatal
fisher.test(data.frame(c(nrow(Fpres[which(rownames(Fpres) %in% nups$gencodeID & Fpres$padj<=0.05),]),
                         nrow(Fpres[which(rownames(Fpres) %in% nups$gencodeID & Fpres$padj>0.05),])),
                       c(nrow(Fpres[c(!which(rownames(Fpres) %in% nups$gencodeID),which(Fpres$padj<=0.05)),]),
                         nrow(Fpres[c(!which(rownames(Fpres) %in% nups$gencodeID),which(Fpres$padj>0.05)),]))))
# p-value = 0.4158
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.000000 2.248247
#sample estimates:
#  odds ratio 
#0

# By age in cytosol (more are sig than expected)
fisher.test(data.frame(c(nrow(Cpres[which(rownames(Cpres) %in% nups$gencodeID & Cpres$padj<=0.05),]),
                         nrow(Cpres[which(rownames(Cpres) %in% nups$gencodeID & Cpres$padj>0.05),])),
                       c(nrow(Cpres[c(!which(rownames(Cpres) %in% nups$gencodeID),which(Cpres$padj<=0.05)),]),
                         nrow(Cpres[c(!which(rownames(Cpres) %in% nups$gencodeID),which(Cpres$padj>0.05)),]))))
#p-value = 0.02292
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.058867 2.791952
#sample estimates:
#  odds ratio 
#1.706797

# By age in nucleus (more are sig than expected)
fisher.test(data.frame(c(nrow(Npres[which(rownames(Npres) %in% nups$gencodeID & Npres$padj<=0.05),]),
                         nrow(Npres[which(rownames(Npres) %in% nups$gencodeID & Npres$padj>0.05),])),
                       c(nrow(Npres[c(!which(rownames(Npres) %in% nups$gencodeID),which(Npres$padj<=0.05)),]),
                         nrow(Npres[c(!which(rownames(Npres) %in% nups$gencodeID),which(Npres$padj>0.05)),]))))
#p-value = 0.03913
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9996858 2.5719497
#sample estimates:
#  odds ratio 
#1.599729


## In significantly DE genes, are they expressed in the same direction as other genes?
# By fraction in adult
fisher.test(data.frame(c(nrow(Apres[which(rownames(Apres) %in% nups$gencodeID & Apres$padj<=0.05 & Apres$log2FoldChange>0),]),
                         nrow(Apres[which(rownames(Apres) %in% nups$gencodeID & Apres$padj<=0.05 & Apres$log2FoldChange<0),])),
                       c(nrow(Apres[c(!which(rownames(Apres) %in% nups$gencodeID),which(Apres$padj<=0.05 & Apres$log2FoldChange>0)),]),
                         nrow(Apres[c(!which(rownames(Apres) %in% nups$gencodeID),which(Apres$padj<=0.05 & Apres$log2FoldChange<0)),]))))
#p-value = 0.8228
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2948002 2.2210082
#sample estimates:
#  odds ratio 
#0.833043

# By fraction in prenatal
fisher.test(data.frame(c(nrow(Fpres[which(rownames(Fpres) %in% nups$gencodeID & Fpres$padj<=0.05 & Fpres$log2FoldChange>0),]),
                         nrow(Fpres[which(rownames(Fpres) %in% nups$gencodeID & Fpres$padj<=0.05 & Fpres$log2FoldChange<0),])),
                       c(nrow(Fpres[c(!which(rownames(Fpres) %in% nups$gencodeID),which(Fpres$padj<=0.05 & Fpres$log2FoldChange>0)),]),
                         nrow(Fpres[c(!which(rownames(Fpres) %in% nups$gencodeID),which(Fpres$padj<=0.05 & Fpres$log2FoldChange<0)),]))))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0

# By age in cytosol (more are higher expressed in prenatal than expected)
fisher.test(data.frame(c(nrow(Cpres[which(rownames(Cpres) %in% nups$gencodeID & Cpres$padj<=0.05 & Cpres$log2FoldChange>0),]),
                         nrow(Cpres[which(rownames(Cpres) %in% nups$gencodeID & Cpres$padj<=0.05 & Cpres$log2FoldChange<0),])),
                       c(nrow(Cpres[c(!which(rownames(Cpres) %in% nups$gencodeID),which(Cpres$padj<=0.05 & Cpres$log2FoldChange>0)),]),
                         nrow(Cpres[c(!which(rownames(Cpres) %in% nups$gencodeID),which(Cpres$padj<=0.05 & Cpres$log2FoldChange<0)),]))))
#p-value = 1.258e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.771662 16.511760
#sample estimates:
#  odds ratio 
#6.24588

# By age in nucleus (more are higher expressed in prenatal than expected)
fisher.test(data.frame(c(nrow(Npres[which(rownames(Npres) %in% nups$gencodeID & Npres$padj<=0.05 & Npres$log2FoldChange>0),]),
                         nrow(Npres[which(rownames(Npres) %in% nups$gencodeID & Npres$padj<=0.05 & Npres$log2FoldChange<0),])),
                       c(nrow(Npres[c(!which(rownames(Npres) %in% nups$gencodeID),which(Npres$padj<=0.05 & Npres$log2FoldChange>0)),]),
                         nrow(Npres[c(!which(rownames(Npres) %in% nups$gencodeID),which(Npres$padj<=0.05 & Npres$log2FoldChange<0)),]))))
#p-value = 2.973e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.38324 14.56230
#sample estimates:
#  odds ratio 
#5.451623