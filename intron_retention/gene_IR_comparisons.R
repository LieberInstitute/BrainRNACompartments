library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Load IRFinder Results
names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames

# Filter introns
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings!="NonUniformIntronCover"),])
IRfiltered = lapply(IRfiltered, function(x) x[grep("clean", x$GeneIntronDetails, fixed=T),])
IRfiltered = lapply(IRfiltered, function(x) x[max(x$SplicesExact,x$SplicesRight,x$SplicesLeft)>4 | (x$ExonToIntronReadsLeft>4 & x$ExonToIntronReadsRight>4),])
IRfiltered = Map(cbind, IRfiltered, 
                 genes = lapply(lapply(IRfiltered, function(x) unlist(strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE), recursive = F)), function(x) x[grep("ENSG", x)]), 
                 intronID = lapply(IRfiltered, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
elementNROWS(IRfiltered)
head(IRfiltered[[1]]) # all introns passing QC per sample
irbyGene = lapply(IRfiltered, data.table)
irbyGene = lapply(irbyGene,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="genes"])) # maximum IR for each gene with at least one intron passing QC

#Number of genes with introns passing filtering steps (not overlapping anything, enough coverage, etc.)
elementNROWS(irbyGene)
elementNROWS(lapply(IRfiltered, function(x) unique(x$genes)))
#Number of Genes in sig gene set
elementNROWS(sig)

# Get IR info for genes significantly regulated by fraction 
max = do.call(rbind, irbyGene)
max$SampleID = gsub("\\..*","", rownames(max))
max$rownum = c(1:nrow(max))
max$Fraction = ifelse(max$rownum %in% grep("C", max$SampleID), "Cytosol", "Nucleus")
max$Age = ifelse(max$rownum %in% grep("Br53", max$SampleID), "Prenatal", "Adult")
max$Group = paste(max$Age, max$Fraction, sep=":")
sigIR = lapply(sig, function(x) max[which(max$genes %in% x$ensID),])
unlist(lapply(sigIR, function(x) length(unique(x$genes))))

# genes >50% or >10% retained in at least one sample
perc = lapply(sigIR, function(x) list(sig.more50=as.character(unique(x[which(x$IRratio>=0.50),"genes"])),
                                      sig.less50=as.character(unique(x[which(x$IRratio<0.50),"genes"])),
                                      all.more50=as.character(unique(max[which(max$IRratio>=0.50),"genes"])),
                                      all.less50=as.character(unique(max[which(max$IRratio<0.50),"genes"])),
                                      sig.more10=as.character(unique(x[which(x$IRratio>=0.10),"genes"])),
                                      sig.less10=as.character(unique(x[which(x$IRratio<0.10),"genes"])),
                                      all.more10=as.character(unique(max[which(max$IRratio>=0.10),"genes"])),
                                      all.less10=as.character(unique(max[which(max$IRratio<0.10),"genes"]))))
for (i in 1:length(perc)){
venn.diagram(perc[[i]][1:4], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",
                               names(perc)[i], ".50.percent.IRratio.jpeg"), 
             main=names(perc)[i], col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"), alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 50.percent.IRratio.overlap_byFraction.pdf

for (i in 1:length(perc)){
  venn.diagram(perc[[i]][5:8], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",
                                      names(perc)[i], ".10.percent.IRratio.jpeg"), 
               main=names(perc)[i], col = "transparent", 
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"), alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold", 
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 10.percent.IRratio.overlap_byFraction.pdf


# test whether genes with >50% or 10% intron retention in any sample are preferentially significantly DEG in the pattern listed

moreless = list()
for (i in 1:length(perc)){
  tmp = perc[[i]]
  moreless[[i]] = list("50%" = data.frame(greater = c(length(unique(tmp$sig.more50)),length(unique(tmp$all.more50[!(tmp$all.more50 %in% tmp$sig.more50)]))),
                                          less = c(length(unique(tmp$sig.less50)),length(unique(tmp$all.less50[!(tmp$all.less50 %in% tmp$sig.less50)]))), 
                                          row.names = c("sig","NS")),
                       "10%" = data.frame(greater = c(length(unique(tmp$sig.more10)),length(unique(tmp$all.more10[!(tmp$all.more10 %in% tmp$sig.more10)]))),
                                          less = c(length(unique(tmp$sig.less10)),length(unique(tmp$all.less10[!(tmp$all.less10 %in% tmp$sig.less10)]))), 
                                          row.names = c("sig","NS")))
}
names(moreless) = names(perc)
moreless.fisher = lapply(moreless, function(x) lapply(x,fisher.test))
write.csv(rbind(pval = data.frame(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/fisher_50perc_10perc_IR_sigDEGs_relationship_byFraction.csv")


# Intron retention between retained and exported genes
IR_ret_exp = list(both = t.test(sigIR[["both_retained"]][,"IRratio"], sigIR[["both_exported"]][,"IRratio"]),
                  adult = t.test(c(sigIR[["both_retained"]][,"IRratio"],sigIR[["Ad_retained"]][,"IRratio"]), 
                                 c(sigIR[["both_exported"]][,"IRratio"],sigIR[["Ad_exported"]][,"IRratio"])),
                  prenatal = t.test(c(sigIR[["both_retained"]][,"IRratio"],sigIR[["Fet_retained"]][,"IRratio"]), 
                                    c(sigIR[["both_exported"]][,"IRratio"],sigIR[["Fet_exported"]][,"IRratio"])))
write.csv(rbind(Tstat = data.frame(lapply(IR_ret_exp, function(x) round(x$statistic,3))), pval = data.frame(lapply(IR_ret_exp, function(x) x$p.value)),
                confInt = data.frame(lapply(IR_ret_exp, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(IR_ret_exp, function(x) round(x$estimate,3)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/Ttest_exported_vs_retained_DEG_IRratios.csv")


# What does the distribution of IR Ratios by Sig group look like?
gr = c("Both Retained", "Both Exported", "Retained in Prenatal", "Retained in Adult",
       "Exported in Prenatal", "Exported in Adult", "Retained in Adult/\nExported in Prenatal",
       "Retained in Prenatal/\nExported in Adult", "Interaction")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_byFraction.pdf")
for (i in 1:length(sigIR)){
x = ggplot(sigIR[[i]], aes(x=IRratio)) +
  geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlim(0,0.15) +
  xlab("IR Ratio") +
  ggtitle(paste0("Intron Retention: ", gr[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
print(x)
}
dev.off() 



### Get IR info for genes significantly regulated by Age

sigIR = lapply(age.sig, function(x) max[which(max$genes %in% x$ensID),])
unlist(lapply(sigIR, function(x) length(unique(x$genes))))

# genes >50% or >10% retained in at least one sample
perc = lapply(sigIR, function(x) list(sig.more50=as.character(unique(x[which(x$IRratio>=0.50),"genes"])),
                                      sig.less50=as.character(unique(x[which(x$IRratio<0.50),"genes"])),
                                      all.more50=as.character(unique(max[which(max$IRratio>=0.50),"genes"])),
                                      all.less50=as.character(unique(max[which(max$IRratio<0.50),"genes"])),
                                      sig.more10=as.character(unique(x[which(x$IRratio>=0.10),"genes"])),
                                      sig.less10=as.character(unique(x[which(x$IRratio<0.10),"genes"])),
                                      all.more10=as.character(unique(max[which(max$IRratio>=0.10),"genes"])),
                                      all.less10=as.character(unique(max[which(max$IRratio<0.10),"genes"]))))
for (i in 1:length(perc)){
  venn.diagram(perc[[i]][1:4], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",
                                      names(perc)[i], ".50.percent.IRratio.jpeg"), 
               main=names(perc)[i], col = "transparent", 
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"), alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold", 
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 50.percent.IRratio.overlap_byAge.pdf

for (i in 1:length(perc)){
  venn.diagram(perc[[i]][5:8], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",
                                      names(perc)[i], ".10.percent.IRratio.jpeg"), 
               main=names(perc)[i], col = "transparent", 
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"), alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold", 
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 10.percent.IRratio.overlap_byAge.pdf


# test whether genes with >50% or 10% intron retention in any sample are preferentially significantly DEG in the pattern listed

moreless = list()
for (i in 1:length(perc)){
  tmp = perc[[i]]
  moreless[[i]] = list("50%" = data.frame(greater = c(length(unique(tmp$sig.more50)),length(unique(tmp$all.more50[!(tmp$all.more50 %in% tmp$sig.more50)]))),
                                          less = c(length(unique(tmp$sig.less50)),length(unique(tmp$all.less50[!(tmp$all.less50 %in% tmp$sig.less50)]))), 
                                          row.names = c("sig","NS")),
                       "10%" = data.frame(greater = c(length(unique(tmp$sig.more10)),length(unique(tmp$all.more10[!(tmp$all.more10 %in% tmp$sig.more10)]))),
                                          less = c(length(unique(tmp$sig.less10)),length(unique(tmp$all.less10[!(tmp$all.less10 %in% tmp$sig.less10)]))), 
                                          row.names = c("sig","NS")))
}
names(moreless) = names(perc)
moreless.fisher = lapply(moreless, function(x) lapply(x,fisher.test))
write.csv(rbind(pval = data.frame(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/fisher_50perc_10perc_IR_sigDEGs_relationship_byAge.csv")


# Intron retention between retained and exported genes
IR_decr_incr = list(both = t.test(sigIR[["both_decreasing"]][,"IRratio"], sigIR[["both_increasing"]][,"IRratio"]),
                  adult = t.test(c(sigIR[["both_decreasing"]][,"IRratio"],sigIR[["Ad_decreasing"]][,"IRratio"]), 
                                 c(sigIR[["both_increasing"]][,"IRratio"],sigIR[["Ad_increasing"]][,"IRratio"])),
                  prenatal = t.test(c(sigIR[["both_decreasing"]][,"IRratio"],sigIR[["Fet_decreasing"]][,"IRratio"]), 
                                    c(sigIR[["both_increasing"]][,"IRratio"],sigIR[["Fet_increasing"]][,"IRratio"])))
write.csv(rbind(Tstat = data.frame(lapply(IR_decr_incr, function(x) round(x$statistic,3))), pval = data.frame(lapply(IR_decr_incr, function(x) x$p.value)),
                confInt = data.frame(lapply(IR_decr_incr, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(IR_decr_incr, function(x) round(x$estimate,3)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/Ttest_increasing_vs_decreasing_DEG_IRratios.csv")


# What does the distribution of IR Ratios by Sig group look like?
gr = c("Both Decreasing","Both Increasing","Decreasing in Cytosol","Decreasing in Nucleus","Increasing in Cytosol","Increasing in Nucleus",
       "Decreasing in Nucleus\nIncreasing in Cytosol","Decreasing in Cytosol\nIncreasing in Nucleus","Interaction")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_byAge.pdf")
for (i in 1:length(sigIR)){
  x = ggplot(sigIR[[i]], aes(x=IRratio)) +
    geom_density(aes(group=Group, colour=Group)) +
    ylab("") + 
    xlim(0,0.15) +
    xlab("IR Ratio") +
    ggtitle(paste0("Intron Retention: ", gr[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = c(0.8, 0.55)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(x)
}
dev.off()



# Plot LFC by Fraction and IR
FracByAge = list(Apres = data.frame(Apres), Fpres.down = data.frame(Fpres.down),Ipres.down = data.frame(Ipres.down))
FracByAge = Map(cbind, FracByAge,lapply(FracByAge, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Adult", "Prenatal", "Interaction"))
r = lapply(irbyGene, function(x) x[which(as.character(x$genes) %in% as.character(FracByAge[["Apres"]][,"ensemblID"])),])
FracByAgeIR = list(Map(cbind, r, lapply(r, function(x) FracByAge[["Apres"]][match(x$genes, FracByAge[["Apres"]][,"ensemblID"]),])),
                   Map(cbind, r, lapply(r, function(x) FracByAge[["Fpres.down"]][match(x$genes, FracByAge[["Fpres.down"]][,"ensemblID"]),])))
FracByAgeIR = do.call(rbind, lapply(FracByAgeIR, function(x) do.call(rbind, x)))
FracByAgeIR$SampleID = gsub("\\..*","", rownames(FracByAgeIR))
FracByAgeIR$rownum = c(1:nrow(FracByAgeIR))
FracByAgeIR$Fraction = ifelse((FracByAgeIR$rownum %in% grep("C", FracByAgeIR$SampleID)), "Cytosol", "Nucleus")
FracByAgeIR$Age = ifelse(FracByAgeIR$rownum %in% grep("53", FracByAgeIR$SampleID), "Prenatal", "Adult")
FracByAgeIR$Group = paste(FracByAgeIR$Age, FracByAgeIR$Fraction, sep=":")
FracByAgeIR$FDR = ifelse(FracByAgeIR$padj<=0.05, "FDR<0.05", "FDR>0.05")
FracByAgeIR$IR.50 = factor(ifelse(FracByAgeIR$IRratio>=0.5, ">0.5", "<0.5"))
FracByAgeIR$IR.10 = factor(ifelse(FracByAgeIR$IRratio>=0.1, ">0.1", "<0.1"))
FracByAgeIR = FracByAgeIR[which(FracByAgeIR$padj!="NA"),]
dim(FracByAgeIR[which(FracByAgeIR$genes!=FracByAgeIR$ensemblID),])

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/RNA_localization_byIRratio.pdf",width=7,height=5)
ggplot(FracByAgeIR, aes(x=IR.50, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(FracByAgeIR, aes(x=IR.10, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

# Plot LFC by Age and IR
AgebyFrac = list(Cpres.down = data.frame(Cpres.down), Npres = data.frame(Npres))
AgebyFrac = Map(cbind, AgebyFrac,lapply(AgebyFrac, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Cytosol", "Nucleus"))
r = lapply(irbyGene, function(x) x[which((as.character(x$genes) %in% as.character(AgebyFrac[["Cpres.down"]][,"ensemblID"])) &
                                           (as.character(x$genes) %in% as.character(AgebyFrac[["Npres"]][,"ensemblID"]))),])
AgebyFracIR = list(Map(cbind, r, lapply(r, function(x) AgebyFrac[["Cpres.down"]][match(x$genes, AgebyFrac[["Cpres.down"]][,"ensemblID"]),])),
                   Map(cbind, r, lapply(r, function(x) AgebyFrac[["Npres"]][match(x$genes, AgebyFrac[["Npres"]][,"ensemblID"]),])))
AgebyFracIR = do.call(rbind, lapply(AgebyFracIR, function(x) do.call(rbind, x)))
AgebyFracIR$SampleID = gsub("\\..*","", rownames(AgebyFracIR))
AgebyFracIR$rownum = c(1:nrow(AgebyFracIR))
AgebyFracIR$Fraction = ifelse((AgebyFracIR$rownum %in% grep("C", AgebyFracIR$SampleID)), "Cytosol", "Nucleus")
AgebyFracIR$Age = ifelse(AgebyFracIR$rownum %in% grep("53", AgebyFracIR$SampleID), "Prenatal", "Adult")
AgebyFracIR$Group = paste(AgebyFracIR$Age, AgebyFracIR$Fraction, sep=":")
AgebyFracIR$FDR = ifelse(AgebyFracIR$padj<=0.05, "FDR<0.05", "FDR>0.05")
AgebyFracIR$IR.50 = factor(ifelse(AgebyFracIR$IRratio>=0.5, ">0.5", "<0.5"))
AgebyFracIR$IR.10 = factor(ifelse(AgebyFracIR$IRratio>=0.1, ">0.1", "<0.1"))
AgebyFracIR = AgebyFracIR[which(AgebyFracIR$padj!="NA"),]
dim(AgebyFracIR[which(AgebyFracIR$genes!=AgebyFracIR$ensemblID),])

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/Devel_expression_trajectory_byIRratio.pdf",width=7,height=5)
ggplot(AgebyFracIR, aes(x=IR.50, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("Developmental Expression\nTrajectory by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # Save as Devel_expression_trajectory_byIRratio
ggplot(AgebyFracIR, aes(x=IR.10, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("Developmental Expression\nTrajectory by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # Save as Devel_expression_trajectory_byIRratio
dev.off()


fracDevel = lapply(sig, function(x) AgebyFracIR[which(AgebyFracIR$genes %in% x$ensID),])
gr = c("Both Retained", "Both Exported", "Retained in Prenatal", "Retained in Adult",
       "Exported in Prenatal", "Exported in Adult", "Retained in Adult/\nExported in Prenatal",
       "Retained in Prenatal/\nExported in Adult", "Interaction")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/RetainedbyAge_LFCxFDRxIR.pdf", width=8, height=6)
for (i in c(1:6,9)){
  g = ggplot(fracDevel[[i]], aes(x=IR.50, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("Age Expression Changes by IR Ratio:\n",gr[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
  g = ggplot(fracDevel[[i]], aes(x=IR.10, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("Age Expression Changes by IR Ratio:\n",gr[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
}
dev.off()

### Is there a relationship between direction of expression of significantly DE genes by age and IR being > or < 0.50 or 0.10?
# In cytosol:
unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])))))),function(x) x$p.value))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#0.8136241      1.0000000      0.3500000      0.3243825      1.0000000      0.8094819      1.0000000      1.0000000      0.1006644 

unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])))))),function(x) x$p.value))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#0.87326397     1.00000000     1.00000000     0.94828157     1.00000000     0.59759582     1.00000000     1.00000000     0.02381284

# In nucleus
unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])))))),function(x) x$p.value))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#0.6112488      1.0000000      0.4736842      0.1136210      1.0000000      1.0000000      1.0000000      1.0000000      1.0000000

unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])))))),function(x) x$p.value))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#0.8693460      1.0000000      1.0000000      0.9463244      1.0000000      0.5775031      1.0000000      1.0000000      0.4395440