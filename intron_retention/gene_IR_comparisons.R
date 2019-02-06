library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

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
max = do.call(rbind, Map(cbind, irbyGene, SampleID = as.list(names(irbyGene))))
max$rownum = c(1:nrow(max))
max$Fraction = ifelse(max$rownum %in% grep("C", max$SampleID), "Cytoplasm", "Nucleus")
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
x = data.frame(comp = rep.int(names(moreless),2), percent = c(rep.int("50%", length(moreless)),rep.int("10%", length(moreless))),
               pval = c(unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)[1])),
                        unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)[2]))),
               OddsRatio = c(unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)[1])),
                             unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)[2]))), row.names = NULL)
x$FDR = p.adjust(x$pval, method="fdr")
write.csv(x, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/fisher_50perc_10perc_IR_sigDEGs_relationship_byFraction.csv")


# Intron retention between retained and exported genes
IR_ret_exp = list(both = t.test(sigIR[["both_retained"]][,"IRratio"], sigIR[["both_exported"]][,"IRratio"]),
                  adult = t.test(c(sigIR[["both_retained"]][,"IRratio"],sigIR[["Ad_retained"]][,"IRratio"]), 
                                 c(sigIR[["both_exported"]][,"IRratio"],sigIR[["Ad_exported"]][,"IRratio"])),
                  prenatal = t.test(c(sigIR[["both_retained"]][,"IRratio"],sigIR[["Fet_retained"]][,"IRratio"]), 
                                    c(sigIR[["both_exported"]][,"IRratio"],sigIR[["Fet_exported"]][,"IRratio"])))
x = data.frame(Comp = names(IR_ret_exp), Tstat = unlist(lapply(IR_ret_exp, function(x) round(x$statistic,3))), pval = unlist(lapply(IR_ret_exp, function(x) x$p.value)),
               confInt1 = unlist(lapply(IR_ret_exp, function(x) round(x$conf.int,3)[1])), confInt2 = unlist(lapply(IR_ret_exp, function(x) round(x$conf.int,3)[2])),
               estMeans1 = unlist(lapply(IR_ret_exp, function(x) round(x$estimate,3)[1])), estMeans2 = unlist(lapply(IR_ret_exp, function(x) round(x$estimate,3)[2])),
               row.names = NULL)
x$FDR = p.adjust(x$pval, method="fdr")
write.csv(x, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/Ttest_exported_vs_retained_DEG_IRratios.csv")


# What does the distribution of IR Ratios by Sig group look like?
gr = c("Both Nuclear", "Both Cytoplasmic", "Nuclear in Prenatal", "Nuclear in Adult",
       "Cytoplasmic in Prenatal", "Cytoplasmic in Adult", "Nuclear in Adult/\nCytoplasmic in Prenatal",
       "Nuclear in Prenatal/\nCytoplasmic in Adult", "Interaction")
pdf("./Dropbox/sorted_figures/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_byFraction.pdf", width=5,height=5)
for (i in 1:length(sigIR)){
sigIR[[i]]$Group = factor(sigIR[[i]]$Group, levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm","Adult:Nucleus","Prenatal:Nucleus"))
x = ggplot(sigIR[[i]], aes(x=IRratio)) +
  geom_density(aes(group=Group, colour=Group), size=1) +
  ylab("") + 
  xlim(0,0.15) +
  xlab("IR Ratio") +
  ggtitle(paste0("Intron Retention: ", gr[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.75, 0.65)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
print(x)
}
x = rbind(cbind(sigIR[["both_retained"]], Genes = "Both Nuclear"), cbind(sigIR[["both_exported"]], Genes = "Both Cytoplasmic"))
x$Genes = factor(x$Genes, levels = c("Both Cytoplasmic","Both Nuclear"))
ggplot(x, aes(x=IRratio, linetype = Genes, col = Group)) + labs(fill="") +
  geom_density(size=1) + ylab("Density") + theme_classic() + xlim(0,0.15) + xlab("IR Ratio") +
  theme(title = element_text(size = 28)) +
  theme(text = element_text(size = 28)) +
  theme(legend.position = c(0.6, 0.65)) +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


### Get IR info for genes significantly regulated by Age

sigIR = lapply(age.sig, function(x) max[which(max$genes %in% x$ensID),])
unlist(lapply(sigIR, function(x) length(unique(x$genes))))


# genes >50% or >10% Nuclear in at least one sample
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
x = data.frame(comp = rep.int(names(moreless),2), percent = c(rep.int("50%", length(moreless)),rep.int("10%", length(moreless))),
               pval = c(unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)[1])),
                        unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)[2]))),
               OddsRatio = c(unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)[1])),
                             unlist(lapply(moreless.fisher, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)[2]))), row.names = NULL)
x$FDR = p.adjust(x$pval, method="fdr")
write.csv(x, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/fisher_50perc_10perc_IR_sigDEGs_relationship_byAge.csv")


# Intron retention between Nuclear and Cytoplasmic genes
IR_decr_incr = list(Both = t.test(sigIR[["both_decreasing"]][,"IRratio"], sigIR[["both_increasing"]][,"IRratio"]),
                    Cytoplasm = t.test(c(sigIR[["both_decreasing"]][,"IRratio"],sigIR[["Cyt_decreasing"]][,"IRratio"]), 
                                   c(sigIR[["both_increasing"]][,"IRratio"],sigIR[["Cyt_increasing"]][,"IRratio"])),
                    Nucleus = t.test(c(sigIR[["both_decreasing"]][,"IRratio"],sigIR[["Nuc_decreasing"]][,"IRratio"]), 
                                      c(sigIR[["both_increasing"]][,"IRratio"],sigIR[["Nuc_increasing"]][,"IRratio"])))
x = data.frame(Comp = names(IR_decr_incr), Tstat = unlist(lapply(IR_decr_incr, function(x) round(x$statistic,3))), pval = unlist(lapply(IR_decr_incr, function(x) x$p.value)),
               confInt1 = unlist(lapply(IR_decr_incr, function(x) round(x$conf.int,3)[1])), confInt2 = unlist(lapply(IR_decr_incr, function(x) round(x$conf.int,3)[2])),
               estMeans1 = unlist(lapply(IR_decr_incr, function(x) round(x$estimate,3)[1])), estMeans2 = unlist(lapply(IR_decr_incr, function(x) round(x$estimate,3)[2])),
               row.names = NULL)
x$FDR = p.adjust(x$pval, method="fdr")
write.csv(x, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/gene_IR_comparisons/Ttest_increasing_vs_decreasing_DEG_IRratios.csv")


# What does the distribution of IR Ratios by Sig group look like?
gr = c("Both Decreasing","Both Increasing","Decreasing in Cytoplasm","Decreasing in Nucleus","Increasing in Cytoplasm","Increasing in Nucleus",
       "Decreasing in Nucleus\nIncreasing in Cytoplasm","Decreasing in Cytoplasm\nIncreasing in Nucleus","Interaction")
pdf("./Dropbox/sorted_figures/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_byAge.pdf")
for (i in 1:length(sigIR)){
  sigIR[[i]]$Group = factor(sigIR[[i]]$Group, levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm","Adult:Nucleus","Prenatal:Nucleus"))
  x = ggplot(sigIR[[i]], aes(x=IRratio)) + geom_density(aes(group=Group, colour=Group)) +
    ylab("") + xlim(0,0.15) + xlab("IR Ratio") +
    ggtitle(paste0("Intron Retention: ", gr[i])) +
    theme(title = element_text(size = 20), text = element_text(size = 20)) +
    theme(legend.position = c(0.8, 0.55)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(x)
}
dev.off()
x = do.call(rbind, Map(cbind, sigIR[!names(sigIR) %in% c("decr_Nuc_incr_Cyt","decr_Cyt_incr_Nuc", "interacting")],
                       comp = as.list(c("Both Decreasing","Both Increasing","Decreasing in Cytoplasm",
                                 "Decreasing in Nucleus","Increasing in Cytoplasm","Increasing in Nucleus"))))
x$compLoc = Dir = NA
x[grep("Decreasing", x$comp),"Dir"] = "Decreasing"
x[grep("Increasing", x$comp),"Dir"] = "Increasing"
x[grep("Both", x$comp),"compLoc"] = "In Both"
x[grep("Cytoplasm", x$comp),"compLoc"] = "In Cytoplasm"
x[grep("Nucleus", x$comp),"compLoc"] = "In Nucleus"
x$Group = factor(x$Group, levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm","Adult:Nucleus","Prenatal:Nucleus"))

pdf("./Dropbox/sorted_figures/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_byAge2.pdf")
ggplot(x[which(x$compLoc!="In Both"),], aes(x=IRratio, linetype = Dir, col = Group)) +
  geom_density(size=1) + facet_grid(compLoc~.) + theme_classic() +
  ylab("Density") + xlim(0,0.15) + xlab("IR Ratio") +
  ggtitle("Developmental IR") +
  theme(title = element_text(size = 28), text = element_text(size = 28)) +
  theme(legend.position = c(0.6, 0.8)) + labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(x[which(x$compLoc!="In Both"),], aes(x=Dir,y=IRratio), fill = Group) +
  geom_boxplot() + facet_grid(compLoc~.) + theme_classic() +
  ylab("Density") + xlab("IR Ratio") +
  ggtitle("Developmental IR") +
  theme(title = element_text(size = 28), text = element_text(size = 28)) +
  theme(legend.position = c(0.6, 0.8)) + labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()



# Plot LFC by Fraction and IR
FracByAge = list(Apres = data.frame(Apres), Fpres.down = data.frame(Fpres.down),Ipres.down = data.frame(Ipres.down))
FracByAge = Map(cbind, FracByAge,lapply(FracByAge, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Adult", "Prenatal", "Interaction"))
r = lapply(irbyGene, function(x) x[which(as.character(x$genes) %in% c(as.character(FracByAge[["Apres"]][,"ensemblID"]),as.character(FracByAge[["Fpres.down"]][,"ensemblID"]))),])
FracByAgeIR = list(Map(cbind, r, lapply(r, function(x) FracByAge[["Apres"]][match(x$genes, FracByAge[["Apres"]][,"ensemblID"]),])),
                   Map(cbind, r, lapply(r, function(x) FracByAge[["Fpres.down"]][match(x$genes, FracByAge[["Fpres.down"]][,"ensemblID"]),])))
FracByAgeIR = do.call(rbind, lapply(FracByAgeIR, function(x) do.call(rbind, x)))
FracByAgeIR$SampleID = gsub("\\..*","", rownames(FracByAgeIR))
FracByAgeIR$rownum = c(1:nrow(FracByAgeIR))
head(FracByAgeIR)
FracByAgeIR$Fraction = ifelse((FracByAgeIR$rownum %in% grep("C", FracByAgeIR$SampleID)), "Cytoplasm", "Nucleus")
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
  scale_fill_manual(values=c("red3","gray47")) +
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
  scale_fill_manual(values=c("red3","gray47")) +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(FracByAgeIR, aes(x=IR.50, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  scale_fill_manual(values=c("red3","gray47")) +
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
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  scale_fill_manual(values=c("red3","gray47")) +
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
                Comparison = list("Cytoplasm", "Nucleus"))
r = lapply(irbyGene, function(x) x[which((as.character(x$genes) %in% as.character(AgebyFrac[["Cpres.down"]][,"ensemblID"])) &
                                           (as.character(x$genes) %in% as.character(AgebyFrac[["Npres"]][,"ensemblID"]))),])
AgebyFracIR = list(Map(cbind, r, lapply(r, function(x) AgebyFrac[["Cpres.down"]][match(x$genes, AgebyFrac[["Cpres.down"]][,"ensemblID"]),])),
                   Map(cbind, r, lapply(r, function(x) AgebyFrac[["Npres"]][match(x$genes, AgebyFrac[["Npres"]][,"ensemblID"]),])))
AgebyFracIR = do.call(rbind, lapply(AgebyFracIR, function(x) do.call(rbind, x)))
AgebyFracIR$SampleID = gsub("\\..*","", rownames(AgebyFracIR))
AgebyFracIR$rownum = c(1:nrow(AgebyFracIR))
AgebyFracIR$Fraction = ifelse((AgebyFracIR$rownum %in% grep("C", AgebyFracIR$SampleID)), "Cytoplasm", "Nucleus")
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
  scale_fill_manual(values=c("red3","gray47")) +
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
  scale_fill_manual(values=c("red3","gray47")) +
  ggtitle("Developmental Expression\nTrajectory by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # Save as Devel_expression_trajectory_byIRratio
ggplot(AgebyFracIR, aes(x=IR.50, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  facet_grid(. ~ Comparison) +
  scale_fill_manual(values=c("red3","gray47")) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("Developmental Expression\nTrajectory by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # Save as Devel_expression_trajectory_byIRratio
ggplot(AgebyFracIR, aes(x=IR.10, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  scale_fill_manual(values=c("red3","gray47")) +
  ggtitle("Developmental Expression\nTrajectory by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # Save as Devel_expression_trajectory_byIRratio
dev.off()


fracDevel = lapply(sig, function(x) AgebyFracIR[which(AgebyFracIR$genes %in% x$ensID),])
gr = c("Both Nuclear", "Both Cytoplasmic", "Nuclear in Prenatal", "Nuclear in Adult",
       "Cytoplasmic in Prenatal", "Cytoplasmic in Adult", "Nuclear in Adult/\nCytoplasmic in Prenatal",
       "Nuclear in Prenatal/\nCytoplasmic in Adult", "Interaction")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/RetainedbyAge_LFCxFDRxIR.pdf", width=8, height=6)
for (i in c(1:6,9)){
  g = ggplot(fracDevel[[i]], aes(x=IR.50, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    scale_fill_manual(values=c("red3","gray47")) +
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
    scale_fill_manual(values=c("red3","gray47")) +
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
# In Cytoplasm:
unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])))))),function(x) x$p.value))
# both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#   0.373284578    1.000000000    0.578392622    0.152100669    1.000000000    0.533118379    1.000000000    1.000000000    0.009861048 

unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Cytoplasm" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])))))),function(x) x$p.value))
# both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#   0.283520158    0.719976209    0.736240267    1.000000000    1.000000000    0.676810029    1.000000000    1.000000000    0.001621381 

# In nucleus
unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50=="<0.5"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                           x$padj<=0.05 & x$IR.50==">0.5"),"genes"])))))),function(x) x$p.value))
# both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#    0.57336751     0.68251619     0.09834369     0.27529318     1.00000000     0.47065787     1.00000000     1.00000000     0.56873148 

unlist(lapply(lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10=="<0.1"),"genes"]))),
                                                                   c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])),
                                                                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                                             x$padj<=0.05 & x$IR.10==">0.1"),"genes"])))))),function(x) x$p.value))
# both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#     0.3432458      0.7384542      0.7379863      0.8186349      1.0000000      0.8508929      1.0000000      1.0000000      0.6351638 
