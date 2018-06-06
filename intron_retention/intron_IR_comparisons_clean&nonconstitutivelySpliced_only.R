library(GenomicRanges)
library(data.table)
library(ggplot2)
library(VennDiagram)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

### Differences by intron
## Differential IR by group
# read in results files

comps = c("Adult_PolyA_Zone_cleanIntrons_adultShared","Fetal_PolyA_Zone_cleanIntrons_prenatalShared",
          "Cytosol_PolyA_Age_cleanIntrons_cytosolShared","Nuclear_PolyA_Age_cleanIntrons_nucleusShared")
IRcomp = list()
for (i in 1:length(comps)){
  IRcomp[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/", comps[i], ".tab"), header = TRUE, comment.char="#")
}
names(IRcomp) = c("Adult_byFraction","Fetal_byFraction","Cytosol_byAge","Nuclear_byAge")
IRcomp = Map(cbind, IRcomp,
             intronID = lapply(IRcomp, function(x) paste0(x$Intron.GeneName.GeneID,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)),
             ensID = lapply(lapply(IRcomp, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE)), function(y) y[grep("ENSG", y)]),
             IR.diff = lapply(IRcomp, function(y) y$A.IRratio - y$B.IRratio),
             Sign = lapply(IRcomp, function(y) ifelse((y$A.IRratio - y$B.IRratio) < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")),
             NonUniformIntronCover = lapply(IRcomp, function(x) paste(ifelse(x$A.warnings=="NonUniformIntronCover", "Yes","No"), 
                                                                      ifelse(x$B.warnings=="NonUniformIntronCover", "Yes","No"), sep=":")),
             LowCover = lapply(IRcomp, function(x) paste(ifelse(x$A.warnings=="LowCover", "Yes","No"), 
                                                         ifelse(x$B.warnings=="LowCover", "Yes","No"), sep=":")),
             LowSplicing = lapply(IRcomp, function(x) paste(ifelse(x$A.warnings=="LowSplicing", "Yes","No"), 
                                                            ifelse(x$B.warnings=="LowSplicing", "Yes","No"), sep=":")))
full = list(adult = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-adultShared.txt", header = TRUE),
            prenatal = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt", header = TRUE),
            cytosol = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt", header = TRUE),
            nucleus = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt", header = TRUE))
total = as.list(elementNROWS(full))
names(total) = names(IRcomp)
IRcomp = Map(cbind, IRcomp, padj = mapply(function(p,t) p.adjust(p, method = "fdr", n = t), lapply(IRcomp, function(x) x$p.diff), total))


do.call(rbind, lapply(IRcomp, function(x) table(x$NonUniformIntronCover)))
#                 No:No No:Yes Yes:No
#Adult_byFraction   938    800    305
#Cytosol_byAge      693    411    228
#Fetal_byFraction   693    537    352
#Nuclear_byAge     1133    511    474
do.call(rbind, lapply(IRcomp, function(x) table(x$LowCover)))
#                 No:No No:Yes Yes:No
#Adult_byFraction  1897     28    118
#Cytosol_byAge     1102     43    187
#Fetal_byFraction  1538      8     36
#Nuclear_byAge     1735    147    236
do.call(rbind, lapply(IRcomp, function(x) table(x$LowSplicing)))
#                 No:No No:Yes Yes:No
#Adult_byFraction  2043      0      0
#Cytosol_byAge     1330      2   1330
#Fetal_byFraction  1582      0      0
#Nuclear_byAge     2108      7      3

elementNROWS(lapply(IRcomp, function(y) y[which(y$NonUniformIntronCover=="No:No" & y$p.diff<=0.05),]))
# Adult_byFraction Fetal_byFraction    Cytosol_byAge    Nuclear_byAge 
#              167               99              101              241 


## Add introns that were included in test but not reported

for (i in 1:length(full)) { colnames(full[[i]])[1:7] = c("Chr","Start","End","Intron.GeneName.GeneID","X.","Direction","ExcludedBases") }
full = lapply(full, function(x) data.frame(x[,1:7], "p.diff"=NA,"p.increased"=NA,"p.decreased"=NA,"A.IRratio"=NA,"A.warnings"=NA,"A.IntronCover"=NA,"A.IntronDepth"=NA,"A.SplicesMax"=NA,
                                           "A.SplicesExact"=NA,"B.IRratio"=NA,"B.warnings"=NA,"B.IntronCover"=NA,"B.IntronDepth"=NA,"B.SplicesMax"=NA,"B.SplicesExact"=NA,"replicates"=NA,
                                           "A1.IRratio"=NA,"A2.IRratio"=NA,"A3.IRratio"=NA,"B1.IRratio"=NA,"B2.IRratio"=NA,"B3.IRratio"=NA))
full = Map(cbind, full, intronID = lapply(full, function(x) paste0(x$Intron.GeneName.GeneID,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)),
             ensID = lapply(lapply(full, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE)), function(y) y[grep("ENSG", y)]),
             IR.diff = NA, Sign = NA, NonUniformIntronCover = NA, LowCover = NA, LowSplicing = NA, padj=1)
IRcomp$Adult_byFraction = rbind(IRcomp$Adult_byFraction, full$adult[-which(full$adult$intronID %in% IRcomp$Adult_byFraction$intronID),])
IRcomp$Fetal_byFraction = rbind(IRcomp$Fetal_byFraction, full$prenatal[-which(full$prenatal$intronID %in% IRcomp$Fetal_byFraction$intronID),])
IRcomp$Cytosol_byAge = rbind(IRcomp$Cytosol_byAge, full$cytosol[-which(full$cytosol$intronID %in% IRcomp$Cytosol_byAge$intronID),])
IRcomp$Nuclear_byAge = rbind(IRcomp$Nuclear_byAge, full$nucleus[-which(full$nucleus$intronID %in% IRcomp$Nuclear_byAge$intronID),])


## Explore the results

counts = data.frame(rbind(unlist(total), elementNROWS(IRcomp),
                          elementNROWS(lapply(IRcomp, function(y) y[which(y$padj<=0.05),])),
                          elementNROWS(lapply(IRcomp, function(y) y[which(y$Sign=="MoreIRInNuc.Fetal"),])),
                          elementNROWS(lapply(IRcomp, function(y) y[which(y$padj<=0.05 & y$Sign=="MoreIRInNuc.Fetal"),])),
                          elementNROWS(lapply(IRcomp, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),])),
                          elementNROWS(lapply(IRcomp, function(y) y[which((y$A.IRratio>=0.5 | y$B.IRratio>=0.5) & y$Sign=="MoreIRInNuc.Fetal"),]))),
                    row.names = c("total measured","all reported", "FDR < 0.05", "MoreIRInNuc.Fetal", "FDR < 0.05 and MoreIRInNuc.Fetal",
                                  "IRratio > 0.5","IRratio > 0.5 and MoreIRInNuc.Fetal"))
write.csv(counts, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/numIntrons_clean_nonconstitutivelySpliced.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/numIntrons_clean_nonconstitutivelySpliced.csv")
df[df$X=="FDR < 0.05 and MoreIRInNuc.Fetal",2:5]/df[df$X=="FDR < 0.05",2:5]
#Adult_byFraction Fetal_byFraction Cytosol_byAge Nuclear_byAge
#               1                1           0.9     0.9047619


## Comparison of significantly vs nonsignificantly retained introns by zone/age

tables = list()
for (i in 1:4){
  tables[[i]] = data.frame(MoreIRInNuc.Fetal = c(counts["FDR < 0.05 and MoreIRInNuc.Fetal",i], counts["MoreIRInNuc.Fetal",i]-counts["FDR < 0.05 and MoreIRInNuc.Fetal",i]),
                           MoreIRINCyt.Adult = c(counts["FDR < 0.05",i]-counts["FDR < 0.05 and MoreIRInNuc.Fetal",i],
                                                 counts["all reported",i]-counts["FDR < 0.05",i]-(counts["MoreIRInNuc.Fetal",i]-counts["FDR < 0.05 and MoreIRInNuc.Fetal",i])),
                           row.names = c("Sig","nonSig"))
}
names(tables) = colnames(counts)

df = data.frame(p.value = unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value)), OR = unlist(lapply(lapply(tables, fisher.test), function(x) x$estimate)))
df$FDR = p.adjust(df$p.value, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonSig_IR_byIntron_byFraction_byAge.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonSig_IR_byIntron_byFraction_byAge.csv")
df


# Get the DEG pval and LFC sign by Fraction for differentially retained introns

degs = list(Apres = data.frame(Apres, Sign = ifelse(Apres$log2FoldChange>0, "Pos","Neg")), Fpres = data.frame(Fpres.down, Sign = ifelse(Fpres.down$log2FoldChange>0, "Pos","Neg")), 
            Cpres = data.frame(Cpres.down), Npres = data.frame(Npres))
degs = Map(cbind, degs, lapply(degs, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
IRcomp = Map(cbind, IRcomp, 
             AP.sig = lapply(IRcomp, function(x) degs$Apres[match(x$ensID, degs$Apres$ensemblID),"padj"]),
             AP.LFC = lapply(IRcomp, function(x) degs$Apres[match(x$ensID, degs$Apres$ensemblID),"log2FoldChange"]),
             FP.sig = lapply(IRcomp, function(x) degs$Fpres[match(x$ensID, degs$Fpres$ensemblID),"padj"]),
             FP.LFC = lapply(IRcomp, function(x) degs$Fpres[match(x$ensID, degs$Fpres$ensemblID),"log2FoldChange"]),
             CP.sig = lapply(IRcomp, function(x) degs$Cpres[match(x$ensID, degs$Cpres$ensemblID),"padj"]),
             CP.LFC = lapply(IRcomp, function(x) degs$Cpres[match(x$ensID, degs$Cpres$ensemblID),"log2FoldChange"]),
             NP.sig = lapply(IRcomp, function(x) degs$Npres[match(x$ensID, degs$Npres$ensemblID),"padj"]),
             NP.LFC = lapply(IRcomp, function(x) degs$Npres[match(x$ensID, degs$Npres$ensemblID),"log2FoldChange"]))
lapply(IRcomp, head)



# Are genes with significantly differentially retained introns more likely to be significantly DEG by fraction?

P = c("AP", "FP", "CP", "NP")
tables = list(list(), list(), list(), list())
for (i in 1:length(IRcomp)){
  for (j in 1:length(P)) {
    tables[[i]][[j]] = data.frame(sigIR = c(nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]<=0.05 & IRcomp[[i]][,paste0(P[j],".sig")]<=0.05),]),
                                            nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]<=0.05 & IRcomp[[i]][,paste0(P[j],".sig")]>0.05),])),
                                  nonsigIR = c(nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]>0.05 & IRcomp[[i]][,paste0(P[j],".sig")]<=0.05),]),
                                               nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]>0.05 & IRcomp[[i]][,paste0(P[j],".sig")]>0.05),])), row.names = c("sigDEG","nonsigDEG"))
  }
  names(tables[[i]]) = paste0(P, "-DEG")
}
names(tables) = names(IRcomp)
fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = data.frame(pval = unlist(lapply(fisher, function(x) lapply(x, function(y) y$p.value))), OR = unlist(lapply(fisher, function(x) lapply(x, function(y) y$estimate))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F,file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonSig_IR_byIntron_and_sig_vs_nonSig_DEG_byFraction_byAge.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonSig_IR_byIntron_and_sig_vs_nonSig_DEG_byFraction_byAge.csv")
df[df$FDR<=0.05,]


# Are genes with significantly differentially retained introns more likely to have LFC in one direction by fraction?

for (i in 1:length(IRcomp)){
  for (j in 1:length(P)) {
    tables[[i]][[j]] = data.frame(sigIR = c(nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]<=0.05 & IRcomp[[i]][,paste0(P[j],".LFC")]>0),]),
                                            nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]<=0.05 & IRcomp[[i]][,paste0(P[j],".LFC")]<0),])),
                                  nonsigIR = c(nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]>0.05 & IRcomp[[i]][,paste0(P[j],".LFC")]>0),]),
                                               nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]>0.05 & IRcomp[[i]][,paste0(P[j],".LFC")]<0),])), row.names = c("posLFC","negLFC"))
  }
  names(tables[[i]]) = paste0(P, "-DEG")
}
names(tables) = names(IRcomp)
fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = data.frame(pval = unlist(lapply(fisher, function(x) lapply(x, function(y) y$p.value))), OR = unlist(lapply(fisher, function(x) lapply(x, function(y) y$estimate))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F,file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonSig_IR_byIntron_and_pos_vs_neg_geneLFC_byFraction_byAge.csv")
df[df$FDR<=0.05,]


# Are significantly differentially retained introns by fraction more likely to have a higher IR ratio in nuclear RNA?

for (i in 1:length(IRcomp)){
    tables[[i]] = data.frame(sigIR = c(nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]<=0.05 & IRcomp[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),]),
                                            nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]<=0.05 & IRcomp[[i]][,"Sign"]=="MoreIRInCyt.Adult"),])),
                                  nonsigIR = c(nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]>0.05 & IRcomp[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),]),
                                               nrow(IRcomp[[i]][which(IRcomp[[i]][,"padj"]>0.05 & IRcomp[[i]][,"Sign"]=="MoreIRInCyt.Adult"),])), row.names = c("MoreIRInNuc.Fetal","MoreIRInCyt.Adult"))
}
names(tables) = names(IRcomp)
fisher = lapply(tables, fisher.test)
df = data.frame(pval = unlist(lapply(fisher, function(y) y$p.value)), OR = unlist(lapply(fisher, function(y) y$estimate)))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F,file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonSig_IR_byIntron_and_MoreIRInNuc.Fetal_byFraction_byAge.csv")
df


# IR regulated genes 

IRlist = unlist(lapply(IRcomp, function(x) split(x, x$Sign)), recursive = F)
names(IRlist) = names = c("Adult:Cytoplasm-Enriched","Adult:Nuclear-Enriched","Prenatal:Cytoplasm-Enriched","Prenatal:Nuclear-Enriched",
                          "Cytoplasm:Adult-Enriched","Cytoplasm:Prenatal-Enriched","Nucleus:Adult-Enriched","Nucleus:Prenatal-Enriched")
degs = Map(cbind, degs, Comparison = list("Adult","Prenatal","Cytoplasm","Nucleus"), Collapsed.Comparison = list("Fraction","Fraction","Age","Age"), 
           FDR = lapply(degs, function(x) ifelse(x$padj<=0.05, "FDR<0.05", "FDR>0.05")))
degs = do.call(rbind, degs)
degs = degs[which(degs$padj!="NA"),]
IRlist = lapply(IRlist, function(x) x[which(x$padj <= 0.05),])
IRlist = lapply(IRlist, function(x) degs[which(degs$ensemblID %in% x$ensID),])
elementNROWS(IRlist)
IRlist = do.call(rbind, Map(cbind, IRlist[elementNROWS(IRlist)>0], IRgroup = as.list(names(IRlist[elementNROWS(IRlist)>0]))))
IRlist$IRgroup = gsub(":",":\n", IRlist$IRgroup, fixed = T)
IRlist$IRgroup = factor(IRlist$IRgroup, levels=c("Adult:\nNuclear-Enriched","Prenatal:\nNuclear-Enriched","Cytoplasm:\nAdult-Enriched",
                                                 "Nucleus:\nAdult-Enriched","Cytoplasm:\nPrenatal-Enriched","Nucleus:\nPrenatal-Enriched"))

# Plot DEG by fraction LFC and FDR for genes that contain introns that are differentially retained

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRgenes_byFraction.pdf", width=14, height=5)
ggplot(IRlist[which(IRlist$Collapsed.Comparison=="Fraction"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() + facet_grid(. ~ IRgroup) +
  scale_fill_manual(values=c("red3","gray47")) +
  ylab("Log2 Fold Change") + 
  xlab("") + geom_hline(yintercept=0, linetype="dotted") +
  ggtitle("Fraction Expression Changes in Retained Introns") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.position = "bottom",
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(IRlist[which(IRlist$Collapsed.Comparison=="Fraction"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  scale_fill_manual(values=c("red3","gray47")) +
  ylab("Log2 Fold Change") + 
  xlab("") + facet_grid(. ~ IRgroup) +
  ggtitle("Fraction Expression Changes in Retained Introns") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.position = "bottom",
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


# Plot DEG by age LFC and FDR for genes that contain introns that are differentially retained

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRgenes_byAge.pdf", width=14, height=5)
ggplot(IRlist[which(IRlist$Collapsed.Comparison=="Age"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() + geom_hline(yintercept=0, linetype="dotted") +
    scale_fill_manual(values=c("red3","gray47")) +
    ylab("Log2 Fold Change") + 
    xlab("") + facet_grid(. ~ IRgroup) +
    ggtitle("Age Expression Changes in Retained Introns") + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"), legend.position = "bottom",
          legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(IRlist[which(IRlist$Collapsed.Comparison=="Age"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
    scale_fill_manual(values=c("red3","gray47")) +
    xlab("") + facet_grid(. ~ IRgroup) +
    ggtitle("Age Expression Changes in Retained Introns") + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"), legend.position = "bottom",
          legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRgenes_byFrac_forPaper.pdf", width=11, height=5)
ggplot(IRlist[which(IRlist$Collapsed.Comparison=="Fraction" & IRlist$IRgroup %in% c("Cytoplasm:\nAdult-Enriched",
"Cytoplasm:\nPrenatal-Enriched","Nucleus:\nAdult-Enriched","Nucleus:\nPrenatal-Enriched")),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  scale_fill_manual(values=c("red3","gray47")) +
  ylab("Log2 Fold Change") + 
  xlab("") + facet_grid(. ~ IRgroup) +
  ggtitle("Fraction Expression Changes in dIRs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRgenes_byAge_forPaper.pdf", width=7, height=5)
ggplot(IRlist[which(IRlist$Collapsed.Comparison=="Age" & IRlist$IRgroup %in% c("Adult:\nNuclear-Enriched","Prenatal:\nNuclear-Enriched")),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dotted") +
  scale_fill_manual(values=c("red3","gray47")) +
  xlab("") + facet_grid(. ~ IRgroup) +
  ggtitle("Age Expression Changes in dIRs") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


## Are the Fraction or Age LFC values greater in genes containing a retained intron differentially?

t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Adult"),"log2FoldChange"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Adult"),"log2FoldChange"])
t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Prenatal"),"log2FoldChange"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Prenatal"),"log2FoldChange"])
t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"])
t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Nucleus"),"log2FoldChange"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Nucleus"),"log2FoldChange"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Adult"),"log2FoldChange"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Adult"),"log2FoldChange"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Prenatal"),"log2FoldChange"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Prenatal"),"log2FoldChange"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Nucleus"),"log2FoldChange"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Nucleus"),"log2FoldChange"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Adult"),"log2FoldChange"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Adult"),"log2FoldChange"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Prenatal"),"log2FoldChange"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Prenatal"),"log2FoldChange"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Nucleus"),"log2FoldChange"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Nucleus"),"log2FoldChange"])
tables = list(t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Adult"),"log2FoldChange"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Adult"),"log2FoldChange"]),
              t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Prenatal"),"log2FoldChange"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Prenatal"),"log2FoldChange"]),
              t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Cytoplasm"),"log2FoldChange"]),
              t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Nucleus"),"log2FoldChange"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Nucleus"),"log2FoldChange"]))
names(tables) = paste0("Nucleus_byAge_", unique(degs$Comparison), "_LFC")
df = data.frame(Comp = names(tables), tstat = unlist(lapply(tables, function(x) x$statistic)),pval = unlist(lapply(tables, function(x) x$p.value)),
                ConfInt1 = unlist(lapply(tables, function(x) x$conf.int[1])), ConfInt2 = unlist(lapply(tables, function(x) x$conf.int[2])),
                Mean1 = unlist(lapply(tables, function(x) x$estimate[1])), Mean2 = unlist(lapply(tables, function(x) x$estimate[2])), row.names = NULL)
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ttest_nuclearIntrons_byAge_geneLFC_diff.csv")
df


## Are the Fraction or Age p.values lower in genes containing a retained intron differentially?

t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Adult"),"padj"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Adult"),"padj"])
t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Prenatal"),"padj"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Prenatal"),"padj"])
t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Cytoplasm"),"padj"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Cytoplasm"),"padj"])
t.test(x = IRlist$"Adult:Cytoplasm-Enriched"[which(IRlist$"Adult:Cytoplasm-Enriched"$Comparison=="Nucleus"),"padj"],
       y = IRlist$"Adult:Nuclear-Enriched"[which(IRlist$"Adult:Nuclear-Enriched"$Comparison=="Nucleus"),"padj"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Adult"),"padj"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Adult"),"padj"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Prenatal"),"padj"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Prenatal"),"padj"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Cytoplasm"),"padj"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Cytoplasm"),"padj"])
t.test(x = IRlist$"Prenatal:Cytoplasm-Enriched"[which(IRlist$"Prenatal:Cytoplasm-Enriched"$Comparison=="Nucleus"),"padj"],
       y = IRlist$"Prenatal:Nuclear-Enriched"[which(IRlist$"Prenatal:Nuclear-Enriched"$Comparison=="Nucleus"),"padj"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Adult"),"padj"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Adult"),"padj"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Prenatal"),"padj"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Prenatal"),"padj"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Cytoplasm"),"padj"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Cytoplasm"),"padj"])
t.test(x = IRlist$"Cytoplasm:Adult-Enriched"[which(IRlist$"Cytoplasm:Adult-Enriched"$Comparison=="Nucleus"),"padj"],
       y = IRlist$"Cytoplasm:Prenatal-Enriched"[which(IRlist$"Cytoplasm:Prenatal-Enriched"$Comparison=="Nucleus"),"padj"])
tables = list(t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Adult"),"padj"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Adult"),"padj"]),
              t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Prenatal"),"padj"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Prenatal"),"padj"]),
              t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Cytoplasm"),"padj"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Cytoplasm"),"padj"]),
              t.test(x = IRlist$"Nucleus:Adult-Enriched"[which(IRlist$"Nucleus:Adult-Enriched"$Comparison=="Nucleus"),"padj"],
                     y = IRlist$"Nucleus:Prenatal-Enriched"[which(IRlist$"Nucleus:Prenatal-Enriched"$Comparison=="Nucleus"),"padj"]))
names(tables) = paste0("Nucleus_byAge_", unique(degs$Comparison), "_padj")
df = data.frame(Comp = names(tables), tstat = unlist(lapply(tables, function(x) x$statistic)),pval = unlist(lapply(tables, function(x) x$p.value)),
                ConfInt1 = unlist(lapply(tables, function(x) x$conf.int[1])), ConfInt2 = unlist(lapply(tables, function(x) x$conf.int[2])),
                Mean1 = unlist(lapply(tables, function(x) x$estimate[1])), Mean2 = unlist(lapply(tables, function(x) x$estimate[2])), row.names = NULL)
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F,file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ttest_nuclearIntrons_byAge_padj_diff.csv")
df


# plot IR ratio in Cytoplasmic prenatal samples and LFC by age in Cytoplasm
ggplot(IRcomp$Fetal_byFraction[which(IRcomp$Fetal_byFraction$padj<0.05),], aes(x=IR.diff, y=CP.LFC)) + 
  geom_point() +
  ylab("Log2 Fold Change") + ylim(-2,2) + 
  xlab("Difference in IR Ratio") +
  ggtitle("Age LFC and IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))


# Are genes with greater IR in Cytoplasm in prenatal higher- or lower-expressed in prenatal than adult in Cytoplasm?

for (i in 1:length(IRcomp)){
tables[[i]] = list(inAdult = t.test(IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInCyt.Adult"),"AP.LFC"],
                                    IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"AP.LFC"]),
                   inPrenatal = t.test(IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInCyt.Adult"),"FP.LFC"],
                                       IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"FP.LFC"]),
                   inCytoplasm = t.test(IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInCyt.Adult"),"CP.LFC"],
                                      IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"CP.LFC"]),
                   inNucleus = t.test(IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInCyt.Adult"),"NP.LFC"],
                                      IRcomp[[i]][which(IRcomp[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"NP.LFC"]))
}
names(tables) = names(IRcomp)
df = data.frame(Comp = names(unlist(tables))[grep("alternative",names(unlist(tables)))], 
                tstat = unlist(lapply(tables, function(x) lapply(x, function(y) y$statistic))),
                pval = unlist(lapply(tables, function(x) lapply(x, function(y) y$p.value))),
                ConfInt1 = unlist(lapply(tables, function(x) lapply(x, function(y) y$conf.int[1]))), 
                ConfInt2 = unlist(lapply(tables, function(x) lapply(x, function(y) y$conf.int[2]))),
                Mean1 = unlist(lapply(tables, function(x) lapply(x, function(y) y$estimate[1]))), 
                Mean2 = unlist(lapply(tables, function(x) lapply(x, function(y) y$estimate[2]))), row.names = NULL)
df$Comp = gsub(".alternative", "", df$Comp)
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F,file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ttest_MoreIRInNuc.Fetal_IR_vs_geneLFC_inGroup.csv")
df[df$FDR<=.05,]


## Are dIR introns by age more likely to be dIR introns by fraction?

sigIR.intron = lapply(IRcomp, function(x) as.character(x[which(x$padj<=0.05),"intronID"])) 
path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_"
venn.diagram(sigIR.intron, paste0(path, "sig_byFrac_byAge", ".jpeg"), 
             main="Significantly differentially retains\nintrons by group", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", sub.cex=3,margin=0.2)


nonsigIR.intron = lapply(IRcomp, function(x) as.character(x[which(x$padj>0.05),"intronID"]))
sigIR.intron.combined = list("Significantly IR\nBy Fraction" = c(sigIR.intron$Adult_byFraction, sigIR.intron$Prenatal_byFraction),
                             "Significantly IR\nBy Age" = c(sigIR.intron$Cytosol_byAge, sigIR.intron$Nuclear_PolyA_Age))
nonsigIR.intron.combined = list("Non-Significantly IR\nBy Fraction" = c(nonsigIR.intron$Adult_byFraction, nonsigIR.intron$Prenatal_byFraction),
                                "Non-Significantly IR\nBy Age" = c(nonsigIR.intron$Cytosol_byAge, nonsigIR.intron$Nuclear_PolyA_Age))
comps = list(c(1,3), c(1,4), c(1:2), c(2,3), c(2,4), c(3,4))
ids = c("FractionbyAge_adult_cytosol", "FractionbyAge_adult_nucleus", "Fraction_adult_prenatal",
        "FractionbyAge_prenatal_cytosol","FractionbyAge_prenatal_nucleus","Age_cytosol_nucleus")

for (i in 1:length(comps)) {
venn.diagram(c(sigIR.intron[comps[[i]]], nonsigIR.intron[comps[[i]]]), paste0(path, ids[i], ".jpeg"), 
             main=paste0("dIR_", ids[i]), col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
}
venn.diagram(c(sigIR.intron.combined, nonsigIR.intron.combined), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_combined.jpeg", 
             main="dIR_FractionbyAge_combined", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
# all 7 saved together at ./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_overlap_Fraction_Age.pdf


tables = list()
for (i in 1:length(comps)) {
  tables[[i]] = data.frame(sig1 = c(length(sigIR.intron[[comps[[i]][1]]][sigIR.intron[[comps[[i]][1]]] %in% sigIR.intron[[comps[[i]][2]]]]),
                                    length(sigIR.intron[[comps[[i]][1]]][sigIR.intron[[comps[[i]][1]]] %in% nonsigIR.intron[[comps[[i]][2]]]])),
                      nonsig1 = c(length(nonsigIR.intron[[comps[[i]][1]]][nonsigIR.intron[[comps[[i]][1]]] %in% sigIR.intron[[comps[[i]][2]]]]),
                                  length(nonsigIR.intron[[comps[[i]][1]]][nonsigIR.intron[[comps[[i]][1]]] %in% nonsigIR.intron[[comps[[i]][2]]]])),row.names = c("sig2","nonsig2"))
}
names(tables) = ids
fisher = c(lapply(tables, fisher.test), list(combined_fraction_age = fisher.test(data.frame(
  sig1 = c(length(sigIR.intron.combined[[1]][sigIR.intron.combined[[1]] %in% sigIR.intron.combined[[2]]]),
           length(sigIR.intron.combined[[1]][sigIR.intron.combined[[1]] %in% nonsigIR.intron.combined[[2]]])),
  nonsig1 = c(length(nonsigIR.intron.combined[[1]][nonsigIR.intron.combined[[1]] %in% sigIR.intron.combined[[2]]]),
              length(nonsigIR.intron.combined[[1]][nonsigIR.intron.combined[[1]] %in% nonsigIR.intron.combined[[2]]])),
  row.names = c("sig2","nonsig2")))))
df = data.frame(pval = unlist(lapply(fisher, function(y) y$p.value)), OR = unlist(lapply(fisher, function(y) y$estimate)))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_overlap_sig_vs_nonSig_IR_byIntron_byGroup.csv")
df


## Are genes that have a dIR intron by age more likely to have a dIR intron by fraction?

sigIR.gene = lapply(IRcomp, function(x) unique(as.character(x[which(x$padj<=0.05),"ensID"])))
names(sigIR.gene) = paste0(names(sigIR.gene),"\nsig")
nonsigIR.gene = lapply(IRcomp, function(x) unique(as.character(x[which(x$padj>0.05),"ensID"])))
names(nonsigIR.gene) = paste0(names(nonsigIR.gene),"\nnonsig")
sigIR.gene.combined = list("Significantly IR\nBy Fraction" = c(sigIR.gene[["Adult_byFraction\nsig"]], sigIR.gene[["Prenatal_byFraction\nsig"]]),
                           "Significantly IR\nBy Age" = c(sigIR.gene[["Cytosol_byAge\nsig"]], sigIR.gene[["Nuclear_PolyA_Age\nsig"]]))
nonsigIR.gene.combined = list("Non-Significantly IR\nBy Fraction" = c(nonsigIR.gene[["Adult_byFraction\nnonsig"]], nonsigIR.gene[["Prenatal_byFraction\nnonsig"]]),
                              "Non-Significantly IR\nBy Age" = c(nonsigIR.gene[["Cytosol_byAge\nnonsig"]], nonsigIR.gene[["Nuclear_PolyA_Age\nnonsig"]]))

comps = list(c(1,3), c(1,4), c(1:2), c(2,3), c(2,4), c(3,4))
path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_"
ids = c("FractionbyAge_adult_cytosol", "FractionbyAge_adult_nucleus", "Fraction_adult_prenatal",
        "FractionbyAge_prenatal_cytosol","FractionbyAge_prenatal_nucleus","Age_cytosol_nucleus")

for (i in 1:length(comps)) {
  venn.diagram(c(sigIR.intron[comps[[i]]], nonsigIR.intron[comps[[i]]]), paste0(path, ids[i], ".jpeg"), 
               main=paste0("dIR_", ids[i]), col = "transparent",
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold",
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
}
venn.diagram(c(sigIR.intron.combined, nonsigIR.intron.combined), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_combined.jpeg", 
             main="dIR_FractionbyAge_combined", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
# all 7 saved together at ./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_byGene_overlap_Fraction_Age.pdf


tables = list()
for (i in 1:length(comps)) {
  tables[[i]] = data.frame(sig1 = c(length(sigIR.gene[[comps[[i]][1]]][sigIR.gene[[comps[[i]][1]]] %in% sigIR.gene[[comps[[i]][2]]]]),
                                    length(sigIR.gene[[comps[[i]][1]]][sigIR.gene[[comps[[i]][1]]] %in% nonsigIR.gene[[comps[[i]][2]]]])),
                           nonsig1 = c(length(nonsigIR.gene[[comps[[i]][1]]][nonsigIR.gene[[comps[[i]][1]]] %in% sigIR.gene[[comps[[i]][2]]]]),
                                       length(nonsigIR.gene[[comps[[i]][1]]][nonsigIR.gene[[comps[[i]][1]]] %in% nonsigIR.gene[[comps[[i]][2]]]])),row.names = c("sig2","nonsig2"))
}
names(tables) = ids
fisher = c(lapply(tables, fisher.test), list(combined_fraction_age = fisher.test(data.frame(
  sig1 = c(length(sigIR.gene.combined[[1]][sigIR.gene.combined[[1]] %in% sigIR.gene.combined[[2]]]),
           length(sigIR.gene.combined[[1]][sigIR.gene.combined[[1]] %in% nonsigIR.gene.combined[[2]]])),
  nonsig1 = c(length(nonsigIR.gene.combined[[1]][nonsigIR.gene.combined[[1]] %in% sigIR.gene.combined[[2]]]),
              length(nonsigIR.gene.combined[[1]][nonsigIR.gene.combined[[1]] %in% nonsigIR.gene.combined[[2]]])),
  row.names = c("sig2","nonsig2")))))
df = data.frame(pval = unlist(lapply(fisher, function(y) y$p.value)), OR = unlist(lapply(fisher, function(y) y$estimate)))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) x$p.value)),unlist(lapply(fisher, function(x) x$estimate))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_overlap_genesContaining_sig_vs_nonSig_IR_byGroup.csv")
df


## What about direction of retention in shared introns?

sigIR.dir = lapply(IRcomp, function(x) split(x, x$Sign))
sigIR.dir = unlist(lapply(sigIR.dir, function(x) lapply(x, function(y) unique(as.character(y[which(y$padj<=0.05),"intronID"])))), recursive = F)
sigIR.dir.combined = list("Significantly IR\nBy Fraction\nMore in Cytosol" = c(sigIR.dir$Adult_byFraction.MoreIRInCyt.Adult, sigIR.dir$Fetal_byFraction.MoreIRInCyt.Adult),
                          "Significantly IR\nBy Fraction\nMore in Nucleus" = c(sigIR.dir$Adult_byFraction.MoreIRInNuc.Fetal, sigIR.dir$Fetal_byFraction.MoreIRInNuc.Fetal),
                          "Significantly IR\nBy Age\nMore in Adult" = c(sigIR.dir$Cytosol_byAge.MoreIRInCyt.Adult, sigIR.dir$Nuclear_byAge.MoreIRInCyt.Adult),
                          "Significantly IR\nBy Age\nMore in Prenatal" = c(sigIR.dir$Cytosol_byAge.MoreIRInNuc.Fetal, sigIR.dir$Nuclear_byAge.MoreIRInNuc.Fetal))

comps = list(c(1,3), c(1,4), c(1,2), c(2,3), c(2,4), c(3,4))
nf = list(c(2,6), c(2,8), c(2,4), c(4,6), c(4,8), c(6,8))
ac = list(c(1,5), c(1,7), c(1,3), c(3,5), c(3,7), c(5,7))
path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRdir_"
ids = c("FractionbyAge_adult_cytosol", "FractionbyAge_adult_nucleus", "Fraction_adult_prenatal",
        "FractionbyAge_prenatal_cytosol","FractionbyAge_prenatal_nucleus","Age_cytosol_nucleus")

both = list()
for (i in 1:length(comps)) {
  both[[i]] = sigIR.intron[[comps[[i]][1]]][sigIR.intron[[comps[[i]][1]]] %in% sigIR.intron[[comps[[i]][2]]]]
}
both = c(both, list())
for (i in 1:length(nf)) {
  tables[[i]] = data.frame(MoreIRInNuc.Fetal1 = c(length(sigIR.dir[[nf[[i]][1]]][(sigIR.dir[[nf[[i]][1]]] %in% sigIR.dir[[nf[[i]][2]]]) & (sigIR.dir[[nf[[i]][1]]] %in% both[[i]])]),
                                    length(sigIR.dir[[nf[[i]][1]]][sigIR.dir[[nf[[i]][1]]] %in% sigIR.dir[[ac[[i]][2]]] & (sigIR.dir[[nf[[i]][1]]] %in% both[[i]])])),
                           MoreIRInCyt.Adult1 = c(length(sigIR.dir[[ac[[i]][1]]][sigIR.dir[[ac[[i]][1]]] %in% sigIR.dir[[nf[[i]][2]]] & (sigIR.dir[[ac[[i]][1]]] %in% both[[i]])]),
                                       length(sigIR.dir[[ac[[i]][1]]][sigIR.dir[[ac[[i]][1]]] %in% sigIR.dir[[ac[[i]][2]]] & (sigIR.dir[[ac[[i]][1]]] %in% both[[i]])])),
                           row.names = c("MoreIRInNuc.Fetal2","MoreIRInCyt.Adult2"))
}
names(tables) = names(both) = ids
fisher = c(lapply(tables, fisher.test), 
           list(combined_fraction_age = fisher.test(data.frame(sig1 = c(length(sigIR.dir.combined[[1]][sigIR.dir.combined[[1]] %in% sigIR.dir.combined[[3]]]),
                                                                        length(sigIR.dir.combined[[1]][sigIR.dir.combined[[1]] %in% sigIR.dir.combined[[4]]])),
                                                               nonsig1 = c(length(sigIR.dir.combined[[2]][sigIR.dir.combined[[2]] %in% sigIR.dir.combined[[3]]]),
                                                                           length(sigIR.dir.combined[[2]][sigIR.dir.combined[[2]] %in% sigIR.dir.combined[[4]]])),
                                                               row.names = c("sig2","nonsig2")))))
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) x$p.value)),unlist(lapply(fisher, function(x) x$estimate))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_overlap_sigIntrons_byDirection_byGroup.csv")