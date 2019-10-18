library(ggplot2)
library(GenomicRanges)
library(data.table)

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

#### Making figure 2C: comparing gene expression by age in sets regulated by fraction

fracDevel = lapply(sig[elementNROWS(sig)>0], function(x) data.frame(geneID = x$geneID, 
                                                                    Cyt.LFC = Cpres.down[match(x$geneID, rownames(Cpres.down)),"log2FoldChange"],
                                                                    Cyt.padj = Cpres.down[match(x$geneID, rownames(Cpres.down)),"padj"],
                                                                    Nuc.LFC = Npres[match(x$geneID, rownames(Npres)),"log2FoldChange"],
                                                                    Nuc.padj = Npres[match(x$geneID, rownames(Npres)),"padj"],
                                                                    Ad.LFC = Apres[match(x$geneID, rownames(Apres)),"log2FoldChange"],
                                                                    Ad.padj = Apres[match(x$geneID, rownames(Apres)),"padj"],
                                                                    Fet.LFC = Fpres.down[match(x$geneID, rownames(Fpres.down)),"log2FoldChange"],
                                                                    Fet.padj = Fpres.down[match(x$geneID, rownames(Fpres.down)),"padj"]))
fracDevel = do.call(rbind, Map(cbind, fracDevel, fracReg = as.list(c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nPrenatal Only", "Nuclear:\nAdult Only", 
                                                                     "Cytoplasmic:\nPrenatal Only", "Cytoplasmic:\nAdult Only","Nuclear: Adult/\nCytoplasmic: Prenatal", "Interaction"))))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nAdult Only", "Nuclear:\nPrenatal Only",
                                      "Cytoplasmic:\nAdult Only","Cytoplasmic:\nPrenatal Only","Nuclear: Adult/\nCytoplasmic: Prenatal", "Interaction"))
fracDevel$Diff = fracDevel$Cyt.LFC - fracDevel$Nuc.LFC
fracDevel = fracDevel[-which(fracDevel$fracReg %in% c("Nuclear: Adult/\nCytoplasmic: Prenatal", "Interaction")),]
fracDevel[grep("Nuclear", fracDevel$fracReg),"Dir"] = "Nuclear"
fracDevel[grep("Cytoplasmic", fracDevel$fracReg),"Dir"] = "Cytoplasmic"
fracDevel[grep("Both", fracDevel$fracReg),"Age"] = "In Both Ages"
fracDevel[grep("Adult Only", fracDevel$fracReg),"Age"] = "In Adult Only"
fracDevel[grep("Prenatal Only", fracDevel$fracReg),"Age"] = "In Prenatal Only"
head(fracDevel)

df = fracDevel[which(fracDevel$fracReg=="Nuclear:\nAdult Only" | fracDevel$fracReg=="Cytoplasmic:\nAdult Only"),]
head(df)
df$fracReg = gsub("Nuclear:\nAdult Only", "Nuclear", df$fracReg)
df$fracReg = gsub("Cytoplasmic:\nAdult Only", "Cytoplasmic", df$fracReg)

means = data.frame(group = c("Nuclear","Cytoplasmic"), 
                   diffMean = c(mean(df[grep("Nuclear", df$fracReg),"Diff"]), mean(df[grep("Cytoplasmic", df$fracReg),"Diff"])))

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/AdultRetainednExported_Age-LFC_Diff.pdf", width=5.5, height=5)
ggplot(df, aes(Diff, fill = fracReg, colour = fracReg)) + geom_density(alpha=0.4) + facet_grid(. ~ Age) +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") + 
  geom_vline(data=means, aes(xintercept=diffMean, color=group), linetype="dashed", size=1) +
  ylab("Density") + xlab("(Age LFC in Cyt.) - (Age LFC in Nuc.)") + theme_classic() +
  ggtitle("") + guides(colour=FALSE) + labs(fill="Gene Set") +
  theme(title = element_text(size = 20), text = element_text(size = 20),
        legend.position = c(0.8, 0.75), legend.background = element_rect(fill = "transparent"))
dev.off()


t.test(x = df[grep("Nuclear", df$fracReg),"Diff"], df[grep("Cytoplasmic", df$fracReg),"Diff"])
#t = 100.18, df = 3537.8, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1.370287 1.424995
#sample estimates:
#  mean of x  mean of y 
#0.7099931 -0.6876481


#### Making figure 2D: comparing gene expression by fraction in sets regulated by age

fracDevel = lapply(age.sig[elementNROWS(age.sig)>0], function(x) data.frame(geneID = x$geneID, 
                                                                    Fet.LFC = Fpres.down[match(x$geneID, rownames(Fpres.down)),"log2FoldChange"],
                                                                    Fet.padj = Fpres.down[match(x$geneID, rownames(Fpres.down)),"padj"],
                                                                    Ad.LFC = Apres[match(x$geneID, rownames(Apres)),"log2FoldChange"],
                                                                    Ad.padj = Apres[match(x$geneID, rownames(Apres)),"padj"]))
fracDevel = do.call(rbind, Map(cbind, fracDevel, fracReg = as.list(c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                                                                     "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only",
                                                                     "Decreasing: Nucleus/\nIncreasing: Cytoplasm", 
                                                                     "Decreasing: Cytoplasm/\nIncreasing: Nucleus","Interaction"))))
fracDevel = fracDevel[-which(fracDevel$fracReg %in% c("Decreasing: Nucleus/\nIncreasing: Cytoplasm", 
                                                      "Decreasing: Cytoplasm/\nIncreasing: Nucleus","Interaction")),]
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                                      "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only"))
fracDevel$Diff = fracDevel$Ad.LFC - fracDevel$Fet.LFC
fracDevel[grep("Decreasing", fracDevel$fracReg),"Dir"] = "Decreasing"
fracDevel[grep("Increasing", fracDevel$fracReg),"Dir"] = "Increasing"
fracDevel$Dir = factor(fracDevel$Dir, levels = c("Increasing", "Decreasing"))
fracDevel[grep("Both", fracDevel$fracReg),"Age"] = "In Both Fractions"
fracDevel[grep("Cytoplasm Only", fracDevel$fracReg),"Age"] = "In Cytoplasm Only"
fracDevel[grep("Nucleus Only", fracDevel$fracReg),"Age"] = "In Nucleus Only"
head(fracDevel)

means = data.frame(Dir = c("Decreasing","Increasing","Decreasing","Increasing","Decreasing","Increasing"),
                   Age = c("In Both Fractions","In Both Fractions","In Cytoplasm Only","In Cytoplasm Only","In Nucleus Only","In Nucleus Only"),
                   diffMean = c(mean(na.omit(fracDevel[which(fracDevel$fracReg=="Decreasing: Both"),"Diff"])), 
                                mean(na.omit(fracDevel[which(fracDevel$fracReg=="Increasing: Both"),"Diff"])),
                                mean(na.omit(fracDevel[which(fracDevel$fracReg=="Decreasing:\nCytoplasm Only"),"Diff"])),
                                mean(na.omit(fracDevel[which(fracDevel$fracReg=="Increasing:\nCytoplasm Only"),"Diff"])),
                                mean(na.omit(fracDevel[which(fracDevel$fracReg=="Decreasing:\nNucleus Only"),"Diff"])),
                                mean(na.omit(fracDevel[which(fracDevel$fracReg=="Increasing:\nNucleus Only"),"Diff"]))))

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/AdultIncrDecr_Fraction-LFC_Diff.pdf", width=10, height=4)
ggplot(fracDevel) + geom_density(aes(Diff, fill = Dir, colour = Dir), alpha=0.4) + 
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + 
  geom_vline(data=means, aes(xintercept=diffMean, color=Dir), linetype="dashed", size=1) + facet_grid(. ~ Age) +
  ylab("Density") + xlab("(Fraction LFC in Adult) - (Fraction LFC in Prenatal)") + theme_classic() +
  ggtitle("") + guides(colour=FALSE) + labs(fill="Gene Set") +
  theme(title = element_text(size = 20), text = element_text(size = 20),
        legend.background = element_rect(fill = "transparent"))
dev.off()


## LET'S TRY AGAIN

fracDevel = lapply(sig[elementNROWS(sig)>0], function(x) data.frame(geneID = x$geneID, 
                                                                    Cyt.LFC = Cpres.down[match(x$geneID, rownames(Cpres.down)),"log2FoldChange"],
                                                                    Nuc.LFC = Npres[match(x$geneID, rownames(Npres)),"log2FoldChange"],
                                                                    Ad.LFC = Apres[match(x$geneID, rownames(Apres)),"log2FoldChange"],
                                                                    Fet.LFC = Fpres.down[match(x$geneID, rownames(Fpres.down)),"log2FoldChange"]))
fracDevel = do.call(rbind, Map(cbind, fracDevel[!names(fracDevel) %in% c("ret_Ad_exp_Fet","interacting")], 
                               fracReg = as.list(c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nPrenatal Only", "Nuclear:\nAdult Only", 
                                                   "Cytoplasmic:\nPrenatal Only", "Cytoplasmic:\nAdult Only"))))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nAdult Only", "Nuclear:\nPrenatal Only",
                                      "Cytoplasmic:\nAdult Only","Cytoplasmic:\nPrenatal Only"))
fracDevel[grep("Nuclear", fracDevel$fracReg),"Dir"] = "Nuclear"
fracDevel[grep("Cytoplasmic", fracDevel$fracReg),"Dir"] = "Cytoplasmic"
fracDevel[grep("Both", fracDevel$fracReg),"Frac"] = "In Both Fractions"
fracDevel[grep("Adult", fracDevel$fracReg),"Frac"] = "In Adult Only"
fracDevel[grep("Prenatal", fracDevel$fracReg),"Frac"] = "In Prenatal Only"
fracDevel$Frac = factor(fracDevel$Frac, levels = c("In Both Fractions", "In Adult Only","In Prenatal Only"))
head(fracDevel)

df = fracDevel[which(fracDevel$fracReg=="Nuclear:\nAdult Only" | fracDevel$fracReg=="Cytoplasmic:\nAdult Only"),]
df$fracReg = gsub("Nuclear:\nAdult Only", "Nuclear", df$fracReg)
df$fracReg = gsub("Cytoplasmic:\nAdult Only", "Cytoplasmic", df$fracReg)

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/Fig2C_scatterplot.pdf", width=9, height=4)

ggplot(fracDevel, aes(x=Ad.LFC,y=Cyt.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Dark2") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) +  facet_grid(. ~ Frac) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Cytoplasm)") + xlab("Nucleus vs Cytoplasm (in Adult)") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())

ggplot(fracDevel, aes(x=Fet.LFC,y=Cyt.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Dark2") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) +  facet_grid(. ~ Frac) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Cytoplasm)") + xlab("Nucleus vs Cytoplasm (in Prenatal)") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())

ggplot(fracDevel, aes(x=Ad.LFC,y=Nuc.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Dark2") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) +  facet_grid(. ~ Frac) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Nucleus)") + xlab("Nucleus vs Cytoplasm (in Adult)") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())

ggplot(fracDevel, aes(x=Fet.LFC,y=Nuc.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Dark2") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) +  facet_grid(. ~ Frac) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Nucleus)") + xlab("Nucleus vs Cytoplasm (in Prenatal)") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())
dev.off()


fracDevel = lapply(age.sig[elementNROWS(age.sig)>0], function(x) data.frame(geneID = x$geneID, 
                                                                            Fet.LFC = Fpres.down[match(x$geneID, rownames(Fpres.down)),"log2FoldChange"],
                                                                            Ad.LFC = Apres[match(x$geneID, rownames(Apres)),"log2FoldChange"],
                                                                            Cyt.LFC = Cpres.down[match(x$geneID, rownames(Cpres.down)),"log2FoldChange"],
                                                                            Nuc.LFC = Npres[match(x$geneID, rownames(Npres)),"log2FoldChange"]))
fracDevel = do.call(rbind, Map(cbind, fracDevel[-which(names(fracDevel) %in% c("decr_Nuc_incr_Cyt","decr_Cyt_incr_Nuc","interacting"))], 
                               fracReg = as.list(c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                                                   "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only"))))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                                      "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only"))
fracDevel[grep("Decreasing", fracDevel$fracReg),"Dir"] = "Decreasing"
fracDevel[grep("Increasing", fracDevel$fracReg),"Dir"] = "Increasing"
fracDevel$Dir = factor(fracDevel$Dir, levels = c("Increasing", "Decreasing"))
fracDevel[grep("Both", fracDevel$fracReg),"Age"] = "In Both Fractions"
fracDevel[grep("Cytoplasm Only", fracDevel$fracReg),"Age"] = "In Cytoplasm Only"
fracDevel[grep("Nucleus Only", fracDevel$fracReg),"Age"] = "In Nucleus Only"
head(fracDevel)


pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/Fig2D_scatterplot.pdf", width=9, height=4)

ggplot(fracDevel, aes(x=Ad.LFC,y=Cyt.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Set1") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) + facet_grid(. ~ Age) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Cytoplasm)") + xlab("Cytoplasm vs Nucleus (in Adult)") +
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())

ggplot(fracDevel, aes(x=Fet.LFC,y=Cyt.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Set1") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) + facet_grid(. ~ Age) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Cytoplasm)") + xlab("Cytoplasm vs Nucleus (in Prenatal)") +
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())

ggplot(fracDevel, aes(x=Ad.LFC,y=Nuc.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Set1") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) + facet_grid(. ~ Age) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Nucleus)") + xlab("Cytoplasm vs Nucleus (in Adult)") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())

ggplot(fracDevel, aes(x=Fet.LFC,y=Nuc.LFC)) + geom_point(aes(colour = Dir), alpha=0.5) +
  scale_color_brewer(palette="Set1") + geom_smooth(method="lm",size=1, colour="black", se=FALSE) + facet_grid(. ~ Age) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab("Prenatal vs Adult\n(in Nucleus)") + xlab("Cytoplasm vs Nucleus (in Prenatal)") +
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())
dev.off()


df = cbind(rbind(fracDevel, fracDevel), FracLFC = c(fracDevel$Fet.LFC, fracDevel$Ad.LFC), 
           Frac = c(rep.int("In Prenatal",nrow(fracDevel)),rep.int("In Adult", nrow(fracDevel))))
head(df)

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/Fig2D_scatterplot_grid.pdf", width=7, height=5)
ggplot(df, aes(x=Cyt.LFC,y=FracLFC)) + geom_point(aes(colour = Dir)) +
  scale_color_brewer(palette="Set1") + facet_grid(Frac ~ Age) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") + theme_classic() +
  xlab(expression(paste(log[2], (Prenatal/Adult)))) + 
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  ggtitle("Developmentally-Regulated Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position="top",
        legend.background = element_rect(fill = "transparent"), 
        legend.title = element_blank())
dev.off()

