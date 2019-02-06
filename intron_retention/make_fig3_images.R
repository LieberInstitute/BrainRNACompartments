library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23_nodownsamp.rda")

# Load IRFinder Results
names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
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

## Remake 3D: Plot IR by fraction LFC

head(irbyGene[[1]])
df = do.call(rbind, Map(cbind, irbyGene, sampleID = as.list(names(irbyGene))))
df$rownum = c(1:nrow(df))
df$Fraction = ifelse(df$rownum %in% grep("C", df$sampleID), "Cytoplasm", "Nucleus")
df$Age = ifelse(df$rownum %in% grep("Br53", df$sampleID), "Prenatal", "Adult")
df$genecodeID = geneMap[match(df$genes, geneMap$ensemblID),"gencodeID"]
dt = data.table(df)
dt = dt[,list(meanIR = mean(IRratio), medianIR = median(IRratio)), by=c("Fraction","Age","genecodeID")]
df = data.frame(dt)

df = cbind(rbind(df,df), LFC=c(data.frame(Apres)[match(df$genecodeID, rownames(Apres)),"log2FoldChange"],
                               data.frame(Fpres.down)[match(df$genecodeID, rownames(Fpres.down)),"log2FoldChange"]),
           padj=c(data.frame(Apres)[match(df$genecodeID, rownames(Apres)),"padj"],
                  data.frame(Fpres.down)[match(df$genecodeID, rownames(Fpres.down)),"padj"]),
           LFC.comp = c(rep.int("In Adult", nrow(df)), rep.int("In Prenatal", nrow(df))))
head(df)
df$group = paste0(df$Age, "\n", df$Fraction)
df$group = factor(df$group, levels = c("Adult\nCytoplasm","Prenatal\nCytoplasm","Adult\nNucleus","Prenatal\nNucleus"))

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/FractionLFC_byIRratio.pdf", width=7.75,height=4.25)
ggplot(df, aes(x=meanIR,y=LFC)) + geom_point(aes(colour = LFC.comp)) +
  facet_grid(Fraction ~ Age) + scale_color_brewer(palette="Set1") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") + theme_classic() +
  xlab("IR Ratio") + 
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  ggtitle("IR Ratio and Fraction Expression Changes") +
  theme(plot.title = element_text(hjust = 0.5),title = element_text(size = 20), text = element_text(size = 20), legend.position = "top",
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank())
dev.off()


## Try again: plot max IR per gene by RPKM in each sample

geneRpkm = geneRpkm[,grep("polyA", colnames(geneRpkm))]
colnames(geneRpkm) = gsub("_polyA", "", colnames(geneRpkm))
head(geneRpkm)
df$RPKM = NA

for (i in 1:ncol(geneRpkm)) {
  val = geneRpkm[,i]
  for (j in 1:length(val)) {
    df[which(as.character(df$sampleID)==colnames(geneRpkm)[i] & df$genecodeID==names(val)[j]),"RPKM"] = val[j]
  }
}
save(df, geneMap, file = "./Dropbox/sorted_figures/github_controlled/Intron_retention/data/gene_IR_comparisons/IR_df_geneMap.rda")

load("./Dropbox/sorted_figures/github_controlled/Intron_retention/data/gene_IR_comparisons/IR_df_geneMap.rda")
head(df)
res = list(AC = cor.test(x=df[which(df$group=="Adult\nCytoplasm"),"IRratio"], y=df[which(df$group=="Adult\nCytoplasm"),"RPKM"]),
           AN = cor.test(x=df[which(df$group=="Adult\nNucleus"),"IRratio"], y=df[which(df$group=="Adult\nNucleus"),"RPKM"]),
           PC = cor.test(x=df[which(df$group=="Prenatal\nCytoplasm"),"IRratio"], y=df[which(df$group=="Prenatal\nCytoplasm"),"RPKM"]),
           PN = cor.test(x=df[which(df$group=="Prenatal\nNucleus"),"IRratio"], y=df[which(df$group=="Prenatal\nNucleus"),"RPKM"]))

res = do.call(rbind, Map(cbind, Group = as.list(names(res)), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            cor = x$estimate))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#Group      Tstat         pval         cor          FDR
#1    AC -13.860367 1.223290e-43 -0.04655218 4.893160e-43
#2    AN -13.724945 7.965300e-43 -0.04610403 1.593060e-42
#3    PC  -8.465848 2.581200e-17 -0.02846799 3.441600e-17
#4    PN  -6.716980 1.866488e-11 -0.02259866 1.866488e-11

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/FractionLFC_byIRratio3.pdf", width=6.17,height=3.75)
ggplot(df, aes(x=IRratio,y=log(RPKM+1))) + geom_point(aes(colour = group)) +
  geom_smooth(aes(x=IRratio,y=log(RPKM+1)), method='lm', colour="black", se=TRUE) + facet_grid(Fraction ~ Age) + 
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") +
  ggtitle("IR Ratio and Expression") +
  theme(plot.title = element_text(hjust = 0.5),title = element_text(size = 20), 
        text = element_text(size = 20), legend.position="none")
dev.off()


## Remake Figure 3D Venn Diagram

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
full = list(adult = read.table("./Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-adultShared.txt", header = TRUE),
            prenatal = read.table("./Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt", header = TRUE),
            cytosol = read.table("./Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt", header = TRUE),
            nucleus = read.table("./Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt", header = TRUE))
total = as.list(elementNROWS(full))
names(total) = names(IRcomp)
IRcomp = Map(cbind, IRcomp, padj = mapply(function(p,t) p.adjust(p, method = "fdr", n = t), lapply(IRcomp, function(x) x$p.diff), total))
elementNROWS(IRcomp)

sigIR.gene = lapply(IRcomp, function(x) unique(as.character(x[which(x$padj<=0.05),"ensID"])))
names(sigIR.gene) = c("By Fraction\nIn Adult","By Fraction\nIn Prenatal", 
                      "By Age\nin Cytoplasm", "By Age\nin Nucleus")

venn.diagram(sigIR.gene, 
             "./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_gene_FractionbyAge.jpeg", 
             main="Significantly differentially retained introns by group", col = "transparent",
             fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.3,
             fill = c("#e41a1c", "#377eb8","#1b9e77","#d95f02"), alpha = 0.50,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(sigIR.gene, 
             "./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_gene_FractionbyAge2.jpeg", 
             main="Significantly differentially retained introns by group",
             fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.3,
             cat.fontfamily = "Arial", margin=0.2)


## Try two different options for fig 3F:

# remake it as a scatterplot like fig 2

sigIR = lapply(IRcomp, function(x) x[which(x$padj<=0.05),])
sigIR = do.call(rbind, Map(cbind, sigIR, comp = as.list(c("By Fraction\nIn Adult","By Fraction\nIn Prenatal",
                                                          "By Age\nIn Cytoplasm", "By Age\nIn Nucleus"))))
x = reshape2::melt(sigIR, measure.vars=c("A1.IRratio","A2.IRratio","A3.IRratio","B1.IRratio","B2.IRratio","B3.IRratio"))
x$genecodeID = geneMap[match(x$ensID, geneMap$ensemblID),"gencodeID"]
x$"In Adult.LFC" = Apres[match(x$genecodeID, rownames(Apres)),"log2FoldChange"]
x$"In Prenatal.LFC" = Fpres.down[match(x$genecodeID, rownames(Fpres.down)),"log2FoldChange"]
x$"In Nucleus.LFC" = Npres[match(x$genecodeID, rownames(Npres)),"log2FoldChange"]
x$"In Cytoplasm.LFC" = Cpres.down[match(x$genecodeID, rownames(Cpres.down)),"log2FoldChange"]
x$"In Adult.padj" = Apres[match(x$genecodeID, rownames(Apres)),"padj"]
x$"In Prenatal.padj" = Fpres.down[match(x$genecodeID, rownames(Fpres.down)),"padj"]
x$"In Nucleus.padj" = Npres[match(x$genecodeID, rownames(Npres)),"padj"]
x$"In Cytoplasm.padj" = Cpres.down[match(x$genecodeID, rownames(Cpres.down)),"padj"]
x = reshape2::melt(x, measure.vars=c("In Adult.LFC","In Prenatal.LFC","In Nucleus.LFC","In Cytoplasm.LFC"), variable_name="LFC")
colnames(x)[33:34] = c("IRratio.Sample","IRratio")
colnames(x)[40:41] = c("LFC.location","LFC")
x$LFC.location = gsub(".LFC", "", x$LFC.location)
x$LFC.comp = ifelse(x$LFC.location %in% c("In Adult","In Prenatal"),"By Fraction","By Age")
x$IR.location = x$comp
x$IR.comp = ifelse(x$comp %in% c("By Fraction\nIn Adult","By Fraction\nIn Prenatal"), "By Fraction", "By Age")
x$IR.location = gsub("By Fraction\n","",x$IR.location, fixed=T)
x$IR.location = gsub("By Age\n","",x$IR.location, fixed=T)
x$IR.location = factor(x$IR.location, levels = c("In Adult","In Prenatal","In Cytoplasm","In Nucleus"))
x$LFC.location = factor(x$LFC.location, levels = c("In Adult","In Prenatal","In Cytoplasm","In Nucleus"))
head(x)
unique(x$LFC.location)

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_by_LFC.pdf", width=6.5,height=5)
ggplot(x[which(x$LFC.comp=="By Fraction"),], aes(x=IRratio,y=LFC)) + geom_point(aes(colour = IR.location)) +
  facet_grid(LFC.location ~ IR.comp) +
  scale_color_manual(values = c("#e41a1c", "#377eb8","#1b9e77","#d95f02")) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") + theme_classic() +
  xlab("IR Ratio") + 
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  ggtitle("IR Ratio and Fraction Expression") +
  theme(plot.title = element_text(hjust = 0.5),title = element_text(size = 20), text = element_text(size = 20),
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank(), legend.position = "bottom")

ggplot(x[which(x$LFC.comp=="By Age"),], aes(x=IRratio,y=LFC)) + geom_point(aes(colour = IR.location)) +
  facet_grid(LFC.location ~ IR.comp) +
  scale_color_manual(values = c("#e41a1c", "#377eb8","#1b9e77","#d95f02")) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") + theme_classic() +
  xlab("IR Ratio") + 
  ylab(expression(paste(log[2], (Prenatal/Adult)))) +
  ggtitle("IR Ratio and Age Expression") +
  theme(plot.title = element_text(hjust = 0.5),title = element_text(size = 20), text = element_text(size = 20),
        legend.background = element_rect(fill = "transparent"), legend.title = element_blank(), legend.position = "bottom")
dev.off()


## Try again

load("./Dropbox/sorted_figures/github_controlled/Intron_retention/data/gene_IR_comparisons/IR_df_geneMap.rda")
sigIR = lapply(IRcomp, function(x) x[which(x$padj<=0.05),])
sigIR = do.call(rbind, Map(cbind, sigIR, comp = as.list(c("By Fraction\nIn Adult","By Fraction\nIn Prenatal",
                                                          "By Age\nIn Cytoplasm", "By Age\nIn Nucleus"))))
sigIR$ID = paste0("chr",sigIR$Chr, ":", sigIR$Start, "-", sigIR$End, ":(", sigIR$Direction, ")")
IRres = Map(cbind, IRres, intronID = lapply(IRres, function(x) paste0("chr",x$Chr, ":", x$Start, "-", x$End, ":(", x$Direction, ")")))
head(IRres[[1]])
x = lapply(IRres, function(z) data.frame(sigIR, IR = z[match(sigIR$ID, z$intronID),"IRratio"]))
tail(x)
x = do.call(rbind, Map(cbind, x, sampleID = as.list(names(x))))
x$genecodeID = geneMap[match(x$ensID, geneMap$ensemblID), "gencodeID"]
x$samp.IRID = paste0(x$genecodeID,"/",x$sampleID)
df$samp.IRID = paste0(df$genecodeID,"/",df$sampleID)
x$RPKM = df[match(x$samp.IRID, df$samp.IRID),"RPKM"]
x$rownum = 1:nrow(x)
x$Fraction = ifelse(x$rownum %in% grep("C", x$sampleID), "Cytoplasm", "Nucleus")
x$Age = ifelse(x$rownum %in% grep("Br53", x$sampleID), "Prenatal", "Adult")
x$group = paste0(x$Age, "\n", x$Fraction)
x$group = factor(x$group, levels = c("Adult\nCytoplasm","Prenatal\nCytoplasm","Adult\nNucleus","Prenatal\nNucleus"))

x$group = paste0(x$Age, ":", x$Fraction)
x$group = factor(x$group, levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm","Adult:Nucleus","Prenatal:Nucleus"))

x$compLocation = x$shortComp = NA
x[grep("By Fraction", x$comp),"shortComp"] = "By Fraction"
x[grep("By Age", x$comp),"shortComp"] = "By Age"

x[grep("In Adult", x$comp),"compLocation"] = "In Adult"
x[grep("In Prenatal", x$comp),"compLocation"] = "In Prenatal"
x[grep("In Cytoplasm", x$comp),"compLocation"] = "In Cytoplasm"
x[grep("In Nucleus", x$comp),"compLocation"] = "In Nucleus"

x$quad = NA
x[which(x$comp=="By Fraction\nIn Adult" & x$Age=="Adult"),"quad"] = "IR By Fraction:\nin Adult"
x[which(x$comp=="By Fraction\nIn Prenatal" & x$Age=="Prenatal"),"quad"] = "IR By Fraction:\nin Prenatal"
x[which(x$comp=="By Age\nIn Cytoplasm" & x$Fraction=="Cytoplasm"),"quad"] = "IR By Age:\nin Cytoplasm"
x[which(x$comp=="By Age\nIn Nucleus" & x$Fraction=="Nucleus"),"quad"] = "IR By Age:\nin Nucleus"
xs = x[which(x$quad %in% c("IR By Fraction:\nin Adult","IR By Age:\nin Cytoplasm","IR By Age:\nin Nucleus","IR By Fraction:\nin Prenatal")),]
head(xs)

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_byRPKM_facet_wrap.pdf", width=3,height=3.5)
ggplot(x[which(x$comp=="By Fraction\nIn Adult" & x$Age=="Adult"),], aes(x=IR,y=log(RPKM+1))) + 
  geom_point(aes(colour = group)) + geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) +
  scale_color_manual(values = c("#F8766D","#00BFC4")) + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("IR By Fraction:\nin Adult") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="none")
ggplot(x[which(x$comp=="By Fraction\nIn Prenatal" & x$Age=="Prenatal"),], aes(x=IR,y=log(RPKM+1))) + 
  geom_point(aes(colour = group)) + geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) +
  scale_color_manual(values = c("#7CAE00","#C77CFF")) + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("IR By Fraction:\nin Prenatal") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="none")
ggplot(x[which(x$comp=="By Age\nIn Cytoplasm" & x$Fraction=="Cytoplasm"),], aes(x=IR,y=log(RPKM+1))) + 
  geom_point(aes(colour = group)) + geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) +
  scale_color_manual(values = c("#F8766D", "#7CAE00")) + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("IR By Age:\nin Cytoplasm") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="none")
ggplot(x[which(x$comp=="By Age\nIn Nucleus" & x$Fraction=="Nucleus"),], aes(x=IR,y=log(RPKM+1))) + 
  geom_point(aes(colour = group)) + geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) +
  scale_color_manual(values = c("#00BFC4","#C77CFF")) + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("IR By Age:\nin Nucleus") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="none")
dev.off()

res = list(fracIR.inAdult.AC = cor.test(x = x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Adult\nCytoplasm"),"IR"],
                                   y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Adult\nCytoplasm"),"RPKM"]+1)),
           fracIR.inAdult.AN = cor.test(x = x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Adult\nNucleus"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Adult\nNucleus"),"RPKM"]+1)),
           fracIR.inAdult.PC = cor.test(x = x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Prenatal\nCytoplasm"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Prenatal\nCytoplasm"),"RPKM"]+1)),
           fracIR.inAdult.PN = cor.test(x = x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Prenatal\nNucleus"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Prenatal\nNucleus"),"RPKM"]+1)),
           fracIR.inPrenatal.AC = cor.test(x = x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Adult\nCytoplasm"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Adult\nCytoplasm"),"RPKM"]+1)),
           fracIR.inPrenatal.AN = cor.test(x = x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Adult\nNucleus"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Adult\nNucleus"),"RPKM"]+1)),
           fracIR.inPrenatal.PC = cor.test(x = x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Prenatal\nCytoplasm"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Prenatal\nCytoplasm"),"RPKM"]+1)),
           fracIR.inPrenatal.PN = cor.test(x = x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Prenatal\nNucleus"),"IR"],
                                        y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Prenatal\nNucleus"),"RPKM"]+1)),
           ageIR.inCyt.AC = cor.test(x = x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Adult\nCytoplasm"),"IR"],
                                        y = log(x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Adult\nCytoplasm"),"RPKM"]+1)),
           ageIR.inCyt.AN = cor.test(x = x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Adult\nNucleus"),"IR"],
                                        y = log(x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Adult\nNucleus"),"RPKM"]+1)),
           ageIR.inCyt.PC = cor.test(x = x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Prenatal\nCytoplasm"),"IR"],
                                        y = log(x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Prenatal\nCytoplasm"),"RPKM"]+1)),
           ageIR.inCyt.PN = cor.test(x = x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Prenatal\nNucleus"),"IR"],
                                        y = log(x[which(x$comp=="By Age\nIn Cytoplasm" & x$group=="Prenatal\nNucleus"),"RPKM"]+1)),
           ageIR.inNuc.AC = cor.test(x = x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Adult\nCytoplasm"),"IR"],
                                           y = log(x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Adult\nCytoplasm"),"RPKM"]+1)),
           ageIR.inNuc.AN = cor.test(x = x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Adult\nNucleus"),"IR"],
                                           y = log(x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Adult\nNucleus"),"RPKM"]+1)),
           ageIR.inNuc.PC = cor.test(x = x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Prenatal\nCytoplasm"),"IR"],
                                           y = log(x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Prenatal\nCytoplasm"),"RPKM"]+1)),
           ageIR.inNuc.PN = cor.test(x = x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Prenatal\nNucleus"),"IR"],
                                           y = log(x[which(x$comp=="By Age\nIn Nucleus" & x$group=="Prenatal\nNucleus"),"RPKM"]+1)))
           
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            cor = x$estimate))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res[res$FDR<=0.05,]
#             Comparison     Tstat         pval        cor          FDR
#1     fracIR.inAdult.AC -5.739586 9.809234e-08 -0.4940895 1.569478e-06
#2     fracIR.inAdult.AN -5.016059 2.208717e-06 -0.4430829 1.766973e-05
#3     fracIR.inAdult.PC -4.047956 1.001741e-04 -0.3704752 4.006964e-04
#5  fracIR.inPrenatal.AC -2.549632 2.141917e-02 -0.5375023 4.895809e-02
#11       ageIR.inCyt.PC -5.063597 2.330702e-05 -0.6913766 1.243041e-04
#12       ageIR.inCyt.PN -3.171129 3.662816e-03 -0.5140460 9.767510e-03
#15       ageIR.inNuc.PC -3.516725 8.311205e-04 -0.4105698 2.659586e-03

table(sigIR[which(sigIR$comp=="By Age\nIn Cytoplasm"),"Sign"]=="MoreIRInNuc.Fetal")
#FALSE  TRUE 
#    1     9 
table(sigIR[which(sigIR$comp=="By Age\nIn Nucleus"),"Sign"]=="MoreIRInNuc.Fetal")
#FALSE  TRUE 
#    2    19 
table(sigIR[which(sigIR$comp=="By Fraction\nIn Adult"),"Sign"]=="MoreIRInNuc.Fetal")
#TRUE 
#  35 
table(sigIR[which(sigIR$comp=="By Fraction\nIn Prenatal"),"Sign"]=="MoreIRInNuc.Fetal")
#TRUE 
#   6

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_byRPKM_facet_wrap.pdf", width=4.75,height=5.75)
ggplot(xs, aes(x=IR,y=log(RPKM+1))) + geom_point(aes(colour = group)) +
  geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) + 
  facet_wrap(. ~ quad) + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ggtitle("") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="bottom", legend.title = element_blank()) +
  guides(colour = guide_legend(nrow = 2))
dev.off()

xs$quad = gsub("IR By Fraction:\ni", "I", xs$quad)
xs$quad = gsub("IR By Age:\ni", "I", xs$quad)

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_byRPKM_facet_grid.pdf", width=4.75,height=4)
ggplot(xs[which(xs$shortComp=="By Fraction"),], aes(x=IR,y=log(RPKM+1))) + geom_point(aes(colour = group)) +
  geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) + 
  facet_grid(. ~ quad, scales="free") + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ggtitle("IR By Fraction") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="bottom", legend.title = element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA), nrow=2)) 
ggplot(xs[which(xs$shortComp=="By Age"),], aes(x=IR,y=log(RPKM+1))) + geom_point(aes(colour = group)) +
  geom_smooth(aes(x=IR,y=log(RPKM+1), color=group), method='lm', se=TRUE) + 
  facet_grid(. ~ quad, scales="free") + theme_classic() +
  xlab("IR Ratio") + ylab("log(RPKM+1)") + ggtitle("IR By Age") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="bottom", legend.title = element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA), nrow=2))
dev.off()

## Try again as a boxplot

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_byRPKM2.pdf", width=5,height=3.5)
ggplot(x[which(x$shortComp=="By Fraction"),], aes(x=Fraction,y=log(RPKM+1))) + facet_grid(. ~ compLocation) +
  geom_boxplot(aes(fill = Age)) + scale_fill_brewer(palette="Set1") + theme_classic() + 
  xlab("") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("Differential IR by Fraction") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())
ggplot(x[which(x$shortComp=="By Fraction"),], aes(x=Age,y=log(RPKM+1))) + facet_grid(. ~ compLocation) +
  geom_boxplot(aes(fill = Fraction)) + scale_fill_brewer(palette="Dark2") + theme_classic() + 
  xlab("") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("Differential IR by Fraction") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())
ggplot(x[which(x$shortComp=="By Age"),], aes(x=Age,y=log(RPKM+1))) + facet_grid(. ~ compLocation) +
  geom_boxplot(aes(fill = Fraction)) + scale_fill_brewer(palette="Dark2") + theme_classic() + 
  xlab("") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("Differential IR by Age") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())
ggplot(x[which(x$shortComp=="By Age"),], aes(x=Fraction,y=log(RPKM+1))) + facet_grid(. ~ compLocation) +
  geom_boxplot(aes(fill = Age)) + scale_fill_brewer(palette="Set1") + theme_classic() + 
  xlab("") + ylab("log(RPKM+1)") + ylim(0,6) + ggtitle("Differential IR by Age") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())
dev.off()

head(x)
res = list(fracIR.inAdult = t.test(x = log(x[which(x$comp=="By Fraction\nIn Adult" & x$Age=="Adult"),"RPKM"]+1),
                                   y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$Age=="Prenatal"),"RPKM"]+1)),
           fracIR.inAdult.cyt = t.test(x = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Adult\nCytoplasm"),"RPKM"]+1),
                                       y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Prenatal\nCytoplasm"),"RPKM"]+1)),
           fracIR.inAdult.nuc = t.test(x = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Adult\nNucleus"),"RPKM"]+1),
                                       y = log(x[which(x$comp=="By Fraction\nIn Adult" & x$group=="Prenatal\nNucleus"),"RPKM"]+1)),
           fracIR.inPren = t.test(x = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$Age=="Adult"),"RPKM"]+1),
                                  y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$Age=="Prenatal"),"RPKM"]+1)),
           fracIR.inPren.cyt = t.test(x = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Adult\nCytoplasm"),"RPKM"]+1),
                                       y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Prenatal\nCytoplasm"),"RPKM"]+1)),
           fracIR.inPren.nuc = t.test(x = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Adult\nNucleus"),"RPKM"]+1),
                                         y = log(x[which(x$comp=="By Fraction\nIn Prenatal" & x$group=="Prenatal\nNucleus"),"RPKM"]+1)))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            adultMean = x$estimate[1], prenatalMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#          Comparison     Tstat         pval adultMean prenatalMean          FDR
#1     fracIR.inAdult  6.761614 5.584195e-11  3.723575     3.070832 3.350517e-10
#2 fracIR.inAdult.cyt  5.284258 3.510181e-07  3.756218     2.992468 1.053054e-06
#3 fracIR.inAdult.nuc  4.217061 4.037004e-05  3.691244     3.151500 8.074008e-05
#4      fracIR.inPren -2.666788 9.573762e-03  4.351757     4.740361 1.436064e-02
#5  fracIR.inPren.cyt -2.469187 1.907076e-02  4.178975     4.600104 2.288491e-02
#6  fracIR.inPren.nuc -1.498265 1.435154e-01  4.524540     4.865034 1.435154e-01


## look at the LFC

sigIR = lapply(IRcomp, function(x) x[which(x$padj<=0.05),])
sigIR = do.call(rbind, Map(cbind, sigIR, comp = as.list(c("By Fraction\nIn Adult","By Fraction\nIn Prenatal",
                                                          "By Age\nIn Cytoplasm", "By Age\nIn Nucleus"))))
sigIR$ID = paste0("chr",sigIR$Chr, ":", sigIR$Start, "-", sigIR$End, ":(", sigIR$Direction, ")")
IRres = Map(cbind, IRres, intronID = lapply(IRres, function(x) paste0("chr",x$Chr, ":", x$Start, "-", x$End, ":(", x$Direction, ")")))
head(IRres[[1]])
x = lapply(IRres, function(z) data.frame(sigIR, IR = z[match(sigIR$ID, z$intronID),"IRratio"]))
x = do.call(rbind, Map(cbind, x, sampleID = as.list(names(x))))
x$genecodeID = geneMap[match(x$ensID, geneMap$ensemblID), "gencodeID"]
x$rownum = 1:nrow(x)
x$Fraction = ifelse(x$rownum %in% grep("C", x$sampleID), "Cytoplasm", "Nucleus")
x$Age = ifelse(x$rownum %in% grep("Br53", x$sampleID), "Prenatal", "Adult")
x$group = paste0(x$Age, ":", x$Fraction)
x$group = factor(x$group, levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm","Adult:Nucleus","Prenatal:Nucleus"))
head(x)

dt = data.table(x)
dt =  dt[,list(meanIR = mean(IR), medianIR = median(IR)), by=c("group","genecodeID", "comp", "Age", "Fraction")]
df = data.frame(dt)
df = cbind(rbind(df,df), LFC=c(data.frame(Apres)[match(df$genecodeID, rownames(Apres)),"log2FoldChange"],
                               data.frame(Fpres.down)[match(df$genecodeID, rownames(Fpres.down)),"log2FoldChange"]),
           padj=c(data.frame(Apres)[match(df$genecodeID, rownames(Apres)),"padj"],
                  data.frame(Fpres.down)[match(df$genecodeID, rownames(Fpres.down)),"padj"]),
           LFC.comp = c(rep.int("In Adult", nrow(df)), rep.int("In Prenatal", nrow(df))))
head(df)

df$compLocation = df$shortComp = NA
df[grep("By Fraction", df$comp),"shortComp"] = "By Fraction"
df[grep("By Age", df$comp),"shortComp"] = "By Age"
df[grep("In Adult", df$comp),"compLocation"] = "In Adult"
df[grep("In Prenatal", df$comp),"compLocation"] = "In Prenatal"
df[grep("In Cytoplasm", df$comp),"compLocation"] = "In Cytoplasm"
df[grep("In Nucleus", df$comp),"compLocation"] = "In Nucleus"

df$quad = NA
df[which(df$comp=="By Fraction\nIn Adult" & df$Age=="Adult"),"quad"] = "IR By Fraction:\nin Adult"
df[which(df$comp=="By Fraction\nIn Prenatal" & df$Age=="Prenatal"),"quad"] = "IR By Fraction:\nin Prenatal"
df[which(df$comp=="By Age\nIn Cytoplasm" & df$Fraction=="Cytoplasm"),"quad"] = "IR By Age:\nin Cytoplasm"
df[which(df$comp=="By Age\nIn Nucleus" & df$Fraction=="Nucleus"),"quad"] = "IR By Age:\nin Nucleus"
xs = df[which(df$quad %in% c("IR By Fraction:\nin Adult","IR By Age:\nin Cytoplasm","IR By Age:\nin Nucleus","IR By Fraction:\nin Prenatal")),]
head(xs)

xs$quad = gsub("IR By Fraction:\ni", "I", xs$quad)
xs$quad = gsub("IR By Age:\ni", "I", xs$quad)

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_byLFC_facet_grid.pdf", width=5,height=5)
ggplot(xs[which(xs$shortComp=="By Fraction"),], aes(x=meanIR,y=LFC)) + geom_point(aes(colour = group)) +
  facet_grid(LFC.comp ~ quad) + theme_classic() +
  xlab("mean IR Ratio") + ylab("") + ggtitle("IR By Fraction") +
  geom_hline(yintercept = 0, linetype="dashed") +
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="bottom", legend.title = element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA), nrow=2)) 
ggplot(xs[which(xs$shortComp=="By Age"),], aes(x=meanIR,y=LFC)) + geom_point(aes(colour = group)) +
  facet_grid(LFC.comp ~ quad) + theme_classic() +
  xlab("mean IR Ratio") + ylab("") + ggtitle("IR By Age") +
  geom_hline(yintercept = 0, linetype="dashed") +
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position="bottom", legend.title = element_blank()) +
  guides(color=guide_legend(override.aes=list(fill=NA), nrow=2)) 
dev.off()


res = list(fracmeanIR.inAdult.AC = cor.test(x = xs[which(xs$comp=="By Fraction\nIn Adult" & xs$group=="Adult:Cytoplasm"),"meanIR"],
                                        y = xs[which(xs$comp=="By Fraction\nIn Adult" & xs$group=="Adult:Cytoplasm"),"LFC"]),
           fracmeanIR.inAdult.AN = cor.test(x = xs[which(xs$comp=="By Fraction\nIn Adult" & xs$group=="Adult:Nucleus"),"meanIR"],
                                        y = xs[which(xs$comp=="By Fraction\nIn Adult" & xs$group=="Adult:Nucleus"),"LFC"]),
           fracmeanIR.inPrenatal.PC = cor.test(x = xs[which(xs$comp=="By Fraction\nIn Prenatal" & xs$group=="Prenatal:Cytoplasm"),"meanIR"],
                                           y = xs[which(xs$comp=="By Fraction\nIn Prenatal" & xs$group=="Prenatal:Cytoplasm"),"LFC"]),
           fracmeanIR.inPrenatal.PN = cor.test(x = xs[which(xs$comp=="By Fraction\nIn Prenatal" & xs$group=="Prenatal:Nucleus"),"meanIR"],
                                           y = xs[which(xs$comp=="By Fraction\nIn Prenatal" & xs$group=="Prenatal:Nucleus"),"LFC"]),
           agemeanIR.inCyt.AC = cor.test(x = xs[which(xs$comp=="By Age\nIn Cytoplasm" & xs$group=="Adult:Cytoplasm"),"meanIR"],
                                     y = xs[which(xs$comp=="By Age\nIn Cytoplasm" & xs$group=="Adult:Cytoplasm"),"LFC"]),
           agemeanIR.inCyt.PC = cor.test(x = xs[which(xs$comp=="By Age\nIn Cytoplasm" & xs$group=="Prenatal:Cytoplasm"),"meanIR"],
                                     y = xs[which(xs$comp=="By Age\nIn Cytoplasm" & xs$group=="Prenatal:Cytoplasm"),"LFC"]),
           agemeanIR.inNuc.AN = cor.test(x = xs[which(xs$comp=="By Age\nIn Nucleus" & xs$group=="Adult:Nucleus"),"meanIR"],
                                     y = xs[which(xs$comp=="By Age\nIn Nucleus" & xs$group=="Adult:Nucleus"),"LFC"]),
           agemeanIR.inNuc.PN = cor.test(x = xs[which(xs$comp=="By Age\nIn Nucleus" & xs$group=="Prenatal:Nucleus"),"meanIR"],
                                     y = xs[which(xs$comp=="By Age\nIn Nucleus" & xs$group=="Prenatal:Nucleus"),"LFC"]))

res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            cor = x$estimate))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#                Comparison      Tstat         pval         cor         FDR
#1    fracmeanIR.inAdult.AC  3.8808729 0.0002614566  0.44794196 0.001045826
#2    fracmeanIR.inAdult.AN  3.9865751 0.0001843804  0.45761443 0.001045826
#3 fracmeanIR.inPrenatal.PC  0.8398436 0.4253812671  0.28464634 0.567175023
#4 fracmeanIR.inPrenatal.PN  0.1925380 0.8521186335  0.06791530 0.878171752
#5       agemeanIR.inCyt.AC  0.1557583 0.8781717523  0.03891008 0.878171752
#6       agemeanIR.inCyt.PC  1.9035237 0.0751195173  0.42970572 0.150239035
#7       agemeanIR.inNuc.AN -1.1341802 0.2646525589 -0.19093195 0.423444094
#8       agemeanIR.inNuc.PN  2.3322644 0.0257447717  0.37137470 0.068652725

pdf("./Dropbox/sorted_figures/github_controlled/Intron_retention/figures/gene_IR_comparisons/dIR_byLFC_facet_grid2.pdf", width=4.18,height=4.27)
ggplot(xs[which(xs$shortComp=="By Fraction" & xs$LFC.comp=="In Adult" & xs$quad=="In Adult"),], aes(x=meanIR,y=LFC)) + geom_point(aes(colour = group)) +
  theme_classic() + xlab("mean IR Ratio") + ylab("") + ggtitle("IR By Fraction In Adult") +
  geom_smooth(aes(x=meanIR,y=LFC, color=group), method='lm', se=TRUE) + 
  geom_hline(yintercept = 0, linetype="dashed") +
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank(), legend.position = c(0.7, 0.15)) +
  guides(color=guide_legend(override.aes=list(fill=NA), nrow=2))
dev.off()


