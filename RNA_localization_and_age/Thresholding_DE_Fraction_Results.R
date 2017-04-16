library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

#Apres = DE based on Fraction in adult polyA
plotMA(Apres, alpha=0.05, main="Adult PolyA Samples", ylim=c(-3,3))
#Fpres = DE based on Fraction in Prenatal polyA
plotMA(Fpres, alpha=0.05, main="Prenatal PolyA Samples", ylim=c(-3,3))
#Fpres.down = DE based on Fraction in Prenatal polyA
plotMA(Fpres.down, alpha=0.05, main="Prenatal PolyA Samples", ylim=c(-3,3))
# Arres = DE based on fraction in adult RiboZero
plotMA(Arres, alpha=0.05, main="Adult RiboZero Samples", ylim=c(-3,3))
# Frres = DE based on fraction in Prenatal RiboZero
plotMA(Frres, alpha=0.05, main="Prenatal RiboZero Samples", ylim=c(-3,3))

### DE for Fraction over several LCF and FDR ###
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres), 
                Arres = data.frame(Arres), Frres = data.frame(Frres),
                Fpres.down = data.frame(Fpres.down))
#Upregulated Genes: LFC=1, FDR<0.05
list = list(UP.05.1 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange >=1),]))),
            UP.01.1 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange >=1),]))),
            UP.001.1 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange >=1),]))),
            UP.05.2 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange >=2),]))),
            UP.01.2 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange >=2),]))),
            UP.001.2 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange >=2),]))),
            D.05.1 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange <=-1),]))),
            D.01.1 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange <=-1),]))),
            D.001.1 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange <=-1),]))),
            D.05.2 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange <=-2),]))),
            D.01.2 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange <=-2),]))),
            D.001.2 = do.call(rbind, lapply(FracList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange <=-2),]))))

DEbyFraction = data.frame(Number = do.call(rbind,lapply(list, function(x) data.frame(Number=x[1:4,]))),
                          Direction = factor(c(rep.int("Nuclear",24), rep.int("Cytosolic",24)), levels=c("Nuclear", "Cytosolic")),
                          FDR = factor(c(rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4)), levels=c("FDR: 0.05", "FDR: 0.01","FDR: 0.001")),
                          LFC = factor(c(rep.int("LFC: 1",12), rep.int("LFC: 2",12), rep.int("LFC: 1",12), rep.int("LFC: 2",12)), levels=c("LFC: 1", "LFC: 2")),
                          Group = factor(rep.int(c("Adult\nPolyA", "Prenatal\nPolyA", "Adult\nRiboZero", "Prenatal\nRiboZero"),12), 
                                         levels=c("Adult\nPolyA", "Prenatal\nPolyA", "Adult\nRiboZero", "Prenatal\nRiboZero")))

ggplot(DEbyFraction, aes(x=Group, y=Number, fill=Direction), color=Direction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  facet_grid(LFC ~ FDR) +
  ylab("Gene Count") + 
  xlab("") +
  ggtitle("Number of Genes Enriched Nuclear or Cytosolic RNA") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

DEbyFraction.down = data.frame(Number = do.call(rbind,lapply(list, function(x) data.frame(Number=x[c(1,5,3,4),]))),
                          Direction = factor(c(rep.int("Nuclear",24), rep.int("Cytosolic",24)), levels=c("Nuclear", "Cytosolic")),
                          FDR = factor(c(rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4)), levels=c("FDR: 0.05", "FDR: 0.01","FDR: 0.001")),
                          LFC = factor(c(rep.int("LFC: 1",12), rep.int("LFC: 2",12), rep.int("LFC: 1",12), rep.int("LFC: 2",12)), levels=c("LFC: 1", "LFC: 2")),
                          Group = factor(rep.int(c("Adult\nPolyA", "Prenatal\nPolyA", "Adult\nRiboZero", "Prenatal\nRiboZero"),12), 
                                         levels=c("Adult\nPolyA", "Prenatal\nPolyA", "Adult\nRiboZero", "Prenatal\nRiboZero")))

ggplot(DEbyFraction.down, aes(x=Group, y=Number, fill=Direction), color=Direction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  facet_grid(LFC ~ FDR) +
  ylab("Gene Count") + 
  xlab("") +
  ggtitle("Number of Genes Enriched Nuclear or Cytosolic RNA") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))