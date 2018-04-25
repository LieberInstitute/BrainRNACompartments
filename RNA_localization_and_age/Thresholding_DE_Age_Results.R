library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

#Cpres = DE based on Age in Cytoplasm polyA
plotMA(Cpres, alpha=0.05, main="Cytoplasm PolyA Samples", ylim=c(-10,10))
#Npres = DE based on Age in Nucleus polyA
plotMA(Npres, alpha=0.05, main="Nucleus PolyA Samples", ylim=c(-10,10))
#Cpres.down = DE based on Age in Nucleus polyA
plotMA(Cpres.down, alpha=0.05, main="Cytoplasm PolyA Samples", ylim=c(-10,10))
# Crres = DE based on Age in Cytoplasm RiboZero
plotMA(Crres, alpha=0.05, main="Cytoplasm RiboZero Samples", ylim=c(-10,10))
# Nrres = DE based on Age in Nucleus RiboZero
plotMA(Nrres, alpha=0.05, main="Nucleus RiboZero Samples", ylim=c(-10,10))

### DE for Age over several LCF and FDR ###
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres), 
                Crres = data.frame(Crres), Nrres = data.frame(Nrres),
                Cpres.down = data.frame(Cpres.down))
#Upregulated Genes: LFC=1, FDR<0.05
list = list(UP.05.1 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange >=1),]))),
            UP.01.1 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange >=1),]))),
            UP.001.1 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange >=1),]))),
            UP.05.2 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange >=2),]))),
            UP.01.2 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange >=2),]))),
            UP.001.2 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange >=2),]))),
            D.05.1 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange <=-1),]))),
            D.01.1 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange <=-1),]))),
            D.001.1 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange <=-1),]))),
            D.05.2 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.05 & x$log2FoldChange <=-2),]))),
            D.01.2 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.01 & x$log2FoldChange <=-2),]))),
            D.001.2 = do.call(rbind, lapply(AgeList, function(x) nrow(x[which(x$padj<=0.001 & x$log2FoldChange <=-2),]))))

DEbyAge = data.frame(Number = do.call(rbind,lapply(list, function(x) data.frame(Number=x[1:4,]))),
                          Direction = factor(c(rep.int("Decreasing",24), rep.int("Increasing",24)), levels=c("Decreasing", "Increasing")),
                          FDR = factor(c(rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                         rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4)), 
                                       levels=c("FDR: 0.05", "FDR: 0.01","FDR: 0.001")),
                          LFC = factor(c(rep.int("LFC: 1",12), rep.int("LFC: 2",12), rep.int("LFC: 1",12), rep.int("LFC: 2",12)), 
                                       levels=c("LFC: 1", "LFC: 2")),
                          Group = factor(rep.int(c("Cytoplasm\nPolyA", "Nucleus\nPolyA", "Cytoplasm\nRiboZero", "Nucleus\nRiboZero"),12), 
                                         levels=c("Cytoplasm\nPolyA", "Nucleus\nPolyA", "Cytoplasm\nRiboZero", "Nucleus\nRiboZero")))

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/age_DEG_different_LFC&FDR_downsampled.pdf", height = 6.5, width = 17)
ggplot(DEbyAge, aes(x=Group, y=Number, fill=Direction), color=Direction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  facet_grid(LFC ~ FDR) +
  ylab("Gene Count") + 
  xlab("") +
  ggtitle("Number of Genes Increasing or Decreasing Over Development") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

DEbyAge.down = data.frame(Number = do.call(rbind,lapply(list, function(x) data.frame(Number=x[c(5,2:4),]))),
                     Direction = factor(c(rep.int("Decreasing",24), rep.int("Increasing",24)), levels=c("Decreasing", "Increasing")),
                     FDR = factor(c(rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                    rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                    rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4),
                                    rep.int("FDR: 0.05",4), rep.int("FDR: 0.01",4), rep.int("FDR: 0.001",4)), 
                                  levels=c("FDR: 0.05", "FDR: 0.01","FDR: 0.001")),
                     LFC = factor(c(rep.int("LFC: 1",12), rep.int("LFC: 2",12), rep.int("LFC: 1",12), rep.int("LFC: 2",12)), 
                                  levels=c("LFC: 1", "LFC: 2")),
                     Group = factor(rep.int(c("Cytoplasm\nPolyA", "Nucleus\nPolyA", "Cytoplasm\nRiboZero", "Nucleus\nRiboZero"),12), 
                                    levels=c("Cytoplasm\nPolyA", "Nucleus\nPolyA", "Cytoplasm\nRiboZero", "Nucleus\nRiboZero")))

ggplot(DEbyAge.down, aes(x=Group, y=Number, fill=Direction), color=Direction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  facet_grid(LFC ~ FDR) +
  ylab("Gene Count") + 
  xlab("") +
  ggtitle("Number of Genes Increasing or Decreasing Over Development") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))