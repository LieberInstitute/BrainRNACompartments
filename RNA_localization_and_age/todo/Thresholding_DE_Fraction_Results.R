## Addressing cytosolic vs. nuclear differences by age and library ##
#Apres = DE based on Fraction in adult polyA
plotMA(Apres, alpha=0.05, main="Adult PolyA Samples", ylim=c(-3,3))
#At FDR 0.05
Apres.05 <- Apres[which(Apres$padj<=0.05),]
Apres.05.1 <- Apres.05[which(abs(Apres.05$log2FoldChange) >=1),]
Apres.05.2 <- Apres.05[which(abs(Apres.05$log2FoldChange) >=2),]
#At FDR 0.01
Apres.01 <- Apres[which(Apres$padj<=0.01),]
Apres.01.1 <- Apres.01[which(abs(Apres.01$log2FoldChange) >=1),]
Apres.01.2 <- Apres.01[which(abs(Apres.01$log2FoldChange) >=2),]
#At FDR 0.001
Apres.001 <- Apres[which(Apres$padj<=0.001),]
Apres.001.1 <- Apres.001[which(abs(Apres.001$log2FoldChange) >=1),]
Apres.001.2 <- Apres.001[which(abs(Apres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Apres[order(rownames(Apres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Apres <- data.frame(dat, g)
Apres <- Apres[order(Apres$padj),]
write.csv(Apres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv")

#Fpres = DE based on Fraction in fetal polyA
plotMA(Fpres, alpha=0.05, main="Fetal PolyA Samples", ylim=c(-3,3))
#At FDR 0.05
Fpres.05 <- Fpres[which(Fpres$padj<=0.05),]
Fpres.05.1 <- Fpres.05[which(abs(Fpres.05$log2FoldChange) >=1),]
Fpres.05.2 <- Fpres.05[which(abs(Fpres.05$log2FoldChange) >=2),]
#At FDR 0.01
Fpres.01 <- Fpres[which(Fpres$padj<=0.01),]
Fpres.01.1 <- Fpres.01[which(abs(Fpres.01$log2FoldChange) >=1),]
Fpres.01.2 <- Fpres.01[which(abs(Fpres.01$log2FoldChange) >=2),]
#At FDR 0.001
Fpres.001 <- Fpres[which(Fpres$padj<=0.001),]
Fpres.001.1 <- Fpres.001[which(abs(Fpres.001$log2FoldChange) >=1),]
Fpres.001.2 <- Fpres.001[which(abs(Fpres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Fpres[order(rownames(Fpres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Fpres <- data.frame(dat, g)
Fpres <- Fpres[order(Fpres$padj),]
write.csv(Fpres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv")

# Arres = DE based on fraction in adult ribozero
plotMA(Arres, alpha=0.05, main="Adult Ribozero Samples", ylim=c(-3,3))
#At FDR 0.05
Arres.05 <- Arres[which(Arres$padj<=0.05),]
Arres.05.1 <- Arres.05[which(abs(Arres.05$log2FoldChange) >=1),]
Arres.05.2 <- Arres.05[which(abs(Arres.05$log2FoldChange) >=2),]
#At FDR 0.01
Arres.01 <- Arres[which(Arres$padj<=0.01),]
Arres.01.1 <- Arres.01[which(abs(Arres.01$log2FoldChange) >=1),]
Arres.01.2 <- Arres.01[which(abs(Arres.01$log2FoldChange) >=2),]
#At FDR 0.001
Arres.001 <- Arres[which(Arres$padj<=0.001),]
Arres.001.1 <- Arres.001[which(abs(Arres.001$log2FoldChange) >=1),]
Arres.001.2 <- Arres.001[which(abs(Arres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Arres[order(rownames(Arres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Arres <- data.frame(dat, g)
Arres <- Arres[order(Arres$padj),]
write.csv(Arres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv")

# Frres = DE based on fraction in fetal ribozero
plotMA(Frres, alpha=0.05, main="Fetal Ribozero Samples", ylim=c(-3,3))
#At FDR 0.05
Frres.05 <- Frres[which(Frres$padj<=0.05),]
Frres.05.1 <- Frres.05[which(abs(Frres.05$log2FoldChange) >=1),]
Frres.05.2 <- Frres.05[which(abs(Frres.05$log2FoldChange) >=2),]
#At FDR 0.01
Frres.01 <- Frres[which(Frres$padj<=0.01),]
Frres.01.1 <- Frres.01[which(abs(Frres.01$log2FoldChange) >=1),]
Frres.01.2 <- Frres.01[which(abs(Frres.01$log2FoldChange) >=2),]
#At FDR 0.001
Frres.001 <- Frres[which(Frres$padj<=0.001),]
Frres.001.1 <- Frres.001[which(abs(Frres.001$log2FoldChange) >=1),]
Frres.001.2 <- Frres.001[which(abs(Frres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Frres[order(rownames(Frres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Frres <- data.frame(dat, g)
Frres <- Frres[order(Frres$padj),]
write.csv(Frres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv")

# Lnres = DE based on library in the nucleus
plotMA(Lnres, alpha=0.05, main="Library Differences in the Nucleus", ylim=c(-5,10))
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Lnres[order(rownames(Lnres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Lnres <- data.frame(dat, g)
Lnres <- Lnres[order(Lnres$padj),]
write.csv(Lnres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Lnres.csv")

### DE for Fraction over several LCF and FDR ###
#Upregulated Genes: LFC=1, FDR<0.05
Res.05.1 <- list("Apres.05.1"=Apres.05.1[which(Apres.05.1$log2FoldChange >=1),],"Fpres.05.1"=Fpres.05.1[which(Fpres.05.1$log2FoldChange >=1),],
                 "Arres.05.1"=Arres.05.1[which(Arres.05.1$log2FoldChange >=1),],"Frres.05.1"=Frres.05.1[which(Frres.05.1$log2FoldChange >=1),])
UP.05.1 <- do.call(rbind, lapply(Res.05.1, function(x) nrow(x)))
colnames(UP.05.1) = "UP.05.1"
#Upregulated Genes: LFC=1, FDR<0.01
Res.01.1 <- list("Apres.01.1"=Apres.01.1[which(Apres.01.1$log2FoldChange >=1),],"Fpres.01.1"=Fpres.01.1[which(Fpres.01.1$log2FoldChange >=1),],
                 "Arres.01.1"=Arres.01.1[which(Arres.01.1$log2FoldChange >=1),],"Frres.01.1"=Frres.01.1[which(Frres.01.1$log2FoldChange >=1),])
UP.01.1 <- do.call(rbind, lapply(Res.01.1, function(x) nrow(x)))
colnames(UP.01.1) = "UP.01.1"
#Upregulated Genes: LFC=1, FDR<0.001
Res.001.1 <- list("Apres.001.1"=Apres.001.1[which(Apres.001.1$log2FoldChange >=1),],"Fpres.001.1"=Fpres.001.1[which(Fpres.001.1$log2FoldChange >=1),],
                  "Arres.001.1"=Arres.001.1[which(Arres.001.1$log2FoldChange >=1),],"Frres.001.1"=Frres.001.1[which(Frres.001.1$log2FoldChange >=1),]) 
UP.001.1 <- do.call(rbind, lapply(Res.001.1, function(x) nrow(x)))
colnames(UP.001.1) = "UP.001.1"
#Upregulated Genes: LFC=2, FDR<0.05
Res.05.2 <- list("Apres.05.2"=Apres.05.2[which(Apres.05.2$log2FoldChange >=1),],"Fpres.05.2"=Fpres.05.2[which(Fpres.05.2$log2FoldChange >=1),],
                 "Arres.05.2"=Arres.05.2[which(Arres.05.2$log2FoldChange >=1),],"Frres.05.2"=Frres.05.2[which(Frres.05.2$log2FoldChange >=1),]) 
UP.05.2 <- do.call(rbind, lapply(Res.05.2, function(x) nrow(x)))
colnames(UP.05.2) = "UP.05.2"
#Upregulated Genes: LFC=2, FDR<0.01
Res.01.2 <- list("Apres.01.2"=Apres.01.2[which(Apres.01.2$log2FoldChange >=1),],"Fpres.01.2"=Fpres.01.2[which(Fpres.01.2$log2FoldChange >=1),],
                 "Arres.01.2"=Arres.01.2[which(Arres.01.2$log2FoldChange >=1),],"Frres.01.2"=Frres.01.2[which(Frres.01.2$log2FoldChange >=1),])
UP.01.2 <- do.call(rbind, lapply(Res.01.2, function(x) nrow(x)))
colnames(UP.01.2) = "UP.01.2"
#Upregulated Genes: LFC=2, FDR<0.001
Res.001.2 <- list("Apres.001.2"=Apres.001.2[which(Apres.001.2$log2FoldChange >=1),],"Fpres.001.2"=Fpres.001.2[which(Fpres.001.2$log2FoldChange >=1),],
                  "Arres.001.2"=Arres.001.2[which(Arres.001.2$log2FoldChange >=1),],"Frres.001.2"=Frres.001.2[which(Frres.001.2$log2FoldChange >=1),])
UP.001.2 <- do.call(rbind, lapply(Res.001.2, function(x) nrow(x)))
colnames(UP.001.2) = "UP.001.2"
#Down-regulated Genes: LFC=1, FDR<0.05
Res.05.1 <- list("Apres.05.1"=Apres.05.1[which(Apres.05.1$log2FoldChange <=-1),],"Fpres.05.1"=Fpres.05.1[which(Fpres.05.1$log2FoldChange <=-1),],
                 "Arres.05.1"=Arres.05.1[which(Arres.05.1$log2FoldChange <=-1),],"Frres.05.1"=Frres.05.1[which(Frres.05.1$log2FoldChange <=-1),])
D.05.1 <- do.call(rbind, lapply(Res.05.1, function(x) nrow(x)))
colnames(D.05.1) = "D.05.1"
#Down-regulated Genes: LFC=1, FDR<0.01
Res.01.1 <- list("Apres.01.1"=Apres.01.1[which(Apres.01.1$log2FoldChange <=-1),],"Fpres.01.1"=Fpres.01.1[which(Fpres.01.1$log2FoldChange <=-1),],
                 "Arres.01.1"=Arres.01.1[which(Arres.01.1$log2FoldChange <=-1),],"Frres.01.1"=Frres.01.1[which(Frres.01.1$log2FoldChange <=-1),]) 
D.01.1 <- do.call(rbind, lapply(Res.01.1, function(x) nrow(x)))
colnames(D.01.1) = "D.01.1"
#Down-regulated Genes: LFC=1, FDR<0.001
Res.001.1 <- list("Apres.001.1"=Apres.001.1[which(Apres.001.1$log2FoldChange <=-1),],"Fpres.001.1"=Fpres.001.1[which(Fpres.001.1$log2FoldChange <=-1),],
                  "Arres.001.1"=Arres.001.1[which(Arres.001.1$log2FoldChange <=-1),],"Frres.001.1"=Frres.001.1[which(Frres.001.1$log2FoldChange <=-1),]) 
D.001.1 <- do.call(rbind, lapply(Res.001.1, function(x) nrow(x)))
colnames(D.001.1) = "D.001.1"
#Down-regulated Genes: LFC=2, FDR<0.05
Res.05.2 <- list("Apres.05.2"=Apres.05.2[which(Apres.05.2$log2FoldChange <=-1),],"Fpres.05.2"=Fpres.05.2[which(Fpres.05.2$log2FoldChange <=-1),],
                 "Arres.05.2"=Arres.05.2[which(Arres.05.2$log2FoldChange <=-1),],"Frres.05.2"=Frres.05.2[which(Frres.05.2$log2FoldChange <=-1),])
D.05.2 <- do.call(rbind, lapply(Res.05.2, function(x) nrow(x)))
colnames(D.05.2) = "D.05.2"
#Down-regulated Genes: LFC=2, FDR<0.01
Res.01.2 <- list("Apres.01.2"=Apres.01.2[which(Apres.01.2$log2FoldChange <=-1),],"Fpres.01.2"=Fpres.01.2[which(Fpres.01.2$log2FoldChange <=-1),],
                 "Arres.01.2"=Arres.01.2[which(Arres.01.2$log2FoldChange <=-1),],"Frres.01.2"=Frres.01.2[which(Frres.01.2$log2FoldChange <=-1),])
D.01.2 <- do.call(rbind, lapply(Res.01.2, function(x) nrow(x)))
colnames(D.01.2) = "D.01.2"
#Down-regulated Genes: LFC=2, FDR<0.001
Res.001.2 <- list("Apres.001.2"=Apres.001.2[which(Apres.001.2$log2FoldChange <=-1),],"Fpres.001.2"=Fpres.001.2[which(Fpres.001.2$log2FoldChange <=-1),],
                  "Arres.001.2"=Arres.001.2[which(Arres.001.2$log2FoldChange <=-1),],"Frres.001.2"=Frres.001.2[which(Frres.001.2$log2FoldChange <=-1),]) 
D.001.2 <- do.call(rbind, lapply(Res.001.2, function(x) nrow(x)))
colnames(D.001.2) = "D.001.2"

readableDE = cbind(UP.05.1, UP.01.1, UP.001.1, UP.05.2, UP.01.2, UP.001.2, D.05.1, D.01.1, D.001.1, D.05.2, D.01.2, D.001.2)
rownames(readableDE) = c("Adult PolyA", "Fetal PolyA", "Adult RiboZero", "Fetal Ribozero")

DEbyFraction = data.frame(Number = c(UP.05.1, UP.01.1, UP.001.1, UP.05.2, UP.01.2, UP.001.2, D.05.1, D.01.1, D.001.1, D.05.2, D.01.2, D.001.2),
                          Direction = factor(c("Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched",
                                               "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched", "Enriched",
                                               "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted",
                                               "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted", "Depleted"), levels=c("Enriched", "Depleted")),
                          FDR = factor(c("FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001",
                                         "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001",
                                         "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001",
                                         "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001"), levels=c("FDR: 0.05", "FDR: 0.01","FDR: 0.001")),
                          LFC = factor(c("LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2",
                                         "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2"), levels=c("LFC: 1", "LFC: 2")),
                          Group = factor(c("Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero",
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero",
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero",
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero", 
                                           "Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero"), levels=c("Adult\nPolyA", "Fetal\nPolyA", "Adult\nRiboZero", "Fetal\nRiboZero")))

ggplot(DEbyFraction, aes(x=Group, y=Number, fill=Direction), color=Direction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  facet_grid(LFC ~ FDR) +
  ylab("Gene Count") + 
  xlab("") +
  ggtitle("Number of Genes Enriched or Depleted in Nuclear RNA") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))