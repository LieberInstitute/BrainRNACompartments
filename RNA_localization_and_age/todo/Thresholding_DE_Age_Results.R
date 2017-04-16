## Addressing Age differences by fraction and library ##
plotMA(Cpres, alpha=0.05, main="Cytosolic PolyA Samples", ylim=c(-8,8))
#Cpres = DE based on Age in cytosol polyA
#At FDR 0.05
Cpres.05 <- Cpres[which(Cpres$padj<=0.05),]
Cpres.05.1 <- Cpres.05[which(abs(Cpres.05$log2FoldChange) >=1),]
Cpres.05.2 <- Cpres.05[which(abs(Cpres.05$log2FoldChange) >=2),]
#At FDR 0.01
Cpres.01 <- Cpres[which(Cpres$padj<=0.01),]
Cpres.01.1 <- Cpres.01[which(abs(Cpres.01$log2FoldChange) >=1),]
Cpres.01.2 <- Cpres.01[which(abs(Cpres.01$log2FoldChange) >=2),]
#At FDR 0.001
Cpres.001 <- Cpres[which(Cpres$padj<=0.001),]
Cpres.001.1 <- Cpres.001[which(abs(Cpres.001$log2FoldChange) >=1),]
Cpres.001.2 <- Cpres.001[which(abs(Cpres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Cpres[order(rownames(Cpres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Cpres <- data.frame(dat, g)
Cpres <- Cpres[order(Cpres$padj),]
write.csv(Cpres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Cpres.csv")

#Npres = DE based on Age in nucleus polyA
plotMA(Npres, alpha=0.05, main="Nuclear PolyA Samples", ylim=c(-8,8))
#At FDR 0.05
Npres.05 <- Npres[which(Npres$padj<=0.05),]
Npres.05.1 <- Npres.05[which(abs(Npres.05$log2FoldChange) >=1),]
Npres.05.2 <- Npres.05[which(abs(Npres.05$log2FoldChange) >=2),]
#At FDR 0.01
Npres.01 <- Npres[which(Npres$padj<=0.01),]
Npres.01.1 <- Npres.01[which(abs(Npres.01$log2FoldChange) >=1),]
Npres.01.2 <- Npres.01[which(abs(Npres.01$log2FoldChange) >=2),]
#At FDR 0.001
Npres.001 <- Npres[which(Npres$padj<=0.001),]
Npres.001.1 <- Npres.001[which(abs(Npres.001$log2FoldChange) >=1),]
Npres.001.2 <- Npres.001[which(abs(Npres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Npres[order(rownames(Npres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Npres <- data.frame(dat, g)
Npres <- Npres[order(Npres$padj),]
write.csv(Npres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Npres.csv")

#Crres = DE based on Age in cytosol ribozero
plotMA(Crres, alpha=0.05, main="Cytosolic Ribozero Samples", ylim=c(-8,8))
#At FDR 0.05
Crres.05 <-   Crres[which(Crres$padj<=0.05),]
Crres.05.1 <- Crres.05[which(abs(Crres.05$log2FoldChange) >=1),]
Crres.05.2 <- Crres.05[which(abs(Crres.05$log2FoldChange) >=2),]
#At FDR 0.01
Crres.01 <-   Crres[which(Crres$padj<=0.01),]
Crres.01.1 <- Crres.01[which(abs(Crres.01$log2FoldChange) >=1),]
Crres.01.2 <- Crres.01[which(abs(Crres.01$log2FoldChange) >=2),]
#At FDR 0.001
Crres.001 <-   Crres[which(Crres$padj<=0.001),]
Crres.001.1 <- Crres.001[which(abs(Crres.001$log2FoldChange) >=1),]
Crres.001.2 <- Crres.001[which(abs(Crres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Crres[order(rownames(Crres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Crres <- data.frame(dat, g)
Crres <- Crres[order(Crres$padj),]
write.csv(Crres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Crres.csv")

#Nrres = DE based on Age in nucleus ribozero
plotMA(Nrres, alpha=0.05, main="Nuclear Ribozero Samples", ylim=c(-8,8))
#At FDR 0.05
Nrres.05 <- Nrres[which(Nrres$padj<=0.05),]
Nrres.05.1 <- Nrres.05[which(abs(Nrres.05$log2FoldChange) >=1),]
Nrres.05.2 <- Nrres.05[which(abs(Nrres.05$log2FoldChange) >=2),]
#At FDR 0.01
Nrres.01 <- Nrres[which(Nrres$padj<=0.01),]
Nrres.01.1 <- Nrres.01[which(abs(Nrres.01$log2FoldChange) >=1),]
Nrres.01.2 <- Nrres.01[which(abs(Nrres.01$log2FoldChange) >=2),]
#At FDR 0.001
Nrres.001 <- Nrres[which(Nrres$padj<=0.001),]
Nrres.001.1 <- Nrres.001[which(abs(Nrres.001$log2FoldChange) >=1),]
Nrres.001.2 <- Nrres.001[which(abs(Nrres.001$log2FoldChange) >=2),]
# Significant genes (FDR<0.05, abs(log2FoldChange>=1))
dat <- Nrres[order(rownames(Nrres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Nrres <- data.frame(dat, g)
Nrres <- Nrres[order(Nrres$padj),]
write.csv(Nrres, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nrres.csv")

### DE for Age over several LCF and FDR ###
#Upregulated Genes: LFC=1, FDR<0.05
UP.05.1 <- list("Cpres.05.1"=Cpres.05.1[which(Cpres.05.1$log2FoldChange >=1),],"Npres.05.1"=Npres.05.1[which(Npres.05.1$log2FoldChange >=1),],
                 "Crres.05.1"=Crres.05.1[which(Crres.05.1$log2FoldChange >=1),],"Nrres.05.1"=Nrres.05.1[which(Nrres.05.1$log2FoldChange >=1),])
UP.05.1 <- do.call(rbind, lapply(UP.05.1, function(x) nrow(x)))
colnames(UP.05.1) = "UP.05.1"
#Upregulated Genes: LFC=1, FDR<0.01
UP.01.1 <- list("Cpres.01.1"=Cpres.01.1[which(Cpres.01.1$log2FoldChange >=1),],"Npres.01.1"=Npres.01.1[which(Npres.01.1$log2FoldChange >=1),],
                 "Crres.01.1"=Crres.01.1[which(Crres.01.1$log2FoldChange >=1),],"Nrres.01.1"=Nrres.01.1[which(Nrres.01.1$log2FoldChange >=1),])
UP.01.1 <- do.call(rbind, lapply(UP.01.1, function(x) nrow(x)))
colnames(UP.01.1) = "UP.01.1"
#Upregulated Genes: LFC=1, FDR<0.001
UP.001.1 <- list("Cpres.001.1"=Cpres.001.1[which(Cpres.001.1$log2FoldChange >=1),],"Npres.001.1"=Npres.001.1[which(Npres.001.1$log2FoldChange >=1),],
                  "Crres.001.1"=Crres.001.1[which(Crres.001.1$log2FoldChange >=1),],"Nrres.001.1"=Nrres.001.1[which(Nrres.001.1$log2FoldChange >=1),]) 
UP.001.1 <- do.call(rbind, lapply(UP.001.1, function(x) nrow(x)))
colnames(UP.001.1) = "UP.001.1"
#Upregulated Genes: LFC=2, FDR<0.05
UP.05.2 <- list("Cpres.05.2"=Cpres.05.2[which(Cpres.05.2$log2FoldChange >=1),],"Npres.05.2"=Npres.05.2[which(Npres.05.2$log2FoldChange >=1),],
                 "Crres.05.2"=Crres.05.2[which(Crres.05.2$log2FoldChange >=1),],"Nrres.05.2"=Nrres.05.2[which(Nrres.05.2$log2FoldChange >=1),]) 
UP.05.2 <- do.call(rbind, lapply(UP.05.2, function(x) nrow(x)))
colnames(UP.05.2) = "UP.05.2"
#Upregulated Genes: LFC=2, FDR<0.01
UP.01.2 <- list("Cpres.01.2"=Cpres.01.2[which(Cpres.01.2$log2FoldChange >=1),],"Npres.01.2"=Npres.01.2[which(Npres.01.2$log2FoldChange >=1),],
                 "Crres.01.2"=Crres.01.2[which(Crres.01.2$log2FoldChange >=1),],"Nrres.01.2"=Nrres.01.2[which(Nrres.01.2$log2FoldChange >=1),])
UP.01.2 <- do.call(rbind, lapply(UP.01.2, function(x) nrow(x)))
colnames(UP.01.2) = "UP.01.2"
#Upregulated Genes: LFC=2, FDR<0.001
UP.001.2 <- list("Cpres.001.2"=Cpres.001.2[which(Cpres.001.2$log2FoldChange >=1),],"Npres.001.2"=Npres.001.2[which(Npres.001.2$log2FoldChange >=1),],
                  "Crres.001.2"=Crres.001.2[which(Crres.001.2$log2FoldChange >=1),],"Nrres.001.2"=Nrres.001.2[which(Nrres.001.2$log2FoldChange >=1),])
UP.001.2 <- do.call(rbind, lapply(UP.001.2, function(x) nrow(x)))
colnames(UP.001.2) = "UP.001.2"
#Down-regulated Genes: LFC=1, FDR<0.05
D.05.1 <- list("Cpres.05.1"=Cpres.05.1[which(Cpres.05.1$log2FoldChange <=-1),],"Npres.05.1"=Npres.05.1[which(Npres.05.1$log2FoldChange <=-1),],
                 "Crres.05.1"=Crres.05.1[which(Crres.05.1$log2FoldChange <=-1),],"Nrres.05.1"=Nrres.05.1[which(Nrres.05.1$log2FoldChange <=-1),])
D.05.1 <- do.call(rbind, lapply(D.05.1, function(x) nrow(x)))
colnames(D.05.1) = "D.05.1"
#Down-regulated Genes: LFC=1, FDR<0.01
D.01.1 <- list("Cpres.01.1"=Cpres.01.1[which(Cpres.01.1$log2FoldChange <=-1),],"Npres.01.1"=Npres.01.1[which(Npres.01.1$log2FoldChange <=-1),],
                 "Crres.01.1"=Crres.01.1[which(Crres.01.1$log2FoldChange <=-1),],"Nrres.01.1"=Nrres.01.1[which(Nrres.01.1$log2FoldChange <=-1),]) 
D.01.1 <- do.call(rbind, lapply(D.01.1, function(x) nrow(x)))
colnames(D.01.1) = "D.01.1"
#Down-regulated Genes: LFC=1, FDR<0.001
D.001.1 <- list("Cpres.001.1"=Cpres.001.1[which(Cpres.001.1$log2FoldChange <=-1),],"Npres.001.1"=Npres.001.1[which(Npres.001.1$log2FoldChange <=-1),],
                  "Crres.001.1"=Crres.001.1[which(Crres.001.1$log2FoldChange <=-1),],"Nrres.001.1"=Nrres.001.1[which(Nrres.001.1$log2FoldChange <=-1),]) 
D.001.1 <- do.call(rbind, lapply(D.001.1, function(x) nrow(x)))
colnames(D.001.1) = "D.001.1"
#Down-regulated Genes: LFC=2, FDR<0.05
D.05.2 <- list("Cpres.05.2"=Cpres.05.2[which(Cpres.05.2$log2FoldChange <=-1),],"Npres.05.2"=Npres.05.2[which(Npres.05.2$log2FoldChange <=-1),],
                 "Crres.05.2"=Crres.05.2[which(Crres.05.2$log2FoldChange <=-1),],"Nrres.05.2"=Nrres.05.2[which(Nrres.05.2$log2FoldChange <=-1),])
D.05.2 <- do.call(rbind, lapply(D.05.2, function(x) nrow(x)))
colnames(D.05.2) = "D.05.2"
#Down-regulated Genes: LFC=2, FDR<0.01
D.01.2 <- list("Cpres.01.2"=Cpres.01.2[which(Cpres.01.2$log2FoldChange <=-1),],"Npres.01.2"=Npres.01.2[which(Npres.01.2$log2FoldChange <=-1),],
                 "Crres.01.2"=Crres.01.2[which(Crres.01.2$log2FoldChange <=-1),],"Nrres.01.2"=Nrres.01.2[which(Nrres.01.2$log2FoldChange <=-1),])
D.01.2 <- do.call(rbind, lapply(D.01.2, function(x) nrow(x)))
colnames(D.01.2) = "D.01.2"
#Down-regulated Genes: LFC=2, FDR<0.001
D.001.2 <- list("Cpres.001.2"=Cpres.001.2[which(Cpres.001.2$log2FoldChange <=-1),],"Npres.001.2"=Npres.001.2[which(Npres.001.2$log2FoldChange <=-1),],
                  "Crres.001.2"=Crres.001.2[which(Crres.001.2$log2FoldChange <=-1),],"Nrres.001.2"=Nrres.001.2[which(Nrres.001.2$log2FoldChange <=-1),]) 
D.001.2 <- do.call(rbind, lapply(D.001.2, function(x) nrow(x)))
colnames(D.001.2) = "D.001.2"

readableDE = cbind(UP.05.1, UP.01.1, UP.001.1, UP.05.2, UP.01.2, UP.001.2, D.05.1, D.01.1, D.001.1, D.05.2, D.01.2, D.001.2)
rownames(readableDE) = c("Cytosol PolyA", "Nucleus PolyA", "Cytosol RiboZero", "Nucleus Ribozero")

DEbyAge = data.frame(Number = c(UP.05.1, UP.01.1, UP.001.1, UP.05.2, UP.01.2, UP.001.2, D.05.1, D.01.1, D.001.1, D.05.2, D.01.2, D.001.2),
                     Direction = factor(c("Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing",
                                          "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing", "Increasing",
                                          "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing",
                                          "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing", "Decreasing"), levels=c("Increasing", "Decreasing")),
                     FDR = factor(c("FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001",
                                    "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001",
                                    "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001",
                                    "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.05", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.01", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001", "FDR: 0.001"), levels=c("FDR: 0.05", "FDR: 0.01","FDR: 0.001")),
                     LFC = factor(c("LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2",
                                    "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 1", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2", "LFC: 2"), levels=c("LFC: 1", "LFC: 2")),
                     Group = factor(c("Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero",
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero",
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero",
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero", 
                                      "Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero"), levels=c("Cytosol\nPolyA", "Nucleus\nPolyA", "Cytosol\nRiboZero", "Nucleus\nRiboZero")))

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