library(ggplot2)
library(DESeq2)

path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "QC_section/data/rawCounts_combined_NucVSCyt_n23.rda"))
load(paste0(path, "QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda"))
load(paste0(path, "characterize_fractioned_transcriptome/data/DESeq2_results.rda"))


################################################################################
## Start with some a priori localizing genes

ens = c("ENSG00000075624","ENSG00000102081","ENSG00000229807","ENSG00000251562")
map = data.frame(ids = rownames(geneMap[which(geneMap$ensemblID %in% ens),]), 
                 ens = geneMap[which(geneMap$ensemblID %in% ens),"ensemblID"],
                 sym = geneMap[which(geneMap$ensemblID %in% ens),"Symbol"],
                 Location = c("Cytosplasmic","Nuclear","Nuclear", "Cytosplasmic"))
testgenes <- geneRpkm[which(rownames(geneRpkm) %in% map$ids),]
rownames(testgenes) <- geneMap[match(map$ids,rownames(geneMap)), "Symbol"]
testgenes <- reshape2::melt(testgenes)
testgenes$Label = pd[match(testgenes$Var2, rownames(pd)),"Label"]
testgenes$Age = pd[match(testgenes$Var2, rownames(pd)),"Fetal"]
testgenes$Fraction = pd[match(testgenes$Var2, rownames(pd)),"Zone"]
testgenes$Library = pd[match(testgenes$Var2, rownames(pd)),"Library"]
testgenes$Fraction = gsub("Cytosol", "Cytoplasm", testgenes$Fraction)
testgenes$Location = map[match(testgenes$Var1, map$sym),"Location"]
testgenes$alt = ifelse(testgenes$Location=="Nuclear", "greater","less")
testgenes <- testgenes[which(testgenes$Var2!="Br1113C1_RiboZero"),]
testgenes$Label2 <- paste0(testgenes$Age, "\n", testgenes$Fraction)
testgenes$Label2 <- factor(testgenes$Label2, levels = c("Adult\nCytoplasm",
                                                        "Prenatal\nCytoplasm",
                                                        "Adult\nNucleus",
                                                        "Prenatal\nNucleus"))
head(testgenes)

testgenes <- split(testgenes, testgenes$Var1)
testgenes2 <- lapply(testgenes, function(x) split(x, x$Label2))
testgenes <- lapply(testgenes, function(x) split(x, x$Label))
testgenes <- lapply(testgenes, function(f) lapply(f, function(x) x[order(x$Var2),]))
testgenes2 <- lapply(testgenes2, function(f) lapply(f, function(x) x[order(x$Var2),]))

r2 <- lapply(testgenes2, function(f) 
  list(Adult = t.test(log(f$`Adult\nNucleus`$value+1),
                            log(f$`Adult\nCytoplasm`$value+1), 
                            paired = TRUE, alternative = f[[1]][1,"alt"]),
       Prenatal = t.test(log(f$`Prenatal\nNucleus`$value+1),
                         log(f$`Prenatal\nCytoplasm`$value+1), 
                         paired = TRUE, alternative = f[[1]][1,"alt"])))
r2 <- do.call(rbind, Map(cbind, Gene = as.list(names(r2)),
                      lapply(r2, function(f) 
                        do.call(rbind, Map(cbind, Comparison = as.list(names(f)), 
                                           lapply(f, function(x) 
                                             data.frame(statistic = x$statistic,	
                                                        df = x$parameter,	
                                                        p.value = x$p.value,
                                                        conf.int1 = x$conf.int[1],
                                                        conf.int2 = x$conf.int[2],
                                                        alternative = x$alternative)))))))

r2
#      Gene Comparison   statistic df      p.value   conf.int1   conf.int2 alternative
#t     ACTB      Adult -10.3975863  4 0.0002415871        -Inf -0.91111535        less
#t1    ACTB   Prenatal  -4.4260550  5 0.0034266028        -Inf -0.24638690        less
#t2  MALAT1      Adult  12.8300966  4 0.0001063689  0.83466830         Inf     greater
#t11 MALAT1   Prenatal   2.8981706  5 0.0169313924  0.15985981         Inf     greater
#t3    XIST      Adult   0.9975009  4 0.1874876073 -0.16071092         Inf     greater
#t12   XIST   Prenatal   1.6090538  5 0.0842590310 -0.06717991         Inf     greater
#t4    FMR1      Adult  -3.5373985  4 0.0120348429        -Inf -0.08925008        less
#t13   FMR1   Prenatal  -0.6162072  5 0.2823633376        -Inf  0.04516074        less


res <- lapply(testgenes, function(f) 
  list(Adult.PolyA = t.test(log(f$`Adult\nNucleus\npolyA`$value+1),
                            log(f$`Adult\nCytosol\npolyA`$value+1), 
                            paired = TRUE, alternative = f[[1]][1,"alt"]),
       Prenatal.PolyA = t.test(log(f$`Prenatal\nNucleus\npolyA`$value+1),
                               log(f$`Prenatal\nCytosol\npolyA`$value+1), 
                               paired = TRUE, alternative = f[[1]][1,"alt"]),
       Adult.RiboZero = t.test(log(f$`Adult\nNucleus\nRiboZero`$value+1),
                               log(f$`Adult\nCytosol\nRiboZero`$value+1),
                               paired = TRUE, alternative = f[[1]][1,"alt"]),
       Prenatal.RiboZero = t.test(log(f$`Prenatal\nNucleus\nRiboZero`$value+1),
                                  log(f$`Prenatal\nCytosol\nRiboZero`$value+1), 
                                  paired = TRUE, alternative = f[[1]][1,"alt"])))

res <- do.call(rbind, Map(cbind, Gene = as.list(names(res)),
                      lapply(res, function(f) 
                       do.call(rbind, Map(cbind, Comparison = as.list(names(f)), 
                                          lapply(f, function(x) 
                                            data.frame(statistic = x$statistic,	
                                            df = x$parameter,	
                                            p.value = x$p.value,
                                            conf.int1 = x$conf.int[1],
                                            conf.int2 = x$conf.int[2],
                                            alternative = x$alternative)))))))

write.csv(res, quote=FALSE, 
          file = paste0(path, "QC_section/data/known_localizing_genes.csv"))

res
res[res$p.value>0.05,]
head(testgenes)


# Plot the Log2 Fold Change and SE by library for these genes

FracList <- list(Apres= Apres, Fpres = Fpres.down, 
                 Arres = Arres, Frres = Frres)
FracList <- lapply(FracList, function(y) data.frame(y))

lfc <- lapply(FracList, function(x) 
  x[which(rownames(x) %in% map[which(map$sym=="ACTB" | map$sym=="MALAT1"),"ids"]),])
lfc <- do.call(rbind, Map(cbind, Group = list("Adult: poly(A)", "Prenatal: poly(A)",
                                              "Adult: Ribo-Zero","Prenatal: Ribo-Zero"),
                          Symbol = list(c("ACTB","MALAT1"),c("ACTB","MALAT1"),
                                        c("ACTB","MALAT1"),c("ACTB","MALAT1")), lfc))
lfc[grep("poly",lfc$Group),"Library"] <- "poly(A)"
lfc[grep("Ribo",lfc$Group),"Library"] <- "Ribo-Zero"
lfc[grep("Adult",lfc$Group),"Age"] <- "Adult"
lfc[-grep("Adult",lfc$Group),"Age"] <- "Prenatal"

head(lfc)

pdf(paste0(path, "QC_section/data/known_localizing_genes_ACTB_MALAT1_LFC.pdf"), 
    width = 4.5, height = 4.5)

dodge <- position_dodge(width=0.9)
limits <- aes(ymax = log2FoldChange + lfcSE, ymin=log2FoldChange - lfcSE)

ggplot(lfc, aes(x=Symbol, y=log2FoldChange, fill=Age)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(mapping = limits, position = dodge, width=0.25) +
  ylim(-2,2) +
  facet_grid(. ~ Library) +
  scale_fill_brewer(palette="Set1") +
  ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  xlab("") +
  theme(title = element_text(size = 18),
        text = element_text(size = 18),
        legend.position = "bottom",
        legend.title = element_blank())

dev.off()


################################################################################
## Expand to the top differentially expressed fraction genes in the ENCODE samples

load(paste0(path, "other_datasets/data/ENCODE_DESeq2_output.rda"))

frac <- sigres$Fraction
frac <- lapply(frac, function(x) x[which(rownames(x) %in% rownames(Ires.down)),])
frac <- lapply(frac, function(x) x[which(abs(x$log2FoldChange)>=3),])
elementNROWS(frac)
frac$Cytoplasmic <- frac$Cytoplasmic[order(frac$Cytoplasmic$padj),]
frac$Cytoplasmic <- data.frame(ids = rownames(frac$Cytoplasmic[1:100,]),
                               Tstat = frac$Cytoplasmic$stat[1:100])
frac$Nuclear <- frac$Nuclear[which(abs(frac$Nuclear$log2FoldChange)>=5.5),]
frac$Nuclear <- frac$Nuclear[order(frac$Nuclear$padj),]
frac$Nuclear <- data.frame(ids = rownames(frac$Nuclear), 
                           Tstat = frac$Nuclear$stat)
frac <- do.call(rbind, Map(cbind, Location = as.list(names(frac)), frac))

map <- cbind(frac, sym = geneMap[match(frac$ids, rownames(geneMap)),"Symbol"])

testgenes <- geneRpkm[which(rownames(geneRpkm) %in% map$ids),]
testgenes <- reshape2::melt(testgenes)
testgenes$Label = pd[match(testgenes$Var2, rownames(pd)),"Label"]
testgenes$Age = pd[match(testgenes$Var2, rownames(pd)),"Fetal"]
testgenes$Fraction = pd[match(testgenes$Var2, rownames(pd)),"Zone"]
testgenes$Library = pd[match(testgenes$Var2, rownames(pd)),"Library"]
testgenes$Fraction = gsub("Cytosol", "Cytoplasm", testgenes$Fraction)
testgenes$Location = map[match(testgenes$Var1, map$ids),"Location"]
testgenes$alt = ifelse(testgenes$Location=="Nuclear", "greater","less")
testgenes <- testgenes[which(testgenes$Var2!="Br1113C1_RiboZero"),]

head(testgenes)
testgenes <- split(testgenes, testgenes$Var1)
testgenes <- lapply(testgenes, function(x) split(x, x$Label))
testgenes <- lapply(testgenes, function(f) lapply(f, function(x) x[order(x$Var2),]))

# paired one-tailed
res <- lapply(testgenes, function(f) 
  list(Adult.PolyA = t.test(log(f$`Adult\nNucleus\npolyA`$value+1),
                            log(f$`Adult\nCytosol\npolyA`$value+1), 
                            paired = TRUE, alternative = f[[1]][1,"alt"]),
       Prenatal.PolyA = t.test(log(f$`Prenatal\nNucleus\npolyA`$value+1),
                               log(f$`Prenatal\nCytosol\npolyA`$value+1), 
                               paired = TRUE, alternative = f[[1]][1,"alt"]),
       Adult.RiboZero = t.test(log(f$`Adult\nNucleus\nRiboZero`$value+1),
                               log(f$`Adult\nCytosol\nRiboZero`$value+1),
                               paired = TRUE, alternative = f[[1]][1,"alt"]),
       Prenatal.RiboZero = t.test(log(f$`Prenatal\nNucleus\nRiboZero`$value+1),
                                  log(f$`Prenatal\nCytosol\nRiboZero`$value+1), 
                                  paired = TRUE, alternative = f[[1]][1,"alt"])))

# paired 2 tailed
res <- lapply(testgenes, function(f) 
  list(Adult.PolyA = t.test(log(f$`Adult\nNucleus\npolyA`$value+1),
                            log(f$`Adult\nCytosol\npolyA`$value+1), 
                            paired = TRUE),
       Prenatal.PolyA = t.test(log(f$`Prenatal\nNucleus\npolyA`$value+1),
                               log(f$`Prenatal\nCytosol\npolyA`$value+1), 
                               paired = TRUE),
       Adult.RiboZero = t.test(log(f$`Adult\nNucleus\nRiboZero`$value+1),
                               log(f$`Adult\nCytosol\nRiboZero`$value+1),
                               paired = TRUE),
       Prenatal.RiboZero = t.test(log(f$`Prenatal\nNucleus\nRiboZero`$value+1),
                                  log(f$`Prenatal\nCytosol\nRiboZero`$value+1), 
                                  paired = TRUE)))

res <- do.call(rbind, Map(cbind, Gene = as.list(names(res)),
                          lapply(res, function(f) 
                            do.call(rbind, Map(cbind, Comparison = as.list(names(f)), 
                                               lapply(f, function(x) 
                                                 data.frame(statistic = x$statistic,	
                                                            df = x$parameter,	
                                                            p.value = x$p.value,
                                                            conf.int1 = x$conf.int[1],
                                                            conf.int2 = x$conf.int[2],
                                                            alternative = x$alternative)
                                                 ))))))

res$ENCstat = map[match(res$Gene, map$ids),"Tstat"]
head(res)

write.csv(res, quote=FALSE, 
          file = paste0(path, "QC_section/data/ENCODE_topfractiongenes.csv"))


pdf(paste0("./Dropbox/sorted_figures/github_controlled/QC_section/",
           "data/ENCODE_topfractiongenes_tstats_inBrain.pdf"),
    width = 3.75, height = 4)
ggplot(res[res$Comparison=="Adult.PolyA",], aes(x = ENCstat, y = statistic)) + 
  geom_point() + theme_bw() +
  xlim(-25,25) + ylim(-25,25) +
  labs(fill="") + ylab("Brain") + xlab("ENCODE") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  ggtitle("Adult: poly(A)") +
  theme(title = element_text(size = 20),
        text = element_text(size = 20))

ggplot(res[res$Comparison=="Prenatal.PolyA",], aes(x = ENCstat, y = statistic)) + 
  geom_point() + theme_bw() +
  xlim(-25,25) + ylim(-25,25) +
  labs(fill="") + ylab("Brain") + xlab("ENCODE") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  ggtitle("Prenatal: poly(A)") +
  theme(title = element_text(size = 20),
        text = element_text(size = 20))

ggplot(res[res$Comparison=="Adult.RiboZero",], aes(x = ENCstat, y = statistic)) + 
  geom_point() + theme_bw() +
  xlim(-25,25) + ylim(-25,25) +
  labs(fill="") + ylab("Brain") + xlab("ENCODE") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  ggtitle("Adult: Ribo-Zero") +
  theme(title = element_text(size = 20),
        text = element_text(size = 20))

ggplot(res[res$Comparison=="Prenatal.RiboZero",], aes(x = ENCstat, y = statistic)) + 
  geom_point() + theme_bw() +
  xlim(-25,25) + ylim(-25,25) +
  labs(fill="") + ylab("Brain") + xlab("ENCODE") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  ggtitle("Prenatal: Ribo-Zero") +
  theme(title = element_text(size = 20),
        text = element_text(size = 20))

dev.off()


res <- split(res, res$Comparison)
df <- lapply(res, function(x) fisher.test(table(x$statistic>0, x$ENCstat>0)))
df <- do.call(rbind, Map(cbind, Comparison = as.list(names(df)),
                         lapply(df, function(x) data.frame(OR = x$estimate,
                                                           pval = x$p.value))))
df
#        Comparison        OR         pval
#       Adult.PolyA 24.578987 9.415994e-21
#    Prenatal.PolyA  5.256085 5.960059e-08
#    Adult.RiboZero 45.471876 3.169586e-34
# Prenatal.RiboZero 46.707170 2.887643e-36