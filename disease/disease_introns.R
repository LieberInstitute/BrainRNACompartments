library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)


path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "characterize_fractioned_transcriptome/data/DESeq2_results.rda"))
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))
load(paste0(path, "updated_gene_sets.rda"))

# Load IRFinder Results

names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],
                                 "/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames
IRres = Map(cbind, IRres, genes = lapply(lapply(IRres, function(x) 
  unlist(strsplit(as.character(x$GeneIntronDetails), "/", 
                  fixed = TRUE), recursive = F)), 
  function(x) x[grep("ENSG", x)]), 
  intronID = lapply(IRres, function(x) 
    paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
elementNROWS(IRres)
head(IRres[[1]])
introns = do.call(rbind, Map(cbind, 
                             lapply(IRres, function(x) 
                               data.frame(x[,which(colnames(x) %in% 
                                         c("Chr","Start","End","Direction",
                                           "IRratio","Warnings","genes",
                                           "intronID"))])), 
                             SampleID = as.list(names(IRres))))
introns$Fraction = ifelse(introns$SampleID %in% 
                            c("Br1113N1","Br2046N","Br2074N","Br5339N1",
                              "Br5340N1","Br5341N1"), "Nucleus","Cytoplasm")
introns$Age = ifelse(introns$SampleID %in% c("Br5339C1","Br5339N1","Br5340C1",
                                             "Br5340N1","Br5341C1","Br5341N1"), 
                     "Prenatal","Adult")
introns$SampleGroup = paste(introns$Age, introns$Fraction, sep=":")
head(introns)


## Get disease gene introns

geneuniverse <- na.omit(unique(as.character(
  geneMap[which(geneMap$gencodeID %in% 
                  rownames(Ipres.down)),"gencodeID"])))
splitSets <- lapply(updated, function(f) 
  f[which(f$gencodeID %in% geneuniverse), ])

elementNROWS(splitSets)

di = lapply(splitSets, function(x) 
  introns[which(introns$genes %in% x$ensemblID),])
elementNROWS(di)
names(di) = c("Intellectual\nDisability","Neuro-\ndevel.","Neuro-\ndegen.",
              "SCZ\n(GWAS)","BPAD\n(GWAS)","ASD\n(SFARI)","ASD\n(CNV)",
              "SCZ\n(CNV)","SCZ\n(SNV)")
di = do.call(rbind, Map(cbind, Disease = as.list(names(di)), di))
di$Disease <- factor(di$Disease, levels = c("ASD\n(CNV)","SCZ\n(CNV)",
                                            "ASD\n(SFARI)","BPAD\n(GWAS)",
                                            "SCZ\n(SNV)","Neuro-\ndegen.",
                                            "SCZ\n(GWAS)", "Neuro-\ndevel.",
                                            "Intellectual\nDisability"))
di$SampleGroup = factor(di$SampleGroup, 
                        levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm",
                                   "Adult:Nucleus","Prenatal:Nucleus"))
head(di)
di$NuclearEnriched = NA
di[which(di$Disease %in% c("ASD\n(CNV)","SCZ\n(CNV)")),
   "NuclearEnriched"] = "Nuclear in Both"
di[which(di$Disease %in% c("ASD\n(SFARI)","BPAD\n(GWAS)","SCZ\n(SNV)",
                           "Neuro-\ndegen.")), "NuclearEnriched"] = "Nuclear in Adult Only"
di[which(di$Disease %in% c("SCZ\n(GWAS)", "Neuro-\ndevel.",
                           "Intellectual\nDisability")), 
   "NuclearEnriched"] = "Not Enriched"
di$NuclearEnriched = factor(di$NuclearEnriched, 
                            levels = c("Nuclear in Both", "Nuclear in Adult Only",
                                       "Not Enriched"))

res = list("Nuclear in Both" = t.test(di[which(di$NuclearEnriched=="Nuclear in Both"),
                                         "IRratio"], 
                                      di[which(di$NuclearEnriched=="Not Enriched"),"IRratio"]),
           "Nuclear in Adult Only" = t.test(di[which(di$NuclearEnriched=="Nuclear in Adult Only"),
                                               "IRratio"], 
                                            di[which(di$NuclearEnriched=="Not Enriched"),
                                               "IRratio"]),
           "Nuclear in Both or Adult Only" = t.test(di[which(di$NuclearEnriched!="Not Enriched"),
                                                       "IRratio"], 
                                                    di[which(di$NuclearEnriched=="Not Enriched"),
                                                       "IRratio"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) 
                           data.frame(Tstat = x$statistic, mean.NucGroup = x$estimate[1], 
                                      mean.Other = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")

res
#                           Group     Tstat mean.NucGroup mean.Other          pval           FDR
#t                Nuclear in Both  10.12710    0.04330846 0.03402363  4.380129e-24  4.380129e-24
#t1         Nuclear in Adult Only -22.73725    0.02326721 0.03402363 2.869546e-114 8.608637e-114
#t2 Nuclear in Both or Adult Only -18.04867    0.02550163 0.03402363  9.492504e-73  1.423876e-72

pdf(paste0("./Dropbox/sorted_figures/github_controlled/disease/figures/",
           "diseaseGene_IRratios.pdf"), width = 11, height = 10.5)
ggplot(di, aes(x = Disease, y = IRratio, fill = NuclearEnriched)) + 
  geom_jitter(size=0.5) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  facet_grid(SampleGroup ~ .) + 
  labs(fill="") + ylab("IR Ratio") + xlab("") +
  ggtitle("IR ratios of disease-associated gene introns") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "bottom")
dev.off()

di$NuclearEnriched <- gsub("Nuclear in Both", "Nuclear\nin Both", di$NuclearEnriched)
di$NuclearEnriched <- gsub("Nuclear in Adult Only", "Nuclear\nin Adult Only", di$NuclearEnriched)
di$NuclearEnriched <- factor(di$NuclearEnriched, 
                             levels = c("Nuclear\nin Both", "Nuclear\nin Adult Only", 
                                        "Not Enriched"))

pdf(paste0("./Dropbox/sorted_figures/github_controlled/disease/figures/",
           "diseaseGene_IRratios_collapsed.pdf"),width=6,height=3.5)
ggplot(di, aes(x = NuclearEnriched, y = IRratio, fill = NuclearEnriched)) + 
  geom_jitter(size=0.5) +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  labs(fill="") + ylab("IR Ratio") + xlab("") +
  ggtitle("Disease-Associated Gene Introns") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "none")
dev.off()


## intron number 

di$NuclearEnriched <- gsub("Nuclear\nin Both", "Nuclear in Both", di$NuclearEnriched)
di$NuclearEnriched <- gsub("Nuclear\nin Adult Only", "Nuclear in Adult Only", 
                           di$NuclearEnriched)
di$NuclearEnriched = factor(di$NuclearEnriched, 
                            levels = c("Nuclear in Both", 
                                       "Nuclear in Adult Only", "Not Enriched"))

library(data.table)
num = data.table(di)[,length(unique(intronID)), by = c("NuclearEnriched","genes")]
num=data.frame(num)

res = list("Nuclear in Both" = t.test(num[which(num$NuclearEnriched=="Nuclear in Both"),"V1"], 
                                      num[which(num$NuclearEnriched=="Not Enriched"),"V1"]),
           "Nuclear in Adult Only" = t.test(num[which(num$NuclearEnriched=="Nuclear in Adult Only"),"V1"], 
                                            num[which(num$NuclearEnriched=="Not Enriched"),"V1"]),
           "Nuclear in Both or Adult Only" = t.test(num[which(num$NuclearEnriched!="Not Enriched"),"V1"], 
                                                    num[which(num$NuclearEnriched=="Not Enriched"),"V1"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, 
                                                            mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], 
                                                            pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group      Tstat mean.NucGroup mean.Other         pval          FDR
#t                Nuclear in Both -0.9323997      13.49444   14.46242 3.517945e-01 3.517945e-01
#t1         Nuclear in Adult Only  6.5216689      19.45469   14.46242 9.377900e-11 2.813370e-10
#t2 Nuclear in Both or Adult Only  5.7775149      18.69110   14.46242 9.240154e-09 1.386023e-08


## intron length

di$width = width(makeGRangesFromDataFrame(di, strand="Direction"))
head(di)

res = list("Nuclear in Both" = t.test(di[which(di$NuclearEnriched=="Nuclear in Both"),"width"],
                                      di[which(di$NuclearEnriched=="Not Enriched"),"width"]),
           "Nuclear in Adult Only" = t.test(di[which(di$NuclearEnriched=="Nuclear in Adult Only"),"width"], 
                                            di[which(di$NuclearEnriched=="Not Enriched"),"width"]),
           "Nuclear in Both or Adult Only" = t.test(di[which(di$NuclearEnriched!="Not Enriched"),"width"], 
                                                    di[which(di$NuclearEnriched=="Not Enriched"),"width"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, 
                                                            mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], 
                                                            pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group      Tstat mean.NucGroup mean.Other          pval           FDR
#t                Nuclear in Both -29.976144      8413.701    13469.1 1.257264e-196 3.771791e-196
#t1         Nuclear in Adult Only  -3.375955     13001.721    13469.1  7.357380e-04  7.357380e-04
#t2 Nuclear in Both or Adult Only  -7.250872     12490.198    13469.1  4.156991e-13  6.235486e-13


## maximum IR for each gene

irbyGene = data.frame(data.table(di)[, list(IRratio=max(IRratio)), 
                                     by=c("genes","SampleID","Disease",
                                          "NuclearEnriched","SampleGroup")]) 

res = list("Nuclear in Both" = 
             t.test(irbyGene[which(irbyGene$NuclearEnriched=="Nuclear in Both"),"IRratio"], 
                    irbyGene[which(irbyGene$NuclearEnriched!="Nuclear in Both"),"IRratio"]),
           "Nuclear in Adult Only" = 
             t.test(irbyGene[which(irbyGene$NuclearEnriched=="Nuclear in Adult Only"),"IRratio"], 
                    irbyGene[which(irbyGene$NuclearEnriched!="Nuclear in Adult Only"),"IRratio"]),
           "Nuclear in Both or Adult Only" = 
             t.test(irbyGene[which(irbyGene$NuclearEnriched!="Not Enriched"),"IRratio"], 
                    irbyGene[which(irbyGene$NuclearEnriched=="Not Enriched"),"IRratio"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, 
                                                            mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], 
                                                            pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group      Tstat mean.NucGroup mean.Other         pval          FDR
#t                Nuclear in Both  11.774232     0.3162994  0.2290954 1.977336e-31 5.932007e-31
#t1         Nuclear in Adult Only -11.443056     0.2181164  0.2698305 3.117899e-30 4.676849e-30
#t2 Nuclear in Both or Adult Only  -3.595372     0.2340065  0.2515110 3.250631e-04 3.250631e-04

pdf(paste0("./Dropbox/sorted_figures/github_controlled/disease/figures/",
           "diseaseGene_maxIRratios.pdf"), width = 11, height = 10.5)
ggplot(irbyGene, aes(x = Disease, y = IRratio, fill = NuclearEnriched)) + 
  geom_jitter(size=0.5) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  facet_grid(SampleGroup ~ .) + 
  labs(fill="") + ylab("IR Ratio") + xlab("") +
  ggtitle("IR ratios of disease-associated gene introns") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "bottom")
dev.off()

