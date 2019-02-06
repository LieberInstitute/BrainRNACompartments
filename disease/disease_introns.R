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
IRres = Map(cbind, IRres, genes = lapply(lapply(IRres, function(x) unlist(strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE), recursive = F)), 
                                         function(x) x[grep("ENSG", x)]), 
                 intronID = lapply(IRres, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
elementNROWS(IRres)
head(IRres[[1]])
introns = do.call(rbind, Map(cbind, lapply(IRres, function(x) data.frame(x[,which(colnames(x) %in% 
                                         c("Chr","Start","End","Direction","IRratio","Warnings","genes","intronID"))])), SampleID = as.list(names(IRres))))
introns$Fraction = ifelse(introns$SampleID %in% c("Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1"), "Nucleus","Cytoplasm")
introns$Age = ifelse(introns$SampleID %in% c("Br5339C1","Br5339N1","Br5340C1","Br5340N1","Br5341C1","Br5341N1"), "Prenatal","Adult")
introns$SampleGroup = paste(introns$Age, introns$Fraction, sep=":")
head(introns)


## Get disease gene introns

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

AEJmap = geneMap[which(geneMap$Symbol %in% as.character(aej_sets$Gene.Symbol)),]
AEJmap = cbind(AEJmap, Gene.Set = aej_sets[match(AEJmap$Symbol, aej_sets$Gene.Symbol),"Gene.Set"])
AEJmap = split(AEJmap, AEJmap$Gene.Set)
elementNROWS(AEJmap)

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
pgc2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]

AEJmap = c(AEJmap[names(AEJmap)!="SCZ PGC GWAS"], list(PGC2 = pgc2))

di = lapply(AEJmap, function(x) introns[which(introns$genes %in% x$ensemblID),])
names(di) = c("ASD\n(CNV)","ASD\n(Database)","BPAD\n(GWAS)","Intellectual\nDisability","Neuro-\ndevel.","Neuro-\ndegen.",
              "SCZ\n(CNV)","SCZ\n(Meta\nanalysis)","SCZ\n(SNV)","SCZ\n(PGC2)")
di = do.call(rbind, Map(cbind, Disease = as.list(names(di)), di))
di$Disease = factor(di$Disease, levels = 
                        c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)","BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)","SCZ\n(Meta\nanalysis)",
                          "Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))
di$SampleGroup = factor(di$SampleGroup, levels = c("Adult:Cytoplasm","Prenatal:Cytoplasm","Adult:Nucleus","Prenatal:Nucleus"))
head(di)
di$NuclearEnriched = NA
di[which(di$Disease %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)")), "NuclearEnriched"] = "Nuclear in Both"
di[which(di$Disease %in% c("BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")), "NuclearEnriched"] = "Nuclear in Adult Only"
di[-which(di$Disease %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)","BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")), "NuclearEnriched"] = "Not Enriched"
di$NuclearEnriched = factor(di$NuclearEnriched, levels = c("Nuclear in Both", "Nuclear in Adult Only", "Not Enriched"))

res = list("Nuclear in Both" = t.test(di[which(di$NuclearEnriched=="Nuclear in Both"),"IRratio"], di[which(di$NuclearEnriched!="Nuclear in Both"),"IRratio"]),
           "Nuclear in Adult Only" = t.test(di[which(di$NuclearEnriched=="Nuclear in Adult Only"),"IRratio"], di[which(di$NuclearEnriched!="Nuclear in Adult Only"),"IRratio"]),
           "Nuclear in Both or Adult Only" = t.test(di[which(di$NuclearEnriched!="Not Enriched"),"IRratio"], di[which(di$NuclearEnriched=="Not Enriched"),"IRratio"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group     Tstat mean.NucGroup mean.Other         pval          FDR
#                 Nuclear in Both  4.857449    0.02956443 0.02680487 1.190265e-06 1.190265e-06
#           Nuclear in Adult Only  7.510002    0.02971083 0.02580183 5.932854e-14 8.899282e-14
#   Nuclear in Both or Adult Only 19.879623    0.02965155 0.01837064 1.020760e-87 3.062279e-87

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/diseaseGene_IRratios.pdf",width=11,height=10.5)
ggplot(di, aes(x = Disease, y = IRratio, fill = NuclearEnriched)) + geom_jitter(size=0.5) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  facet_grid(SampleGroup ~ .) + labs(fill="") + ylab("IR Ratio") + xlab("") +
  ggtitle("IR ratios of disease-associated gene introns") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "bottom")
dev.off()

## maximum IR for each gene

irbyGene = data.frame(data.table(di)[, list(IRratio=max(IRratio)), by=c("genes","SampleID","Disease","NuclearEnriched","SampleGroup")]) 

res = list("Nuclear in Both" = t.test(irbyGene[which(irbyGene$NuclearEnriched=="Nuclear in Both"),"IRratio"], 
                                      irbyGene[which(irbyGene$NuclearEnriched!="Nuclear in Both"),"IRratio"]),
           "Nuclear in Adult Only" = t.test(irbyGene[which(irbyGene$NuclearEnriched=="Nuclear in Adult Only"),"IRratio"], 
                                            irbyGene[which(irbyGene$NuclearEnriched!="Nuclear in Adult Only"),"IRratio"]),
           "Nuclear in Both or Adult Only" = t.test(irbyGene[which(irbyGene$NuclearEnriched!="Not Enriched"),"IRratio"], 
                                                    irbyGene[which(irbyGene$NuclearEnriched=="Not Enriched"),"IRratio"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                          Group     Tstat mean.NucGroup mean.Other         pval          FDR
#                Nuclear in Both 0.8810247     0.2414038  0.2357259 3.783281e-01 3.783281e-01
#          Nuclear in Adult Only 4.1566971     0.2500462  0.2248785 3.248700e-05 4.873050e-05
#  Nuclear in Both or Adult Only 7.4050341     0.2466148  0.1905496 1.652441e-13 4.957323e-13

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/diseaseGene_maxIRratios.pdf",width=11,height=10.5)
ggplot(irbyGene, aes(x = Disease, y = IRratio, fill = NuclearEnriched)) + geom_jitter(size=0.5) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  facet_grid(SampleGroup ~ .) + labs(fill="") + ylab("IR Ratio") + xlab("") +
  ggtitle("IR ratios of disease-associated gene introns") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "bottom")
dev.off()

