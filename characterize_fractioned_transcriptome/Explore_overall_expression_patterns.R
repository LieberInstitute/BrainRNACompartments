library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(data.table)


load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")


### Check how comparable gene expression is between the fractions

degs = list("By Fraction In Adult PolyA" = Apres, "By Fraction In Prenatal PolyA" = Fpres.down, "By Age In Cytosol PolyA" = Cpres.down, "By Age In Nucleus PolyA" = Npres)
elementNROWS(degs)
#AP    FP    CP    NP 
#44124 43610 43610 44124 
write.csv(data.frame(rbind(elementNROWS(degs), elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05),])),
                           elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05 & x$log2FoldChange<0),])),
                           elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05 & x$log2FoldChange>0),])),
                           paste0(round(elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05),]))/elementNROWS(degs)*100,2), "%"),
                           paste0(round(elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05 & x$log2FoldChange<0),]))/
                           elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05),]))*100,2), "%"),
                           paste0(round(elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05 & x$log2FoldChange>0),]))/
                           elementNROWS(lapply(degs, function(x) x[which(x$padj <=0.05),]))*100,2), "%")), 
                    row.names = c("Total genes measured", "Total DEGs", "DEGs with (-)LFC", "DEGs with (+)LFC","Percent of Genes that are DEGs",
                                  "Percent of DEGs with (-)LFC","Percent of DEGs with (+)LFC")), quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/TotalGeneNumbers_thatAreDEG.csv")
#                               By.Fraction.In.Adult.PolyA By.Fraction.In.Prenatal.PolyA By.Age.In.Cytosol.PolyA By.Age.In.Nucleus.PolyA
#Total genes measured                                44124                         43610                   43610                   44124
#Total DEGs                                           4643                           343                   10870                   10206
#DEGs with (-)LFC                                     2579                            71                    5681                    5324
#DEGs with (+)LFC                                     2064                           272                    5189                    4882
#Percent of Genes that are DEGs                     10.52%                         0.79%                  24.93%                  23.13%
#Percent of DEGs with (-)LFC                        55.55%                         20.7%                  52.26%                  52.17%
#Percent of DEGs with (+)LFC                        44.45%                         79.3%                  47.74%                  47.83%



## Assess gene expression by fraction

counts = melt(geneCounts.down/(colSums(geneCounts.down)/1000000))
counts = data.table(counts[grep("polyA",counts$Var2),])
counts[grep("C", counts$Var2),"Fraction"] = "Cytosol"
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("53", counts$Var2),"Age"] = "Prenatal"
counts[-grep("53", counts$Var2),"Age"] = "Adult"
counts$group = paste(counts$Age, counts$Fraction, sep=":")
colnames(counts) = c("Gene","SampleID","Count","Fraction","Age","Group")

means = counts[,list(mean(Count)),by=c("Gene","Fraction","Age","Group")]
df = t(data.frame(means[V1==0,length(unique(Gene)), by="Group"][,list(round(V1/43610,2)),],
                  means[,mean(V1), by="Group"][,list(round(V1,2)),],
                  means[,median(V1), by="Group"][,list(round(V1,2)),],
                  means[,sd(V1), by="Group"][,list(round(V1,2)),],
                  means[,max(V1), by="Group"][,list(round(V1,2)),],
                  means[,min(V1), by="Group"][,list(round(V1,2)),],
                  means[V1<183380,mean(V1), by="Group"][,list(round(V1,2)),],
                  means[V1<183380,median(V1), by="Group"][,list(round(V1,2)),],
                  means[V1<183380,sd(V1), by="Group"][,list(round(V1,2)),],
                  means[V1<183380,max(V1), by="Group"][,list(round(V1,2)),],
                  means[V1<183380,min(V1), by="Group"][,list(round(V1,2)),]))
colnames(df) = as.character(as.data.frame(means[,mean(V1), by="Group"][,list(as.character(Group)),])[,"V1"])
df = data.frame(Value = c("Proportion total genes not expressed in each group","Mean Gene RPM","Median Gene RPM","SD Gene RPM","Max Gene RPM","Min Gene RPM",
                          "Mean RPM (excluding AN outlier)","Median RPM (excluding AN outlier)","SD RPM (excluding AN outlier)",
                          "Max RPM (excluding AN outlier)","Min RPM (excluding AN outlier)"),df)
write.csv(df, quote = F,row.names = F,file = "./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/Gene_RPM_stats_byGroup.csv")
           


ggplot(means, aes(x = Fraction, y = V1)) + geom_boxplot() +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("RPM") + ylim(0,1) + 
  xlab("") +
  ggtitle("Gene Read Counts Per Million Mapped") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
head(means$V1)











