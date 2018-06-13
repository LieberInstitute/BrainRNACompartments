library(scales)
library(ggplot2)
library(data.table)

# Create R objects
names = scan("/media/DATA/Amanda/fullnames.txt", what = "character")
dist = total = list()
for (i in 1:length(names)){
  dist[[i]] = read.table(paste0("/media/DATA/Amanda/read_distribution/", names[i],".dist.txt"),
                         header = TRUE, 
                         skip = 4, comment.char = "=")
}
for (i in 1:length(names)){
  total[[i]] = scan(paste0("/media/DATA/Amanda/read_distribution/", names[i],".dist.txt"), nlines = 3, what = "character")
}
names(dist) = names(total) = names
dist[["Br5340C1_downsamp"]] = read.table("/media/DATA/Amanda/read_distribution/Br5340C1_downsamp.dist.txt",header = TRUE,skip=4,comment.char = "=")
dist[["Br5339C1_downsamp"]] = read.table("/media/DATA/Amanda/read_distribution/Br5339C1_downsamp.dist.txt",header = TRUE,skip=4,comment.char = "=")
total[["Br5340C1_downsamp"]] = scan("/media/DATA/Amanda/read_distribution/Br5340C1_downsamp.dist.txt", nlines = 3, what = "character")
total[["Br5339C1_downsamp"]] = scan("/media/DATA/Amanda/read_distribution/Br5339C1_downsamp.dist.txt", nlines = 3, what = "character")

values = lapply(total, function(x) as.numeric(x[c(3,6,10)]))
total = lapply(values, function(x) data.frame(total = c("Total.Reads","Total.Tags","Total.Assigned.Tags"), values = x))
for (i in 1:length(total)){total[[i]]["SampleID"] = names(total)[i]}
for (i in 1:length(dist)){dist[[i]]["SampleID"] = names(dist)[i]}
dist = do.call(rbind, dist)
totals = do.call(rbind, total)
save(dist, totals, file="/media/DATA/Amanda/read_distribution/read_distribution_data.rda")

# load and format objects
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/read_distribution_data.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

head(dist)
names = unique(totals$SampleID)
dist[grep("C", dist$SampleID), "Fraction"] = "Cytoplasm"
dist[grep("N", dist$SampleID), "Fraction"] = "Nucleus"
dist[grep("53", dist$SampleID), "Age"] = "Prenatal"
dist[-grep("53", dist$SampleID), "Age"] = "Adult"
dist[c(grep("poly", dist$SampleID),grep("down", dist$SampleID)), "Library"] = "polyA"
dist[grep("Ribo", dist$SampleID), "Library"] = "RiboZero"
dist$label = factor(paste(dist$Age, dist$Fraction, dist$Library, sep="\n"), 
                    levels = c("Prenatal\nCytoplasm\npolyA", "Prenatal\nNucleus\npolyA", "Adult\nCytoplasm\npolyA","Adult\nNucleus\npolyA",   
                               "Prenatal\nCytoplasm\nRiboZero","Prenatal\nNucleus\nRiboZero","Adult\nCytoplasm\nRiboZero","Adult\nNucleus\nRiboZero")) 
for (i in 1:nrow(dist)){
  dist[i,"Percent"] = dist[i,"Tag_count"] / totals[which(totals$total=="Total.Assigned.Tags" & totals$SampleID==dist[i,"SampleID"]),"values"] * 100   
}
for (i in 1:nrow(dist)){
  dist[i,"Percent.Kb"] = dist[i,"Tags.Kb"] / sum(dist[which(dist$SampleID==dist[i,"SampleID"]),"Tags.Kb"]) * 100   
}
dist$Group = gsub("3'UTR_Exons","3'UTR", dist$Group)
dist$Group = gsub("5'UTR_Exons","5'UTR", dist$Group)
dist$Group = gsub("CDS_Exons","CDS Exons", dist$Group)
dist$Group = gsub("TES_down_10kb","TES (10kb downstream)", dist$Group)
dist$Group = gsub("TSS_up_10kb","TSS (10kb upstream)", dist$Group)
dist = data.table(dist)
dist = dist[Group!="TSS_up_1kb" & Group!="TSS_up_5kb" & Group!="TES_down_1kb" & Group!="TES_down_5kb" &
              SampleID!="Br5340C1_polyA" & SampleID!="Br5339C1_polyA",,]
dist$Group = factor(dist$Group, levels=c("TES (10kb downstream)","TSS (10kb upstream)",
                                         "Introns","3'UTR","5'UTR","CDS Exons"))

perc = dist[, mean(Percent), by=c("Group", "label")]
Tags.Kb = dist[, mean(Tags.Kb), by=c("Group", "label")]

# Plot read distribution (using 2 downsampled counts)
pdf("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/figures/read_distribution_6_features.pdf", width = 12, height = 6)
ggplot(perc, aes(x=label, y=V1, fill=Group), color=Group) + 
    geom_bar(position = "fill",stat = "identity", width=0.75) + 
    scale_y_continuous(labels = percent_format()) +
    ylab("Percent") + 
    xlab("") + ggtitle("Percent of Reads Mapping to Six Genomic Features") +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(Tags.Kb, aes(x=label, y=V1, fill=Group), color=Group) + 
  geom_bar(position = "fill",stat = "identity", width=0.75) + 
  scale_y_continuous(labels = percent_format()) +
  ylab("Percent") + 
  xlab("") + ggtitle("Percent of Tags Per Kb of Genomic Feature ") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


## Calculate T stats

# by intron counts

ttests = list(introns.perc.both = t.test(x = dist[Fraction=="Nucleus" & Group =="Introns",list(Percent),], 
                                         y = dist[Fraction=="Cytoplasm" & Group =="Introns",list(Percent),]),
              introns.perc.polyA = t.test(x = dist[Library=="polyA" & Fraction=="Nucleus" & Group =="Introns",list(Percent),], 
                                          y = dist[Library=="polyA" & Fraction=="Cytoplasm" & Group =="Introns",list(Percent),]),
              introns.perc.ribo = t.test(x = dist[Library=="RiboZero" & Fraction=="Nucleus" & Group =="Introns",list(Percent),], 
                                         y = dist[Library=="RiboZero" & Fraction=="Cytoplasm" & Group =="Introns",list(Percent),]),
              introns.tags.kb.both = t.test(x = dist[Fraction=="Nucleus" & Group =="Introns",list(Tags.Kb),], 
                                         y = dist[Fraction=="Cytoplasm" & Group =="Introns",list(Tags.Kb),]),
              introns.tags.kb.polyA = t.test(x = dist[Library=="polyA" & Fraction=="Nucleus" & Group =="Introns",list(Tags.Kb),], 
                                          y = dist[Library=="polyA" & Fraction=="Cytoplasm" & Group =="Introns",list(Tags.Kb),]),
              introns.tags.kb.ribo = t.test(x = dist[Library=="RiboZero" & Fraction=="Nucleus" & Group =="Introns",list(Tags.Kb),], 
                                         y = dist[Library=="RiboZero" & Fraction=="Cytoplasm" & Group =="Introns",list(Tags.Kb),]))
df = data.frame(test = names(ttests), p.value = unlist(lapply(ttests, function(x) x$p.value)),
                t.stat = unlist(lapply(ttests, function(x) x$statistic)),
                mean.x = unlist(lapply(ttests, function(x) x$estimate[1])),
                mean.y = unlist(lapply(ttests, function(x) x$estimate[2])), row.names = NULL)
df$FDR = p.adjust(df$p.value, method = "fdr")
write.csv(df, file = "./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/read_distribution_intron_tstat.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/read_distribution_intron_tstat.csv")
df
