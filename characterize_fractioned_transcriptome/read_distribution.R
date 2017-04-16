library(scales)
library(ggplot2)

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
load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/read_distribution_data.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

head(dist)
names = unique(totals$SampleID)
dist$label = pd[match(dist$SampleID,pd$SampleID),"Label"]
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
dist$Group = factor(dist$Group, levels=c("TES (10kb downstream)","TSS (10kb upstream)",
                                         "Introns","3'UTR","5'UTR","CDS Exons"))

# Plot Percent Tags across 6 features
ggplot(dist[which(dist$Group!="TSS_up_1kb" & dist$Group!="TSS_up_5kb" & dist$Group!="TES_down_1kb" & dist$Group!="TES_down_5kb" &
                    dist$SampleID!="Br5340C1_downsamp" & dist$SampleID!="Br5339C1_downsamp"),], 
                                     aes(x=label, y=Percent, fill=Group), color=Group) + 
    geom_bar(position = "fill",stat = "identity", width=0.75) + 
    scale_y_continuous(labels = percent_format()) +
    ylab("Percent") + 
    xlab("") + ggtitle("Percent of Reads Mapping to Six Genomic Features") +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))

# Plot Percent Tags across 6 features normalized to bases covered by each feature
ggplot(dist[which(dist$Group!="TSS_up_1kb" & dist$Group!="TSS_up_5kb" & dist$Group!="TES_down_1kb" & dist$Group!="TES_down_5kb" &
                    dist$SampleID!="Br5340C1_downsamp" & dist$SampleID!="Br5339C1_downsamp"),], 
       aes(x=label, y=Tags.Kb, fill=Group), color=Group) + 
  geom_bar(position = "fill",stat = "identity", width=0.75) + 
  scale_y_continuous(labels = percent_format()) +
  ylab("Percent") + 
  xlab("") + ggtitle("Percent of Tags Per Kb of Genomic Feature ") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Plot Percent Tags across 6 features (downsampled)
ggplot(dist[which(dist$Group!="TSS_up_1kb" & dist$Group!="TSS_up_5kb" & dist$Group!="TES_down_1kb" & dist$Group!="TES_down_5kb" &
                    dist$SampleID!="Br5340C1_polyA" & dist$SampleID!="Br5339C1_polyA"),], 
       aes(x=label, y=Percent, fill=Group), color=Group) + 
  geom_bar(position = "fill",stat = "identity", width=0.75) + 
  scale_y_continuous(labels = percent_format()) +
  ylab("Percent") + 
  xlab("") + ggtitle("Percent of Reads Mapping to Six Genomic Features") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Plot Percent Tags across 6 features normalized to bases covered by each feature
ggplot(dist[which(dist$Group!="TSS_up_1kb" & dist$Group!="TSS_up_5kb" & dist$Group!="TES_down_1kb" & dist$Group!="TES_down_5kb" &
                    dist$SampleID!="Br5340C1_polyA" & dist$SampleID!="Br5339C1_polyA"),], 
       aes(x=label, y=Tags.Kb, fill=Group), color=Group) + 
  geom_bar(position = "fill",stat = "identity", width=0.75) + 
  scale_y_continuous(labels = percent_format()) +
  ylab("Percent") + 
  xlab("") + ggtitle("Percent of Tags Per Kb of Genomic Feature ") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

