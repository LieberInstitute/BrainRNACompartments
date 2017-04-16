library(pheatmap)
library(RColorBrewer)
library(DESeq2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

rlog  = rlog(geneCounts)
rlog.down = rlog(geneCounts.down)
save(rlog, rlog.down, 
     file="./Dropbox/sorted_figures/new/github_controlled/characterize_transcriptome/data/rlog_transformed_counts.rda")

# Euclidean distance between samples heatmap
sampleDists <- dist(t(rlog))
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
match(colnames(geneCounts),pd$SampleID)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = pd$WorkingID
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples")

sampleDists <- dist(t(rlog.down))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = pd$WorkingID
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples (Downsampled)")
