library(GenomicRanges)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Calculate the correlation coefficients between gene log2(count+1) between ribozero and polyA 
head(pd)
names = unique(gsub("_.*","",colnames(geneCounts)))
names = names[c(1,3:12)]
geneCounts = as.data.frame(geneCounts[,c(1:2,4:23)])

pdf("./Dropbox/sorted_figures/new/github_controlled/QC_section/figures/geneCounts_by_library.pdf", width=8, height=8)
cor = cov = counts = list()
for (i in 1:length(names)){
  cor[[i]] = cor(log2(geneCounts[,paste0(names[i],"_polyA")]+1), log2(geneCounts[,paste0(names[i],"_RiboZero")]+1))
  cov[[i]] = cov(log2(geneCounts[,paste0(names[i],"_polyA")]+1), log2(geneCounts[,paste0(names[i],"_RiboZero")])+1)
  counts[[i]] = log2(geneCounts[,grep(names[i], colnames(geneCounts))]+1)
  g = ggplot(counts[[i]], aes(x=counts[[i]][,1], y=counts[[i]][,2])) + 
    geom_point(alpha=1/4) +
    geom_smooth(method=lm) +
    ylab("RiboZero") + 
    xlab("PolyA") +
    ggtitle(paste0("Gene-Level log2(count+1): ",names[i], "\ncor = ",round(cor[[i]],digits =3),"\ncov = ",round(cov[[i]],digits =3))) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="")
  print(g)
}
dev.off()

# Calculate the correlation coefficients between exon log2(count+1) between ribozero and polyA
head(exonCounts)
exonCounts = as.data.frame(exonCounts)
pdf("./Dropbox/sorted_figures/new/github_controlled/QC_section/figures/exonCounts_by_library.pdf", width=8, height=8)
cor = cov = counts = list()
for (i in 1:length(names)){
  cor[[i]] = cor(log2(exonCounts[,paste0(names[i],"_polyA")]+1), log2(exonCounts[,paste0(names[i],"_RiboZero")]+1))
  cov[[i]] = cov(log2(exonCounts[,paste0(names[i],"_polyA")]+1), log2(exonCounts[,paste0(names[i],"_RiboZero")])+1)
  counts[[i]] = log2(exonCounts[,grep(names[i], colnames(exonCounts))]+1)
  g = ggplot(counts[[i]], aes(x=counts[[i]][,1], y=counts[[i]][,2])) + 
    geom_point(alpha=1/4) +
    geom_smooth(method=lm) +
    ylab("RiboZero") + 
    xlab("PolyA") +
    ggtitle(paste0("Exon-Level log2(count+1): ",names[i], "\ncor = ",round(cor[[i]],digits =3),"\ncov = ",round(cov[[i]],digits =3))) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="")
  print(g)
}
dev.off()