library("readxl")
library('DESeq2')
library("ggplot2")

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")

percrRNA = df$rRNAreads / df$totalMapped * 100
percMapped = df$totalMapped / df$TotalReads * 100
percMito = df$mitoMapped / df$totalMapped * 100
Unmapped = 100 - percMapped
MeanRpkm = colMeans(geneRpkm)

d = data.frame("SampleID" =df$SampleID, percMapped, Unmapped, percrRNA, percMito, 
               MeanRpkm, "Label" = df$Label, TotalReads = df$TotalReads, totalMapped = df$totalMapped)
VSD = varianceStabilizingTransformation(geneCounts)
VSD = data.frame(VSD)
pVSD = varianceStabilizingTransformation(polya.counts)
pVSD = data.frame(pVSD)
rVSD = varianceStabilizingTransformation(ribozero.counts)
rVSD = data.frame(rVSD)
geneRpkm = data.frame(geneRpkm)

##Plotting by Sample
library(ggplot2)

# Variance-stabilized Count Density Plots
ggplot(pVSD) + geom_density(aes(x=Br1113C1_polyA)) + geom_density(aes(x=Br1113N1_polyA)) + 
  geom_density(aes(x=Br2046C_polyA)) + geom_density(aes(x=Br2046N_polyA)) + 
  geom_density(aes(x=Br2074C_polyA)) + geom_density(aes(x=Br2074N_polyA)) +
  geom_density(aes(x=Br5339C1_polyA)) + geom_density(aes(x=Br5339N1_polyA)) +
  geom_density(aes(x=Br5340C1_polyA)) + geom_density(aes(x=Br5340N1_polyA)) +
  geom_density(aes(x=Br5341C1_polyA)) + geom_density(aes(x=Br5341N1_polyA)) + 
  ggtitle("Variance-stabilized counts: PolyA") + theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) + xlab("") + ylim(0,.675)
  
ggplot(rVSD) + geom_density(aes(x=Br1113C1_RiboZero)) + 
  geom_density(aes(x=Br2046C_RiboZero)) + geom_density(aes(x=Br2046N_RiboZero)) + 
  geom_density(aes(x=Br2074C_RiboZero)) + geom_density(aes(x=Br2074N_RiboZero)) +
  geom_density(aes(x=Br5339C1_RiboZero)) + geom_density(aes(x=Br5339N1_RiboZero)) +
  geom_density(aes(x=Br5340C1_RiboZero)) + geom_density(aes(x=Br5340N1_RiboZero)) +
  geom_density(aes(x=Br5341C1_RiboZero)) + geom_density(aes(x=Br5341N1_RiboZero)) + 
  ggtitle("Variance-stabilized counts: RiboZero") + theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) + xlab("") + ylim(0,.675)
  
ggplot(geneRpkm) + geom_density(aes(x=Br1113C1_polyA)) + geom_density(aes(x=Br1113N1_polyA)) +
  geom_density(aes(x=Br2046C_polyA)) + geom_density(aes(x=Br2046N_polyA)) + 
  geom_density(aes(x=Br2074C_polyA)) + geom_density(aes(x=Br2074N_polyA)) +
  geom_density(aes(x=Br5339C1_polyA)) + geom_density(aes(x=Br5339N1_polyA)) +
  geom_density(aes(x=Br5340C1_polyA)) + geom_density(aes(x=Br5340N1_polyA)) +
  geom_density(aes(x=Br5341C1_polyA)) + geom_density(aes(x=Br5341N1_polyA)) + xlim(0,5)
  
  ggplot(geneRpkm) + geom_density(aes(x=Br1113C1_RiboZero)) + 
  geom_density(aes(x=Br2046C_RiboZero)) + geom_density(aes(x=Br2046N_RiboZero)) + 
  geom_density(aes(x=Br2074C_RiboZero)) + geom_density(aes(x=Br2074N_RiboZero)) +
  geom_density(aes(x=Br5339C1_RiboZero)) + geom_density(aes(x=Br5339N1_RiboZero)) +
  geom_density(aes(x=Br5340C1_RiboZero)) + geom_density(aes(x=Br5340N1_RiboZero)) +
  geom_density(aes(x=Br5341C1_RiboZero)) + geom_density(aes(x=Br5341N1_RiboZero))

# Total reads
ggplot(d, aes(x= SampleID, y=TotalReads)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Total Reads") + 
  xlab("") +
  ggtitle("QC Metrics: Total Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Total Mapped Reads
ggplot(d, aes(x= SampleID, y=totalMapped)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Total Mapped Reads") + 
  xlab("") +
  ggtitle("QC Metrics: Total Mapped Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# % Unmapped reads
ggplot(d, aes(x= SampleID, y=Unmapped)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("QC Metrics: Percent Unmapped Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Label)

# % Mitochondrial reads
ggplot(d, aes(x= SampleID, y=percMito)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Reads Mapping\nTo Mitochondrial Genome") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# % rRNA reads
ggplot(d, aes(x= SampleID, y=percrRNA)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Reads Mapping\nTo Ribosomal RNA") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

##Plotting by group

# Total Reads
ggplot(d, aes(x=Label, y=TotalReads)) + geom_boxplot() + 
  geom_jitter() +
  ylab("Total Reads") + 
  xlab("") +
  ggtitle("QC Metrics: Total Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Total Mapped Reads
ggplot(d, aes(x=Label, y=totalMapped)) + geom_boxplot() + 
  geom_jitter() +
  ylab("Total Mapped Reads") + 
  xlab("") +
  ggtitle("QC Metrics: Total Mapped Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# % Unmapped reads
ggplot(d, aes(x=Label, y=Unmapped)) + geom_boxplot() + 
  geom_jitter() +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Unmapped Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# % Mitochondrial reads
ggplot(d, aes(x=Label, y=percMito)) + geom_boxplot() + 
  geom_jitter() +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Reads Mapping\nTo Mitochondrial Genome") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# % rRNA reads
ggplot(d, aes(x=Label, y=percrRNA)) + geom_boxplot() + 
  geom_jitter() +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Reads Mapping\nTo Ribosomal RNA") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Mean RPKM of Expressed Genes by Group
ggplot(d, aes(x=Label, y=MeanRpkm)) + geom_boxplot() + 
  geom_jitter() + 
  ylab("Mean RPKM") + 
  xlab("") +
  ggtitle("Mean RPKM of Expressed Genes By Group") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))