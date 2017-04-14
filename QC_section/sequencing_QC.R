library("readxl")
library('DESeq2')
library("ggplot2")

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")

MeanRpkm = c(colMeans(geneRpkm), colMeans(geneRpkm.down[,c(12,16)]))
names(MeanRpkm) = c(names(MeanRpkm[1:23]), "Br5339C1_downsamp", "Br5340C1_downsamp")
MeanRpkm = MeanRpkm[order(names(MeanRpkm))]
metrics = metrics[order(metrics$SAMPLE_ID),]
metrics = metrics[which(metrics$SAMPLE_ID!="Br1113N1_RiboZero"),]
pd = rbind(pd[1:11,], pd[12,], pd[12:15,], pd[16,], pd[16:23,])
pheno.table = cbind(pd, metrics)
pheno.table$SampleID = pheno.table$SAMPLE_ID
write.csv(pheno.table[,c(1:18,20:65)], 
          file="/Users/amanda/Dropbox/sorted_figures/new/github_controlled/QC_section/data/pheno.seq.data.csv",
          row.names = F, quote = F)

d = data.frame(SampleID =metrics$SAMPLE_ID, 
               percMapped = metrics$overallMapRate * 100,
               Unmapped = 100 - metrics$overallMapRate * 100,
               percrRNA = metrics$rRNA_rate * 100,
               percMito = metrics$mitoRate * 100, 
               TotalReads = metrics$numReads,
               totalMapped = metrics$numMapped,
               MeanRpkm, 
               Label = pd$Label)
geneRpkm = data.frame(geneRpkm)
geneRpkm.down = data.frame(geneRpkm.down)

##Plotting by Sample

# Variance-stabilized Count Density Plots
ggplot(geneRpkm) + geom_density(aes(x=Br1113C1_polyA)) + geom_density(aes(x=Br1113N1_polyA)) +
  geom_density(aes(x=Br2046C_polyA)) + geom_density(aes(x=Br2046N_polyA)) + 
  geom_density(aes(x=Br2074C_polyA)) + geom_density(aes(x=Br2074N_polyA)) +
  geom_density(aes(x=Br5339C1_polyA)) + geom_density(aes(x=Br5339N1_polyA)) +
  geom_density(aes(x=Br5340C1_polyA)) + geom_density(aes(x=Br5340N1_polyA)) +
  geom_density(aes(x=Br5341C1_polyA)) + geom_density(aes(x=Br5341N1_polyA)) + xlim(0,5) +
  ggtitle("Density Plot (RPKM): PolyA") + xlab("RPKM")

ggplot(geneRpkm.down) + geom_density(aes(x=Br1113C1_polyA)) + geom_density(aes(x=Br1113N1_polyA)) +
  geom_density(aes(x=Br2046C_polyA)) + geom_density(aes(x=Br2046N_polyA)) + 
  geom_density(aes(x=Br2074C_polyA)) + geom_density(aes(x=Br2074N_polyA)) +
  geom_density(aes(x=Br5339C1_polyA)) + geom_density(aes(x=Br5339N1_polyA)) +
  geom_density(aes(x=Br5340C1_polyA)) + geom_density(aes(x=Br5340N1_polyA)) +
  geom_density(aes(x=Br5341C1_polyA)) + geom_density(aes(x=Br5341N1_polyA)) + xlim(0,5) +
  ggtitle("Density Plot (RPKM): Downsampled PolyA") + xlab("RPKM")
  
ggplot(geneRpkm) + geom_density(aes(x=Br1113C1_RiboZero)) + ggtitle("Density Plot (RPKM): RiboZero") +
  geom_density(aes(x=Br2046C_RiboZero)) + geom_density(aes(x=Br2046N_RiboZero)) + 
  geom_density(aes(x=Br2074C_RiboZero)) + geom_density(aes(x=Br2074N_RiboZero)) +
  geom_density(aes(x=Br5339C1_RiboZero)) + geom_density(aes(x=Br5339N1_RiboZero)) +
  geom_density(aes(x=Br5340C1_RiboZero)) + geom_density(aes(x=Br5340N1_RiboZero)) +
  geom_density(aes(x=Br5341C1_RiboZero)) + geom_density(aes(x=Br5341N1_RiboZero)) + xlim(0,5)  + xlab("RPKM")

# Total reads
ggplot(d, aes(x= SampleID, y=TotalReads/1000000)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Total Reads (Million)") + 
  xlab("") +
  ggtitle("QC Metrics: Total Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Label, space = "free_x", scales = "free_x")

# Total Mapped Reads
ggplot(d, aes(x= SampleID, y=totalMapped/1000000)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Total Mapped Reads (Million)") + 
  xlab("") +
  ggtitle("QC Metrics: Total Mapped Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Label, space = "free_x", scales = "free_x")

# % Unmapped reads
ggplot(d, aes(x= SampleID, y=Unmapped)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("QC Metrics: Percent Unmapped Reads") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Label, space = "free_x", scales = "free_x")

# % Mitochondrial reads
ggplot(d, aes(x= SampleID, y=percMito)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Reads Mapping\nTo Mitochondrial Genome") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Label, space = "free_x", scales = "free_x")

# % rRNA reads
ggplot(d, aes(x= SampleID, y=percrRNA)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 70, hjust=1)) + 
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Reads Mapping\nTo Ribosomal RNA") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Label, space = "free_x", scales = "free_x")

# Mean RPKM of Expressed Genes by Group
ggplot(d, aes(x=Label, y=MeanRpkm)) + geom_boxplot() + 
  geom_jitter() + 
  ylab("Mean RPKM") + 
  xlab("") +
  ggtitle("Mean RPKM of Expressed Genes By Group") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
