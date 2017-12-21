library("readxl")
library('DESeq2')
library("ggplot2")

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")

MeanRpkm = round(c(colMeans(geneRpkm), colMeans(geneRpkm.down[,colnames(geneRpkm.down) %in% c("Br5339C1_polyA","Br5340C1_polyA")])),2)
names(MeanRpkm) = c(names(MeanRpkm[1:23]), "Br5339C1_downsamp", "Br5340C1_downsamp")
MeanRpkm = MeanRpkm[order(names(MeanRpkm))]
metrics = metrics[which(metrics$SAMPLE_ID!="Br1113N1_RiboZero"),]
metrics = metrics[order(metrics$SAMPLE_ID),]
pd = rbind(pd, pd[pd$SampleID %in% c("Br5339C1_polyA","Br5340C1_polyA"),])
rownames(pd[24:25,]) = pd$SampleID[24:25] = c("Br5339C1_downsamp","Br5340C1_downsamp")
pd = pd[order(pd$SampleID),]
pheno.table = cbind(pd, metrics[match(pd$SampleID,metrics$SAMPLE_ID),])
match(pheno.table$SampleID, pheno.table$SAMPLE_ID)
write.csv(pheno.table[,colnames(pheno.table)!="SAMPLE_ID"], 
          file="/Users/amanda/Dropbox/sorted_figures/new/github_controlled/QC_section/data/pheno.seq.data.csv",
          row.names = F, quote = F)

d = data.frame(SampleID =pheno.table$SampleID, 
               percMapped = round(pheno.table$overallMapRate * 100,2),
               Unmapped = round(100 - pheno.table$overallMapRate * 100,2),
               percrRNA = pheno.table$rRNA_rate * 100,
               percMito = pheno.table$mitoRate * 100, 
               TotalReads = pheno.table$numReads,
               totalMapped = pheno.table$numMapped,
               MeanRpkm[match(pheno.table$SampleID, names(MeanRpkm))], 
               Label = pheno.table$Label)
geneRpkm = data.frame(geneRpkm)
geneRpkm.down = data.frame(geneRpkm.down)

##Plotting by Sample

# RPKM Density Plots
pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/QC_section/figures/rpkm_density_plot_perSample.pdf")
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
dev.off()

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/QC_section/figures/QC_metrics_perSample.pdf",width=14,height=7)
# Total Reads
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
dev.off()

# Mean RPKM of Expressed Genes by Group
pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/QC_section/figures/mean_rpkm_byGroup.pdf", width=10,height=6)
ggplot(d, aes(x=Label, y=MeanRpkm)) + geom_boxplot() + 
  geom_jitter() + 
  ylab("Mean RPKM") + 
  xlab("") + 
  ggtitle("Mean RPKM of Expressed Genes By Group") + 
  theme(title = element_text(size = 20)) + 
  theme(text = element_text(size = 20))
dev.off()