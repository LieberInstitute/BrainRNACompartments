library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

jcountsP = as.data.frame(jCountsP)
jcountsR = as.data.frame(jCountsR)
jcountsD = as.data.frame(jCountsD)

match(rownames(metrics[1:24,]),
      names(c(colSums(jcountsP),colSums(jcountsR))))
df = data.frame(ID = metrics$SAMPLE_ID[1:24], 
                sum = c(colSums(jcountsP),colSums(jcountsR)),
                numMapped = metrics$numMapped[1:24])
df$propJunc = df$sum / df$numMapped
df = df[order(df$ID),]
df = df[which(df$ID!="Br1113N1_RiboZero"),]
match(df$ID, pd$SampleID)
df$Label = pd[match(df$ID, pd$SampleID),"Label"]
df$Label = gsub("Cytosol", "Cytoplasm", df$Label)
df$Label = factor(df$Label, levels = c("Prenatal\nCytoplasm\npolyA", "Prenatal\nNucleus\npolyA", "Adult\nCytoplasm\npolyA", 
                                       "Adult\nNucleus\npolyA", "Prenatal\nCytoplasm\nRiboZero", "Prenatal\nNucleus\nRiboZero",
                                       "Adult\nCytoplasm\nRiboZero", "Adult\nNucleus\nRiboZero"))
df$Library = pd[match(df$ID, pd$SampleID),"Library"]
df$Age = pd[match(df$ID, pd$SampleID),"Fetal"]
df$Age = factor(df$Age, levels = c("Prenatal", "Adult"))
df$Fraction = pd[match(df$ID, pd$SampleID),"Zone"]
df$Fraction = gsub("Cytosol", "Cytoplasm", df$Fraction)
df$Fraction = factor(df$Fraction)


pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/figures/proportion_junction_reads_full.pdf", 
    width = 8.5, height = 4)
ggplot(df, aes(x=Age, y=propJunc, fill = Fraction)) + geom_boxplot() + 
  geom_jitter() + facet_grid(. ~ Library) +
  scale_fill_brewer(palette="Dark2") +
  ylim(0, 0.4) + xlab("") +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  theme(legend.position = c(.9, 0.85)) +
  ggtitle("Proportion of Reads Overlapping Splice Junctions") + 
  theme(title = element_text(size = 20), legend.title=element_blank(), legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) +
  theme(text = element_text(size = 20))
dev.off()

sum.down = c(colSums(jcountsP),colSums(jcountsR),colSums(jcountsD))
match(rownames(metrics),names(sum.down))
df.down = data.frame(ID = metrics$SAMPLE_ID,
                     sum = sum.down,
                     numMapped = metrics$numMapped)
df.down$propJunc = df.down$sum / df.down$numMapped
df.down = df.down[which(df.down$ID!="Br5339C1_polyA" & 
                        df.down$ID!="Br5340C1_polyA" &
                        df.down$ID!="Br1113N1_RiboZero"),]
df.down$ID = gsub("downsamp", "polyA", df.down$ID)
rownames(df.down) = gsub("downsamp", "polyA", rownames(df.down))
df.down = df.down[order(df.down$ID),]
match(df.down$ID, pd$SampleID)
df.down$Label = pd[match(df.down$ID, pd$SampleID),"Label"]
df.down$Label = gsub("Cytosol", "Cytoplasm", df.down$Label)
df.down$Label = factor(df.down$Label, levels = c("Prenatal\nCytoplasm\npolyA", "Prenatal\nNucleus\npolyA", "Adult\nCytoplasm\npolyA", 
                                       "Adult\nNucleus\npolyA", "Prenatal\nCytoplasm\nRiboZero", "Prenatal\nNucleus\nRiboZero",
                                       "Adult\nCytoplasm\nRiboZero", "Adult\nNucleus\nRiboZero"))
df.down$Library = pd[match(df.down$ID, pd$SampleID),"Library"]
df.down$Age = pd[match(df.down$ID, pd$SampleID),"Fetal"]
df.down$Age = factor(df.down$Age, levels = c("Prenatal", "Adult"))
df.down$Fraction = pd[match(df.down$ID, pd$SampleID),"Zone"]
df.down$Fraction = gsub("Cytosol", "Cytoplasm", df.down$Fraction)
df.down$Fraction = factor(df.down$Fraction)

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/figures/proportion_junction_reads_downsampled.pdf", 
    width = 8.5, height = 4)
ggplot(df.down, aes(x=Age, y=propJunc, fill = Fraction)) + geom_boxplot() + 
  geom_jitter() + facet_grid(. ~ Library) +
  scale_fill_brewer(palette="Dark2") +
  ylim(0, 0.4) + xlab("") +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  theme(legend.position = c(.9, 0.85)) +
  ggtitle("Proportion of Reads Overlapping Splice Junctions") + 
  theme(title = element_text(size = 20), legend.title=element_blank(), legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) +
  theme(text = element_text(size = 20))
dev.off()


# t tests

ttests = list(junctions.frac.both = t.test(x = df.down[c(grep("Prenatal\nNucleus", df.down$Label),grep("Adult\nNucleus", df.down$Label)),"propJunc"], 
                                           y = df.down[c(grep("Prenatal\nCytoplasm", df.down$Label), grep("Adult\nCytoplasm", df.down$Label)),"propJunc"]),
              junctions.frac.polyA = t.test(x = df.down[c(grep("Prenatal\nNucleus\npolyA", df.down$Label),grep("Adult\nNucleus\npolyA", df.down$Label)), "propJunc"], 
                                            y = df.down[c(grep("Prenatal\nCytoplasm\npolyA", df.down$Label), grep("Adult\nCytoplasm\npolyA", df.down$Label)),"propJunc"]),
              junctions.frac.ribo = t.test(x = df.down[c(grep("Prenatal\nNucleus\nRiboZero", df.down$Label),grep("Adult\nNucleus\nRiboZero", df.down$Label)),"propJunc"], 
                                           y = df.down[c(grep("Prenatal\nCytoplasm\nRiboZero", df.down$Label), grep("Adult\nCytoplasm\nRiboZero", df.down$Label)),"propJunc"]),
              junctions.frac.adult.polyA = t.test(x = df.down[df.down$Label=="Adult\nNucleus\npolyA", "propJunc"], 
                                            y = df.down[df.down$Label=="Adult\nCytoplasm\npolyA","propJunc"]),
              junctions.frac.adult.ribo = t.test(x = df.down[df.down$Label=="Adult\nNucleus\nRiboZero","propJunc"], 
                                           y = df.down[df.down$Label=="Adult\nCytoplasm\nRiboZero","propJunc"]),
              junctions.frac.prenatal.polyA = t.test(x = df.down[df.down$Label=="Prenatal\nNucleus\npolyA", "propJunc"], 
                                            y = df.down[df.down$Label=="Prenatal\nCytoplasm\npolyA","propJunc"]),
              junctions.frac.prenatal.ribo = t.test(x = df.down[df.down$Label=="Prenatal\nNucleus\nRiboZero","propJunc"], 
                                           y = df.down[df.down$Label=="Prenatal\nCytoplasm\nRiboZero","propJunc"]),
              junctions.age.both = t.test(x = df.down[c(grep("Prenatal\nNucleus", df.down$Label),grep("Prenatal\nCytoplasm", df.down$Label)),"propJunc"], 
                                          y = df.down[c(grep("Adult\nNucleus", df.down$Label), grep("Adult\nCytoplasm", df.down$Label)),"propJunc"]),
              junctions.age.polyA = t.test(x = df.down[c(grep("Prenatal\nNucleus\npolyA", df.down$Label),grep("Prenatal\nCytoplasm\npolyA", df.down$Label)),"propJunc"], 
                                           y = df.down[c(grep("Adult\nNucleus\npolyA", df.down$Label), grep("Adult\nCytoplasm\npolyA", df.down$Label)),"propJunc"]),
              junctions.age.ribo = t.test(x = df.down[c(grep("Prenatal\nNucleus\nRiboZero", df.down$Label),grep("Prenatal\nCytoplasm\nRiboZero", df.down$Label)),"propJunc"], 
                                          y = df.down[c(grep("Adult\nNucleus\nRiboZero", df.down$Label), grep("Adult\nCytoplasm\nRiboZero", df.down$Label)),"propJunc"]))
df = data.frame(Group = names(ttests), t.stat = unlist(lapply(ttests, function(x) x$statistic)),
           mean.x = unlist(lapply(ttests, function(x) x$estimate[1])),
           mean.y = unlist(lapply(ttests, function(x) x$estimate[2])),
           p.value = unlist(lapply(ttests, function(x) x$p.value)))
df$FDR = p.adjust(df$p.value, method="fdr")
rownames(df) = NULL
write.csv(df, file = "/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/junction_proportions_tstat.csv")

