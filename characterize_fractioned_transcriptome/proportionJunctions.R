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
df$Label = pd$Label
levels(df$Label) = levels = c("Adult\nCytosol\npolyA", "Prenatal\nCytosol\npolyA", "Adult\nNucleus\npolyA",
                              "Prenatal\nNucleus\npolyA", "Adult\nCytosol\nRiboZero", "Prenatal\nCytosol\nRiboZero",
                              "Adult\nNucleus\nRiboZero", "Prenatal\nNucleus\nRiboZero")


pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/figures/proportion_junction_reads_full.pdf", 
    width = 12, height = 6)
ggplot(df, aes(x=Label, y=propJunc)) + geom_boxplot() + 
  geom_jitter() + 
  ylim(0, 0.4) +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  xlab("") +
  ggtitle("Proportion of Reads Mapping to Splice Junctions") + 
  theme(title = element_text(size = 20)) +
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
df.down$Label = pd$Label
levels(df$Label) = levels = c("Adult\nCytosol\npolyA", "Prenatal\nCytosol\npolyA", "Adult\nNucleus\npolyA",
                              "Prenatal\nNucleus\npolyA", "Adult\nCytosol\nRiboZero", "Prenatal\nCytosol\nRiboZero",
                              "Adult\nNucleus\nRiboZero", "Prenatal\nNucleus\nRiboZero")


pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/figures/proportion_junction_reads_downsampled.pdf", 
    width = 12, height = 6)
ggplot(df.down, aes(x=Label, y=propJunc)) + geom_boxplot() + 
  geom_jitter() + 
  ylim(0, 0.4) +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  xlab("") +
  ggtitle("Proportion of Reads Mapping to Splice Junctions") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# t tests

ttests = list(junctions.frac.both = t.test(x = df.down[c(grep("Prenatal\nNucleus", df.down$Label),grep("Adult\nNucleus", df.down$Label)),"propJunc"], 
                                           y = df.down[c(grep("Prenatal\nCytosol", df.down$Label), grep("Adult\nCytosol", df.down$Label)),"propJunc"]),
              junctions.frac.polyA = t.test(x = df.down[c(grep("Prenatal\nNucleus\npolyA", df.down$Label),grep("Adult\nNucleus\npolyA", df.down$Label)), "propJunc"], 
                                            y = df.down[c(grep("Prenatal\nCytosol\npolyA", df.down$Label), grep("Adult\nCytosol\npolyA", df.down$Label)),"propJunc"]),
              junctions.frac.ribo = t.test(x = df.down[c(grep("Prenatal\nNucleus\nRiboZero", df.down$Label),grep("Adult\nNucleus\nRiboZero", df.down$Label)),"propJunc"], 
                                           y = df.down[c(grep("Prenatal\nCytosol\nRiboZero", df.down$Label), grep("Adult\nCytosol\nRiboZero", df.down$Label)),"propJunc"]),
              junctions.frac.adult.polyA = t.test(x = df.down[df.down$Label=="Adult\nNucleus\npolyA", "propJunc"], 
                                            y = df.down[df.down$Label=="Adult\nCytosol\npolyA","propJunc"]),
              junctions.frac.adult.ribo = t.test(x = df.down[df.down$Label=="Adult\nNucleus\nRiboZero","propJunc"], 
                                           y = df.down[df.down$Label=="Adult\nCytosol\nRiboZero","propJunc"]),
              junctions.frac.prenatal.polyA = t.test(x = df.down[df.down$Label=="Prenatal\nNucleus\npolyA", "propJunc"], 
                                            y = df.down[df.down$Label=="Prenatal\nCytosol\npolyA","propJunc"]),
              junctions.frac.prenatal.ribo = t.test(x = df.down[df.down$Label=="Prenatal\nNucleus\nRiboZero","propJunc"], 
                                           y = df.down[df.down$Label=="Prenatal\nCytosol\nRiboZero","propJunc"]),
              junctions.age.both = t.test(x = df.down[c(grep("Prenatal\nNucleus", df.down$Label),grep("Prenatal\nCytosol", df.down$Label)),"propJunc"], 
                                          y = df.down[c(grep("Adult\nNucleus", df.down$Label), grep("Adult\nCytosol", df.down$Label)),"propJunc"]),
              junctions.age.polyA = t.test(x = df.down[c(grep("Prenatal\nNucleus\npolyA", df.down$Label),grep("Prenatal\nCytosol\npolyA", df.down$Label)),"propJunc"], 
                                           y = df.down[c(grep("Adult\nNucleus\npolyA", df.down$Label), grep("Adult\nCytosol\npolyA", df.down$Label)),"propJunc"]),
              junctions.age.ribo = t.test(x = df.down[c(grep("Prenatal\nNucleus\nRiboZero", df.down$Label),grep("Prenatal\nCytosol\nRiboZero", df.down$Label)),"propJunc"], 
                                          y = df.down[c(grep("Adult\nNucleus\nRiboZero", df.down$Label), grep("Adult\nCytosol\nRiboZero", df.down$Label)),"propJunc"]))
df = rbind(p.value = unlist(lapply(ttests, function(x) x$p.value)),
           t.stat = unlist(lapply(ttests, function(x) x$statistic)),
           mean.x = unlist(lapply(ttests, function(x) x$estimate[1])),
           mean.y = unlist(lapply(ttests, function(x) x$estimate[2])))
write.csv(df, file = "/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/junction_proportions_tstat.csv")

