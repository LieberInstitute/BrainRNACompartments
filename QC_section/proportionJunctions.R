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

ggplot(df, aes(x=Label, y=propJunc)) + geom_boxplot() + 
  geom_jitter() + 
  ylim(0, 0.4) +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  xlab("") +
  ggtitle("Proportion of Reads Mapping to Splice Junctions") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

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

ggplot(df.down, aes(x=Label, y=propJunc)) + geom_boxplot() + 
  geom_jitter() + 
  ylim(0, 0.4) +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  xlab("") +
  ggtitle("Proportion of Reads Mapping to Splice Junctions") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
