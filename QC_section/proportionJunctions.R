jcounts = as.data.frame(jCounts)
propJunc = colSums(jcounts) / pd$totalMapped
df = data.frame(propJunc, "Label" = pd$Label)

library(ggplot2)
ggplot(df, aes(x=Label, y=propJunc)) + geom_boxplot() + 
  geom_jitter() + 
  ylim(0, 0.4) +
  ylab("Junction Read Count/\nTotal Mapped Reads") + 
  xlab("") +
  ggtitle("Proportion of Reads Mapping to Splice Junctions") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))