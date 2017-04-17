
# Disease Ontology
DO = read.csv("/Users/amandaprice/Dropbox/sorted_figures/new/GO_terms/interaction_DO.csv")
head(DO)
colnames(DO)
az = as.character(DO$geneID[1])
tao = as.character(DO$geneID[2]) # same as az
az = strsplit(as.character(az), "/", fixed=TRUE)
az = az[[1]]

res = read.csv("/Users/amandaprice/Dropbox/sorted_figures/new/interactionRes.csv")
load("/Users/amandaprice/Dropbox/sorted_figures/new/dds_interaction.rda")
az = res[which(res$Symbol %in% az),]
pdf("/Users/amandaprice/Dropbox/sorted_figures/new/AZ_gene_expression.pdf")
plots = list()
for (i in 1:nrow(az)){
  plots[[i]] = plotCounts(dds, as.character(res$EnsID[i]), 
                          intgroup = c("Fetal", "Zone", "Library"), returnData =TRUE)
  tmp = plots[[i]]
  tmp$Group = paste(tmp$Fetal,tmp$Zone, sep = "\n")
  x = ggplot(tmp, aes(x=Group, y=count)) + geom_boxplot() + 
    geom_jitter() +
    scale_y_log10() +
    ylab("Normalized Count") + 
    xlab("") +
    ggtitle(as.character(az$Symbol[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(x)
}
dev.off()



x 
dim(x)
head(geneMap)
