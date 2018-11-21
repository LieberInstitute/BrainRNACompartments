library(ggplot2)
library(GenomicRanges)
library(data.table)


load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/rawCounts_darmanisSingleCell.rda")
load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/retained.byAge.downsampled.rda")

## Compare expression in Darmanis data to see what cell types may be contributing to this effect

gc = geneCounts/(colSums(geneCounts)/1000000)

gcFrac = lapply(sig[elementNROWS(sig)>0], function(x) data.frame(gc[which(rownames(gc) %in% x$ensID),]))
gcFrac = do.call(rbind, Map(cbind, gcFrac[names(gcFrac)!="ret_Ad_exp_Fet"], 
                               fracReg = as.list(c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nPrenatal Only", "Nuclear:\nAdult Only", 
                                                   "Cytoplasmic:\nPrenatal Only", "Cytoplasmic:\nAdult Only", "Interaction"))))
gcFrac$fracReg = factor(gcFrac$fracReg, 
                           levels = c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nAdult Only", "Nuclear:\nPrenatal Only",
                                      "Cytoplasmic:\nAdult Only","Cytoplasmic:\nPrenatal Only", "Interaction"))
gcFrac$geneID = rownames(gcFrac)
gcFrac$geneID = gsub("^.*\\.","", gcFrac$geneID)
gcFrac = reshape2::melt(gcFrac)
gcFrac$AgeGroup = pd[match(gcFrac$variable, rownames(pd)),"AgeGroup"]
gcFrac$Cell_type = pd[match(gcFrac$variable, rownames(pd)),"Cell_type"]



## Plot the cell type specific expression patterns in fraction genes

load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/CTS_expression_fractionDEGs.rda")


pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/CTS_expression_fractionDEGs.pdf", width=8.5, height=4.5)
for (i in 1:length(unique(gcFrac$fracReg))) {
  g = ggplot(gcFrac[which(gcFrac$fracReg==unique(gcFrac$fracReg)[i]),], 
             aes(x = Cell_type, y = value, fill=Cell_type)) + geom_boxplot() +
    scale_fill_brewer(palette="Accent") + 
    labs(fill="") + ylab("RPM") + xlab("") +
    ggtitle(paste0("Cell Type Expression of Fraction-Regulated Genes\n", unique(gcFrac$fracReg)[i])) +
    theme(title = element_text(size = 20), text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
}
dev.off()


dt = data.table(gcFrac)
dt = dt[dt[, .I[value == max(value)], by = c("geneID", "fracReg")]$V1]
dt = dt[,length(unique(geneID)), by = c("Cell_type", "fracReg")]
dt$fracReg = gsub("Nuclear: Both", "Nuclear:\nBoth", dt$fracReg)
dt$fracReg = gsub("Cytoplasmic: Both", "Cytoplasmic:\nBoth", dt$fracReg)
dt$Cell_type = gsub("Fetal_replicating", "Fetal (replicating)", dt$Cell_type)
dt$Cell_type = gsub("Fetal_quiescent", "Fetal (quiescent)", dt$Cell_type)
dt$fracReg = factor(dt$fracReg, levels = c("Nuclear:\nBoth", "Cytoplasmic:\nBoth", "Nuclear:\nAdult Only", "Nuclear:\nPrenatal Only",
                                           "Cytoplasmic:\nAdult Only","Cytoplasmic:\nPrenatal Only", "Interaction"))
dt$Cell_type = factor(dt$Cell_type, levels = c("Neurons", "Oligodendrocytes", "OPC", "Astrocytes", 
                                               "Microglia", "Endothelial", "Fetal (quiescent)", "Fetal (replicating)", "Hybrid"))
dt = dt[which(dt$fracReg!="Interaction"),,]

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/CTS_expression_max.pdf", width=11.5, height=3.75)
ggplot(dt, aes(x = fracReg, y = V1, fill=Cell_type)) + geom_bar(position = "fill",stat = "identity") +
    scale_fill_brewer(palette="Set1") + 
    labs(fill="") + ylab("Proportion") + xlab("") +
    ggtitle("Cell Type With Highest Expression") +
    theme(title = element_text(size = 20), text = element_text(size = 20))
dev.off()


dt = gcFrac[which(gcFrac$fracReg!="Interaction"),]
dt$Dir = NA
dt[grep("Nuclear", dt$fracReg),"Dir"] = "Nuclear"
dt[grep("Cytoplasmic", dt$fracReg),"Dir"] = "Cytoplasmic"
dt = data.table(dt)

dt = dt[dt[, .I[value == max(value)], by = c("geneID", "Dir")]$V1]
dt = dt[,length(unique(geneID)), by = c("Cell_type", "Dir")]
dt$Cell_type = gsub("Fetal_replicating", "Fetal (replicating)", dt$Cell_type)
dt$Cell_type = gsub("Fetal_quiescent", "Fetal (quiescent)", dt$Cell_type)
dt$Cell_type = factor(dt$Cell_type, levels = c("Neurons", "Oligodendrocytes", "OPC", "Astrocytes", 
                                               "Microglia", "Endothelial", "Fetal (quiescent)", "Fetal (replicating)", "Hybrid"))

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/CTS_expression_max_byDirOnly.pdf", width=6, height=3.75)
ggplot(dt, aes(x = Dir, y = V1, fill=Cell_type)) + geom_bar(position = "fill",stat = "identity") +
  scale_fill_brewer(palette="Set1") + 
  labs(fill="") + ylab("Proportion") + xlab("") +
  ggtitle("Cell Type With Highest Expression") +
  theme(title = element_text(size = 20), text = element_text(size = 20))
dev.off()

dt[,sum(V1), by=Dir]
#           Dir   V1
#1:     Nuclear 2369
#2: Cytoplasmic 2650

dt$perc = dt$V1 / c(rep.int(2369, 9), rep.int(2650, 9)) * 100
dt
#Cell_type         Dir  V1      perc
# 1:             Neurons     Nuclear 693 29.252849
# 2:          Astrocytes     Nuclear 428 18.066695
# 3:              Hybrid     Nuclear 248 10.468552
# 4:           Microglia     Nuclear  94  3.967919
# 5:    Oligodendrocytes     Nuclear 189  7.978050
# 6:   Fetal (quiescent)     Nuclear 330 13.929928
# 7:                 OPC     Nuclear  75  3.165893
# 8:         Endothelial     Nuclear 145  6.120726
# 9: Fetal (replicating)     Nuclear 167  7.049388
#10:    Oligodendrocytes Cytoplasmic 260  9.811321
#11:   Fetal (quiescent) Cytoplasmic 341 12.867925
#12:             Neurons Cytoplasmic 686 25.886792
#13:              Hybrid Cytoplasmic 338 12.754717
#14:           Microglia Cytoplasmic 119  4.490566
#15:                 OPC Cytoplasmic  76  2.867925
#16:          Astrocytes Cytoplasmic 521 19.660377
#17: Fetal (replicating) Cytoplasmic 188  7.094340
#18:         Endothelial Cytoplasmic 121  4.566038

