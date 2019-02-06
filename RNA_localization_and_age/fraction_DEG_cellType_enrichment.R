library(ggplot2)
library(GenomicRanges)
library(data.table)
library(RColorBrewer)


load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/rawCounts_darmanisSingleCell.rda", verbose=T)
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/rpkmCounts_darmanisSingleCell.rda", verbose=T)
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


## Compare expression between all cell types

library(DESeq2)

geneCounts = geneCounts[which(rowSums(geneCounts) > 0),]
pd$isNeuron = ifelse(pd$Cell_type =="Neurons", "Neuron", "nonNeuron")
pd$isAstrocytes = ifelse(pd$Cell_type =="Astrocytes", "Astrocytes", "nonAstrocytes")
pd$isHybrid = ifelse(pd$Cell_type =="Hybrid", "Hybrid", "nonHybrid")
pd$isMicroglia = ifelse(pd$Cell_type =="Microglia", "Microglia", "nonMicroglia")
pd$isOligodendrocytes = ifelse(pd$Cell_type =="Oligodendrocytes", "Oligodendrocytes", "nonOligodendrocytes")
pd$isOPC = ifelse(pd$Cell_type =="OPC", "OPC", "nonOPC")
pd$isEndothelial = ifelse(pd$Cell_type =="Endothelial", "Endothelial", "nonEndothelial")

cellTypes = c("Neurons", "Astrocytes","Hybrid","Microglia","Oligodendrocytes","Fetal (quiescent)","OPC","Endothelial","Fetal (replicating)")

dds.all = list(Neurons = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))], 
                                                colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isNeuron),
               Astrocytes = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))], 
                                                   colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isAstrocytes),
               Hybrid = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))],
                                               colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isHybrid),
               Microglia = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))],
                                                  colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isMicroglia),
               Oligodendrocytes = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))],
                                                         colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isOligodendrocytes),
               OPC = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))],
                                            colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isOPC),
               Endothelial = DESeqDataSetFromMatrix(countData = geneCounts[,-which(colnames(geneCounts) %in% unique(pd[which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"]))],
                                                    colData = pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),], design = ~ isEndothelial))

dds.all = lapply(dds.all, DESeq)
res = lapply(dds.all, results)
res = lapply(res, data.frame)
res = Map(cbind, res, lapply(res, function(x) geneMap[match(rownames(x), rownames(geneMap)),]))
save(res, dds.all, geneMap, pd, file = "/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/Darmanis_DEGs_byCellType.rda")


nup170 = do.call(rbind, Map(cbind, Comparison = as.list(names(res)), lapply(res, function(x) x[which(x$Symbol %in% c("NUP133","NUP107")),])))
nup170[which(nup170$padj<=0.05),]
#        Comparison  baseMean log2FoldChange    lfcSE     stat        pvalue         padj
#NUP133   Microglia 118.59605       6.515826 1.169095 5.573391  2.498278e-08 1.244803e-07
#NUP107   Microglia  67.73954       5.905580 1.422545 4.151419  3.304198e-05 1.079473e-04
#NUP133 Endothelial 117.50304       5.259165 1.042857 5.043037  4.581996e-07 3.821288e-06


genes = do.call(rbind, Map(cbind, CellType = as.list(names(res)), 
                           lapply(res, function(x) cbind(x[which(x$Symbol %in% c("MVP","RANGAP1","TMEM33","EIF5A2","SENP2","RANBP3L","MX2","SEH1L","NPIPA1")),],
                                                              ensID = rownames(x[which(x$Symbol %in% c("MVP","RANGAP1","TMEM33","EIF5A2","SENP2","RANBP3L","MX2","SEH1L","NPIPA1")),])))))
genes = cbind(genes, geneRpkm[match(genes$ensID, rownames(geneRpkm)), which(colnames(geneRpkm) %in% pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"])])

genes = reshape2::melt(genes, measure.vars = colnames(genes)[grep("SRR", colnames(genes))])
genes$Cell_type = pd[match(genes$variable, rownames(pd)),"Cell_type"]
save(genes,file = "/dcl01/lieber/ajaffe/Amanda/NPCgenes.rda")
load("./Desktop/BAMS/NPCgenes.rda")

brewer.pal(n = 9, name = "Accent")
# "#FDC086" "#FFFF99" were the fetal cell types earlier

#pdf("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/RNA_localization_and_age/NPC_increasing_genes_byCellType.pdf", width = 13.5, height=9)
pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/NPC_increasing_genes_byCellType.pdf", width = 13.5, height=9)

ggplot(genes[which(genes$CellType=="Neurons"),], aes(x = Cell_type, y = log(value+1), colour = Cell_type)) + geom_jitter() + geom_boxplot() +
  facet_wrap(Symbol ~ .) +
  scale_color_manual(values=c("#7FC97F","#BEAED4","#386CB0","#F0027F","#BF5B17","#666666","ivory3")) +
  labs(fill="") + ylab("log(RPKM+1)") + xlab("") +
  ggtitle("NPC Genes with Increasing Expression") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = .8))
dev.off()


npc170 = data.frame(t(geneRpkm[which(rownames(geneRpkm) %in% rownames(geneMap[which(geneMap$Symbol %in% c("NUP133","NUP107")),])), 
                               which(colnames(geneRpkm) %in% pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"])]),
                    ID = rownames(t(geneRpkm[which(rownames(geneRpkm) %in% rownames(geneMap[which(geneMap$Symbol %in% c("NUP133","NUP107")),])), 
                                             which(colnames(geneRpkm) %in% pd[-which(pd$Cell_type %in% c("Fetal_quiescent","Fetal_replicating")),"RunName"])])))
npc170 = reshape2::melt(npc170, measure.vars = colnames(npc170)[grep("ENSG", colnames(npc170))])
npc170$Cell_type = pd[match(npc170$ID, rownames(pd)),"Cell_type"]
npc170$Symbol = geneMap[match(npc$variable, rownames(geneMap)),"Symbol"]


pdf("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/RNA_localization_and_age/Nup133_Nup107.pdf", width = 8.3, height=3.5)
ggplot(npc170, aes(x = Cell_type, y = log(value+1), colour = Cell_type)) + geom_jitter() + geom_boxplot() +
  facet_grid(. ~ Symbol) +
  scale_color_manual(values=c("#7FC97F","#BEAED4","#386CB0","#F0027F","#BF5B17","#666666","ivory3")) +
  labs(fill="") + ylab("log(RPKM+1)") + xlab("") +
  ggtitle("Nup107-160 Nucleoporin Subcomplex Expression") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = .8))
dev.off()


## which genes are statistically significantly different?

unique(genes[which(genes$padj<=0.05),1:15])
#     CellComparison   baseMean log2FoldChange     lfcSE      stat       pvalue         padj   Chr     Start       End Strand Length  Symbol  EntrezID           ensID
#4           Neurons 381.834093      5.2191854 0.6583361  7.927843 2.229855e-15 3.295761e-13  chr5  36248536  36302216      -   3363 RANBP3L    202151 ENSG00000164188
#5           Neurons  24.477348     -0.9702069 0.3230074 -3.003668 2.667466e-03 1.426248e-02 chr16  15016659  15045915      +   7275  NPIPA1 101059953 ENSG00000183426
#8           Neurons  46.051886     -2.0328962 0.5551573 -3.661838 2.504124e-04 2.111105e-03 chr22  41641615  41682255      -   5463 RANGAP1      5905 ENSG00000100401
#13       Astrocytes 390.094352     -3.7239811 0.8428475 -4.418333 9.946492e-06 8.704158e-05  chr5  36248536  36302216      -   3363 RANBP3L    202151 ENSG00000164188
#16       Astrocytes  96.061878      1.9781513 0.6270263  3.154814 1.606005e-03 7.461328e-03 chr18  12947132  12987535      +  10481   SEH1L     81929 ENSG00000085415
#17       Astrocytes  32.778991      2.8063901 0.6866232  4.087235 4.365451e-05 3.276295e-04 chr22  41641615  41682255      -   5463 RANGAP1      5905 ENSG00000100401
#28        Microglia  63.361382     30.0000000 1.7408324 17.233135 1.498070e-66 7.550987e-65  chr3 170606204 170626482      -   5593  EIF5A2     56648 ENSG00000163577
#29        Microglia 100.601127      6.9683288 0.9838834  7.082474 1.416032e-12 1.182540e-11  chr3 185300284 185351339      +   7169   SENP2     59343 ENSG00000163904
#30        Microglia 939.511061      2.9043486 0.8689623  3.342318 8.308171e-04 2.173084e-03  chr4  41937137  41962589      +   8240  TMEM33     55161 ENSG00000109133
#31        Microglia 269.746406     25.8526548 1.6282947 15.877135 9.125202e-57 3.115815e-55  chr5  36248536  36302216      -   3363 RANBP3L    202151 ENSG00000164188
#32        Microglia  24.040649      6.6286831 0.8545420  7.757001 8.696102e-15 8.012905e-14 chr16  15016659  15045915      +   7275  NPIPA1 101059953 ENSG00000183426
#33        Microglia   4.793893      4.3983283 1.8086850  2.431782 1.502473e-02 2.924814e-02 chr16  29831715  29859355      +   6136     MVP      9961 ENSG00000013364
#34        Microglia 110.207246      3.5910341 1.1592202  3.097801 1.949620e-03 4.703527e-03 chr18  12947132  12987535      +  10481   SEH1L     81929 ENSG00000085415
#35        Microglia  32.598035      9.7961564 1.3293612  7.369070 1.718219e-13 1.496132e-12 chr22  41641615  41682255      -   5463 RANGAP1      5905 ENSG00000100401
#37 Oligodendrocytes  62.162996     -3.2006445 1.1438239 -2.798197 5.138876e-03 1.481298e-02  chr3 170606204 170626482      -   5593  EIF5A2     56648 ENSG00000163577
#40 Oligodendrocytes 304.050938    -10.4128870 1.0625629 -9.799784 1.128268e-22 9.915417e-21  chr5  36248536  36302216      -   3363 RANBP3L    202151 ENSG00000164188
#41 Oligodendrocytes  24.477348     -3.8577003 0.4976638 -7.751620 9.072774e-15 2.762087e-13 chr16  15016659  15045915      +   7275  NPIPA1 101059953 ENSG00000183426
#42 Oligodendrocytes   4.400168     -4.2613987 1.1966650 -3.561062 3.693573e-04 1.548246e-03 chr16  29831715  29859355      +   6136     MVP      9961 ENSG00000013364
#44 Oligodendrocytes  35.738189     -4.9231014 0.8471738 -5.811206 6.202434e-09 7.568497e-08 chr22  41641615  41682255      -   5463 RANGAP1      5905 ENSG00000100401
#46              OPC  63.361382     -9.4193009 1.6474679 -5.717441 1.081401e-08 4.628003e-07  chr3 170606204 170626482      -   5593  EIF5A2     56648 ENSG00000163577
#49              OPC 269.746406     -8.8938096 1.5296358 -5.814331 6.087669e-09 2.753260e-07  chr5  36248536  36302216      -   3363 RANBP3L    202151 ENSG00000164188
#50              OPC  24.477348     -3.0466312 0.7204356 -4.228874 2.348638e-05 5.089982e-04 chr16  15016659  15045915      +   7275  NPIPA1 101059953 ENSG00000183426
#53              OPC  32.598035     -8.4173120 1.2578827 -6.691651 2.206669e-11 1.418292e-09 chr22  41641615  41682255      -   5463 RANGAP1      5905 ENSG00000100401
#56      Endothelial  99.589212      2.6214399 0.8524658  3.075126 2.104136e-03 7.121264e-03  chr3 185300284 185351339      +   7169   SENP2     59343 ENSG00000163904
#59      Endothelial  24.319178      3.3739888 0.6928225  4.869918 1.116444e-06 8.534581e-06 chr16  15016659  15045915      +   7275  NPIPA1 101059953 ENSG00000183426
#61      Endothelial 108.611379      2.6452511 1.0435116  2.534951 1.124630e-02 2.841769e-02 chr18  12947132  12987535      +  10481   SEH1L     81929 ENSG00000085415


## which is the numerator?

#Wald test p-value: isNeuron nonNeuron vs Neuron 
#Wald test p-value: isAstrocytes nonAstrocytes vs Astrocytes
#Wald test p-value: isHybrid nonHybrid vs Hybrid 
#Wald test p-value: isMicroglia nonMicroglia vs Microglia
#Wald test p-value: isOligodendrocytes Oligodendrocytes vs nonOligodendrocytes 
#Wald test p-value: isOPC OPC vs nonOPC
#Wald test p-value: isEndothelial nonEndothelial vs Endothelial

sig = unique(genes[which(genes$padj<=0.05),1:15])
rbind(sig[which(sig$CellType %in% c("Oligodendrocytes","OPC") & sig$log2FoldChange>0),],
      sig[which(sig$CellType!="Oligodendrocytes" & sig$CellType!="OPC" & sig$log2FoldChange<0),]) 

#     CellType  baseMean log2FoldChange     lfcSE      stat       pvalue         padj   Chr    Start      End Strand Length  Symbol  EntrezID           ensID
#5     Neurons  24.47735     -0.9702069 0.3230074 -3.003668 2.667466e-03 1.426248e-02 chr16 15016659 15045915      +   7275  NPIPA1 101059953 ENSG00000183426
#8     Neurons  46.05189     -2.0328962 0.5551573 -3.661838 2.504124e-04 2.111105e-03 chr22 41641615 41682255      -   5463 RANGAP1      5905 ENSG00000100401
#13 Astrocytes 390.09435     -3.7239811 0.8428475 -4.418333 9.946492e-06 8.704158e-05  chr5 36248536 36302216      -   3363 RANBP3L    202151 ENSG00000164188

load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Darmanis_DEGs_byCellType.rda", verbose=T)
head(geneMap)
rownames(geneMap[which(geneMap$Symbol %in% c("NUP133","NUP107")),])
NPCsubcomplex = do.call(rbind, lapply(res, function(x) x[which(rownames(x) %in% rownames(geneMap[which(geneMap$Symbol %in% c("NUP133","NUP107")),])),]))
NPCsubcomplex[which(NPCsubcomplex$padj<=0.05),]
