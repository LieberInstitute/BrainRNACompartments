library(GenomicRanges)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(DESeq2)


load("./Dropbox/sorted_figures/github_controlled/rna_editing/data/rna_editing_results.rda")

## Redo venn diagram (Fig 4B)
editingres = editingres[c(1:6,8,10:14)]
names(editingres) = gsub(".downsampled", "", names(editingres))
editingres = Map(cbind, editingres, conversion = lapply(editingres, function(x) paste0(x$ref, ":", x$alt)),
                 sampleID = lapply(names(editingres), function(x) x), rnum = lapply(editingres, function(x) 1:nrow(x)),
                 editingID = lapply(editingres, function(x) paste0(x$chromosome,":",x$start,"-",x$end,":",x$ref,":",x$alt)))
for (i in 1:length(editingres)){
  tmp = editingres[[i]]
  tmp$Fraction = ifelse((tmp$rnum %in% grep("N", tmp$sampleID)), "Nucleus", "Cytoplasm")
  tmp$Age = ifelse((tmp$rnum %in% grep("53", tmp$sampleID)), "Prenatal", "Adult")
  tmp$Group = paste0(tmp$Age, ":", tmp$Fraction)
  tmp$collapsedconversion = NA
  tmp[which(tmp$conversion=="A:C" | tmp$conversion=="T:G"), "collapsedconversion"] = "A:C / T:G"
  tmp[which(tmp$conversion=="A:G" | tmp$conversion=="T:C"), "collapsedconversion"] = "A:G / T:C"
  tmp[which(tmp$conversion=="A:T" | tmp$conversion=="T:A"), "collapsedconversion"] = "A:T / T:A"
  tmp[which(tmp$conversion=="C:A" | tmp$conversion=="G:T"), "collapsedconversion"] = "C:A / G:T"
  tmp[which(tmp$conversion=="C:G" | tmp$conversion=="G:C"), "collapsedconversion"] = "C:G / G:C"
  tmp[which(tmp$conversion=="C:T" | tmp$conversion=="G:A"), "collapsedconversion"] = "C:T / G:A"
  editingres[[i]] = tmp
}
editingres = data.table(do.call(rbind, editingres))

editingID = editingres[collapsedconversion=="A:G / T:C",list(editingID), by = "Group"]
editingID = split(editingID, f = editingID$Group)
editingID = lapply(editingID, function(x) as.character(x$editingID))
lapply(editingID, head)
names(editingID) = gsub(":","\n", names(editingID))
editingID = list("Adult\nCytoplasm"=editingID[["Adult\nCytoplasm"]], "Prenatal\nCytoplasm"=editingID[["Prenatal\nCytoplasm"]], 
                 "Adult\nNucleus"=editingID[["Adult\nNucleus"]], "Prenatal\nNucleus"=editingID[["Prenatal\nNucleus"]])

venn.diagram(editingID, "./Dropbox/sorted_figures/github_controlled/rna_editing/figures/editing_site_overlap.jpeg", 
             main="Editing Sites Identified in All Four Groups", col = "transparent",
             fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.3,
             fill = c("#F8766D", "#7CAE00","#00BFC4","#C77CFF"), alpha = 0.50,
             cat.fontfamily = "Arial", margin=0.2)

## Make venn diagram of the unique sites

load("./Dropbox/sorted_figures/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")

unique_bySamp_all = list(cytosolOnly = unique_bySamp$cytosolOnly[-grep("no", unique_bySamp$cytosolOnly$cytosolAll),],
                         nucleusOnly = unique_bySamp$nucleusOnly[-grep("no", unique_bySamp$nucleusOnly$nucleusAll),], 
                         adultOnly = unique_bySamp$adultOnly[-grep("no", unique_bySamp$adultOnly$adultAll),], 
                         prenatalOnly = unique_bySamp$prenatalOnly[-grep("no", unique_bySamp$prenatalOnly$prenatalAll),], 
                         ANnotAC = unique_bySamp$ANnotAC[-grep("no", unique_bySamp$ANnotAC$allAN),], 
                         ACnotAN = unique_bySamp$ACnotAN[-grep("no", unique_bySamp$ACnotAN$allAC),], 
                         ANnotPN = unique_bySamp$ANnotPN[-grep("no", unique_bySamp$ANnotPN$allAN),], 
                         PNnotAN = unique_bySamp$PNnotAN[-grep("no", unique_bySamp$PNnotAN$allPN),],
                         ACnotPC = unique_bySamp$ACnotPC[-grep("no", unique_bySamp$ACnotPC$allAC),], 
                         PCnotAC = unique_bySamp$PCnotAC[-grep("no", unique_bySamp$PCnotAC$allPC),], 
                         PCnotPN = unique_bySamp$PCnotPN[-grep("no", unique_bySamp$PCnotPN$allPC),], 
                         PNnotPC = unique_bySamp$PNnotPC[-grep("no", unique_bySamp$PNnotPC$allPN),])

unique_all = lapply(unique_bySamp_all, function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = lapply(unique_all, function(x) as.character(unique(x$editingID)))
elementNROWS(unique_all)

sites = list("Adult\nCytoplasm" = c(as.character(unique_all$ACnotAN$editingID), as.character(unique_all$ACnotPC$editingID)),
             "Prenatal\nCytoplasm" = c(as.character(unique_all$PCnotPN$editingID), as.character(unique_all$PCnotAC$editingID)),
             "Adult\nNucleus" = c(as.character(unique_all$ANnotPN$editingID), as.character(unique_all$ANnotAC$editingID)),
             "Prenatal\nNucleus" = c(as.character(unique_all$PNnotAN$editingID), as.character(unique_all$PNnotPC$editingID)))
sites = lapply(sites, unique)
elementNROWS(sites)

venn.diagram(sites, "./Dropbox/sorted_figures/github_controlled/rna_editing/figures/editing_site_overlap_unique_inAll.jpeg", 
             main="Editing Sites Identified in All Four Groups", col = "transparent",
             fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.3,
             fill = c("#F8766D", "#7CAE00","#00BFC4","#C77CFF"), alpha = 0.50,
             cat.fontfamily = "Arial", margin=0.2)

## Make venn diagram of the edited genes

ov = findOverlaps(makeGRangesFromDataFrame(editingres), makeGRangesFromDataFrame(geneMap))
editingres = cbind(data.frame(editingres)[queryHits(ov),], geneMap[subjectHits(ov),])
geneID = data.table(editingres)[collapsedconversion=="A:G / T:C",list(gencodeID), by = "Group"]
geneID = split(geneID, f = geneID$Group)
geneID = lapply(geneID, function(x) as.character(x$gencodeID))
lapply(geneID, head)
names(geneID) = gsub(":","\n", names(geneID))
geneID = list("Adult\nCytoplasm"=geneID[["Adult\nCytoplasm"]], "Prenatal\nCytoplasm"=geneID[["Prenatal\nCytoplasm"]], 
                 "Adult\nNucleus"=geneID[["Adult\nNucleus"]], "Prenatal\nNucleus"=geneID[["Prenatal\nNucleus"]])
geneID = lapply(geneID, unique)

venn.diagram(geneID, "./Dropbox/sorted_figures/github_controlled/rna_editing/figures/edited_gene_overlap.jpeg", 
             main="Editing Sites Identified in All Four Groups", col = "transparent",
             fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.3,
             fill = c("#F8766D", "#7CAE00","#00BFC4","#C77CFF"), alpha = 0.50,
             cat.fontfamily = "Arial", margin=0.2)


## Redo the LFC plots as RPKM plots

load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23_nodownsamp.rda")

geneRpkm = geneRpkm[,grep("polyA", colnames(geneRpkm))]
colnames(geneRpkm) = gsub("_polyA", "", colnames(geneRpkm))
head(geneRpkm)

unique_all = lapply(unique_bySamp_all, function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = do.call(rbind, Map(cbind, unique_all[elementNROWS(unique_all)>0], group = as.list(names(unique_all[elementNROWS(unique_all)>0]))))
dim(unique_all) # 10390
dim(unique_all[which(unique_all$distToGene==0),]) # 8825
unique_all = unique_all[which(unique_all$distToGene==0),]
df = data.frame(unique_all)[,which(colnames(unique_all) %in% c("seqnames","start","end","width","strand","editingID",
                                                               "collapsedconversion","overlappingGene", "group"))] 
df = cbind(df, data.frame(geneRpkm)[match(df$overlappingGene, rownames(geneRpkm)),])
df = unique(df)
df = reshape2::melt(df, measure.vars = c("Br1113C1","Br1113N1","Br2046C","Br2046N","Br2074C","Br2074N",
                                         "Br5339C1","Br5339N1","Br5340C1","Br5340N1","Br5341C1","Br5341N1"))
df$rnum =  1:nrow(df)
df$Fraction = ifelse(df$rnum %in% grep("C", df$variable), "Cytoplasm", "Nucleus")
df$Age = ifelse(df$rnum %in% grep("Br53", df$variable), "Prenatal", "Adult")
df$group4 = paste0(df$Age, "\n", df$Fraction)
unique(df$group4)
df$group4 = factor(df$group4, levels = c("Adult\nCytoplasm","Prenatal\nCytoplasm","Adult\nNucleus","Prenatal\nNucleus"))
head(df)
df$group = gsub("adultOnly", "Adult\nOnly", df$group)
df$group = gsub("prenatalOnly", "Prenatal\nOnly", df$group)
df$group = gsub("ANnotAC", "AN\nnotAC", df$group)
df$group = gsub("ACnotAN", "AC\nnotAN", df$group)
df$group = gsub("ANnotPN", "AN\nnotPN", df$group)
df$group = gsub("PNnotAN", "PN\nnotAN", df$group)
df$group = gsub("ACnotPC", "AC\nnotPC", df$group)
df$group = gsub("PCnotAC", "PC\nnotAC", df$group)
df$group = gsub("PCnotPN", "PC\nnotPN", df$group)
df$group = gsub("PNnotPC", "PN\nnotPC", df$group)
df$group = factor(df$group, levels = c("Adult\nOnly","Prenatal\nOnly","AN\nnotAC","AC\nnotAN","PN\nnotPC","PC\nnotPN","AC\nnotPC",
                                       "PC\nnotAC","AN\nnotPN","PN\nnotAN"))
unique(df$group)
pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/unique_sites_RPKM.pdf", width=10,height=3.75)
ggplot(df, aes(x=group,y=log(value+1))) + geom_boxplot(aes(fill = group4)) + 
  scale_fill_manual(values=c("#F8766D", "#7CAE00","#00BFC4","#C77CFF")) + 
  xlab("") + ylab("log(RPKM+1)") + ggtitle("Unique Editing Sites and Gene Expression") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())
dev.off()

## Plot just the prenatal-adult only groups

head(df)
df2 = df[which(df$group %in% c("Adult\nOnly","Prenatal\nOnly")),]
df2$group = gsub("Adult\nOnly", "Adult Only", df2$group)
df2$group = gsub("Prenatal\nOnly", "Prenatal Only", df2$group)
df2$group4 = factor(df2$group4, levels = c("Adult\nCytoplasm","Adult\nNucleus","Prenatal\nCytoplasm","Prenatal\nNucleus"))

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/unique_sites_RPKM_prenataladultOnly.pdf", width=6,height=3.5)
ggplot(df2, aes(x=Fraction,y=log(value+1))) + 
  geom_boxplot(aes(fill = Age)) + facet_grid(. ~ group) +
  scale_fill_brewer(palette="Set1") + 
  xlab("") + ylab("log(RPKM+1)") + ggtitle("Unique Editing Sites and Expression") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())
dev.off()


## plot LFC differences in age groups

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

head(df2)
df2 = data.frame(df2, A.LFC = data.frame(Apres)[match(df2$overlappingGene, rownames(Apres)),"log2FoldChange"],
                 F.LFC = data.frame(Fpres.down)[match(df2$overlappingGene, rownames(Fpres.down)),"log2FoldChange"],
                 C.LFC = data.frame(Cpres.down)[match(df2$overlappingGene, rownames(Cpres.down)),"log2FoldChange"],
                 N.LFC = data.frame(Npres)[match(df2$overlappingGene, rownames(Npres)),"log2FoldChange"])
df2 = reshape2::melt(df2, measure.vars = c("A.LFC","F.LFC","C.LFC","N.LFC"))
colnames(df2)[c(10:11,16:17)] = c("sampleID","RPKM","LFC.group","LFC")
(unique(df2$LFC.group))
df2$LFC.group = gsub("A.LFC", "In Adult", df2$LFC.group)
df2$LFC.group = gsub("F.LFC", "In Prenatal", df2$LFC.group)
df2$LFC.group = gsub("C.LFC", "In Cytoplasm", df2$LFC.group)
df2$LFC.group = gsub("N.LFC", "In Nucleus", df2$LFC.group)

df3 = unique(df2[,which(colnames(df2) %in% c("editingID","overlappingGene","group","LFC.group","LFC"))])
head(df3)

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/unique_sites_LFC_byAge.pdf", width=4.3,height=3.25)
ggplot(df3[which(df3$LFC.group %in% c("In Cytoplasm","In Nucleus")),], aes(x=group,y=LFC)) + 
  geom_boxplot(aes(fill=LFC.group)) + ggtitle("Exclusively Edited Sites") +
  scale_fill_brewer(palette="Dark2") + geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Edited in:") + ylab(expression(paste(log[2], (Prenatal/Adult)))) +
  theme(title = element_text(size = 20), text = element_text(size = 18),legend.position = c(0.275, 0.8),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) +
  guides(fill=guide_legend(title="Expression:"))
dev.off()


## Plot fraction differences in RPKM

df2 = df[which(df$group %in% c("AN\nnotAC","AC\nnotAN","PC\nnotPN","PN\nnotPC")),]
df2$group4 = factor(df2$group4, levels = c("Adult\nCytoplasm","Adult\nNucleus","Prenatal\nCytoplasm","Prenatal\nNucleus"))

ggplot(df2, aes(x=Age,y=log(value+1))) + 
  geom_boxplot(aes(fill = Fraction)) + facet_grid(. ~ group) +
  scale_fill_brewer(palette="Dark2") + 
  xlab("") + ylab("log(RPKM+1)") + ggtitle("Unique Editing Sites and Expression") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), 
        legend.position="bottom", legend.title = element_blank())


## Plot fraction differences in LFC

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

head(df2)
df2 = data.frame(df2, A.LFC = data.frame(Apres)[match(df2$overlappingGene, rownames(Apres)),"log2FoldChange"],
                 F.LFC = data.frame(Fpres.down)[match(df2$overlappingGene, rownames(Fpres.down)),"log2FoldChange"],
                 C.LFC = data.frame(Cpres.down)[match(df2$overlappingGene, rownames(Cpres.down)),"log2FoldChange"],
                 N.LFC = data.frame(Npres)[match(df2$overlappingGene, rownames(Npres)),"log2FoldChange"])
df2 = reshape2::melt(df2, measure.vars = c("A.LFC","F.LFC","C.LFC","N.LFC"))
colnames(df2)[c(10:11,16:17)] = c("sampleID","RPKM","LFC.group","LFC")
(unique(df2$LFC.group))
df2$LFC.group = gsub("A.LFC", "In Adult", df2$LFC.group)
df2$LFC.group = gsub("F.LFC", "In Prenatal", df2$LFC.group)
df2$LFC.group = gsub("C.LFC", "In Cytoplasm", df2$LFC.group)
df2$LFC.group = gsub("N.LFC", "In Nucleus", df2$LFC.group)

df3 = unique(df2[,which(colnames(df2) %in% c("editingID","overlappingGene","group","LFC.group","LFC"))])
df3$Age = df3$Fraction = NA
df3[grep("N\nnot", df3$group),"Fraction"] = "Nucleus"
df3[grep("C\nnot", df3$group),"Fraction"] = "Cytoplasm"
df3[grep("notA", df3$group),"Age"] = "Adult"
df3[grep("notP", df3$group),"Age"] = "Prenatal"
df3$group = gsub("AN\nnotAC", "Nucleus", df3$group)
df3$group = gsub("PN\nnotPC", "Nucleus", df3$group)
df3$group = gsub("AC\nnotAN", "Cytoplasm", df3$group)
df3$group = gsub("PC\nnotPN", "Cytoplasm", df3$group)

head(df3)

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/unique_sites_LFC_byFraction2.pdf", width=4.3,height=3.25)
ggplot(df3[which(df3$LFC.group=="In Adult"),], aes(x=group,y=LFC)) + 
  geom_boxplot(aes(fill=Age)) + ggtitle("Exclusively Edited Sites") +
  scale_fill_brewer(palette="Set1") + geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Edited in:") + ylab(expression(paste(log[2], (Nuc./Cyt.)))) +
  theme(title = element_text(size = 20), text = element_text(size = 18),legend.position = c(0.22, 0.8),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) +
  guides(fill=guide_legend(title="Edited in:"))
ggplot(df3[which(df3$LFC.group=="In Prenatal"),], aes(x=group,y=LFC)) + 
  geom_boxplot(aes(fill=Age)) + ggtitle("Exclusively Edited Sites") +
  scale_fill_brewer(palette="Set1") + geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Edited in:") + ylab(expression(paste(log[2], (Nuc./Cyt.)))) +
  theme(title = element_text(size = 20), text = element_text(size = 18),legend.position = c(0.22, 0.8),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) +
  guides(fill=guide_legend(title="Edited in:"))
dev.off()


head(df2)
df2$LFC.group = gsub("In Adult", "Expression\nIn Adult", df2$LFC.group)
df2$LFC.group = gsub("In Prenatal", "Expression\nIn Prenatal", df2$LFC.group)
df2$LFC.group = gsub("In Cytoplasm", "Expression\nIn Cytoplasm", df2$LFC.group)
df2$LFC.group = gsub("In Nucleus", "Expression\nIn Nucleus", df2$LFC.group)

df3 = unique(df2[,which(colnames(df2) %in% c("editingID","overlappingGene","group","LFC.group","LFC"))])
df3$Age = df3$Fraction = NA
df3[grep("N\nnot", df3$group),"Fraction"] = "Nuclear-\nEdited"
df3[grep("C\nnot", df3$group),"Fraction"] = "Cytoplasm-\nEdited"
df3[grep("notA", df3$group),"Age"] = "In Adult"
df3[grep("notP", df3$group),"Age"] = "In Prenatal"
df3$group = gsub("AN\nnotAC", "Nucleus", df3$group)
df3$group = gsub("PN\nnotPC", "Nucleus", df3$group)
df3$group = gsub("AC\nnotAN", "Cytoplasm", df3$group)
df3$group = gsub("PC\nnotPN", "Cytoplasm", df3$group)

head(df3)

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/unique_sites_LFC_byFraction.pdf", width=4,height=4)
ggplot(df3[which(df3$LFC.group %in% c("Expression\nIn Adult","Expression\nIn Prenatal")),], aes(x=group,y=LFC)) + 
  geom_boxplot(aes(fill=Fraction)) + facet_grid(LFC.group ~ Age) +
  scale_fill_brewer(palette="Dark2") + geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Unique Editing Group") + ylab(expression(paste(log[2], (Nucleus/Cytoplasm)))) +
  theme(title = element_text(size = 20), text = element_text(size = 18), 
        legend.position="none", axis.text.x = element_text(angle = 15, hjust = .5, vjust=.9))
dev.off()

