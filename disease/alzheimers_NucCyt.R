library(clusterProfiler)
library(GenomicRanges)
library(ggplot2)
library(reshape2)


load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.polyAonly.sig1.downsampled.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23_nodownsamp.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Disease Ontology
plot(compareDO.noLFC,colorBy="p.adjust",  showCategory = 45, title= "Disease Ontology Enrichment")
DO = as.data.frame(compareDO)

disease = list(az = as.character(DO[which(DO$Description=="Alzheimer's disease" & DO$Cluster=="Interaction"),"geneID"]),
               ls = as.character(DO[which(DO$Description=="lateral sclerosis" & DO$Cluster=="Interaction"),"geneID"]))

disease = lapply(disease, function(x) strsplit(as.character(x), "/", fixed=TRUE))
disease = unlist(disease, recursive=F)
elementNROWS(disease)

disease = lapply(disease, function(x) geneMap[which(geneMap$EntrezID %in% x),"gencodeID"])

# make rpkm object for ggplot
azCounts = geneRpkm.down[rownames(geneRpkm.down) %in% disease$az,grep("poly",colnames(geneRpkm.down))] 
x = t(azCounts)
x = data.frame(x)
match(rownames(x), pd$SampleID)
x$Age = pd[match(rownames(x), pd$SampleID),"Fetal"] 
x$Fraction = pd[match(rownames(x), pd$SampleID),"Zone"]
x = melt(x)
x$sym = as.character(geneMap[match(x$variable, geneMap$gencodeID),"Symbol"])
x$Age = gsub("Prenatal", "P", x$Age)
x$Age = gsub("Adult", "A", x$Age)
head(x)
x[x$variable=="ENSG00000126767.17_1",]

# plot AZ genes

pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/AZ_rpkm_plots.pdf", width=4, height=3)
for (i in 1:length(unique(x$sym))) {
g = ggplot(x[x$sym == unique(x$sym)[i],], aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") +
  ylab("Log(RPKM+1)") + xlab("") +
  ggtitle(unique(x$sym)[i]) + 
  theme(title = element_text(size = 16)) + theme(text = element_text(size = 16)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
print(g)
}
dev.off()

## Plot interaction genes
int = DO[which(DO$Cluster=="Interaction" & DO$Description!="pleural disease"
               & DO$Description!="prion disease" & DO$Description!="tauopathy"),]
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/DO_interactionOnly.pdf", height = 3, width = 5)
ggplot(int, aes(x = Description, y = -log(p.adjust))) + 
  geom_bar(stat = "identity") + scale_fill_brewer(palette="Dark2") +
  geom_text(data=int, aes(x = Description, y = -log(p.adjust), label = GeneRatio), 
            size=4, nudge_y = -1, color="white") +
  coord_flip() + labs(fill="") + ylab("-log(Adjusted P-Value)") + xlab("") +
  ggtitle("Disease Ontology") +
  theme(title = element_text(size = 16)) +
  theme(text = element_text(size = 16))
dev.off()


### plot expression of disease-associated genes
# Prader-Willi
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/NPAP1_prader-willi_byFraction_byAge_byLibrary.pdf", width = 7, height = 5)
ggplot(x, aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=ENSG00000185823.3_1)) +
  facet_grid(. ~ Library) + scale_fill_brewer(palette="Dark2") +
  labs(fill="") +
  ylab("log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("NPAP1 (ENSG00000185823.3_1)\nFDR (Adult)=",round(Ares[which(rownames(Ares)=="ENSG00000185823.3_1"),"padj"], digits=5),
                 "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)=="ENSG00000185823.3_1"),"padj"], digits=5))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/NPAP1_prader-willi_byFraction_byAge.pdf", width = 7, height = 5)
ggplot(x, aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=ENSG00000185823.3_1)) +
  labs(fill="") + scale_fill_brewer(palette="Dark2") +
  ylab("log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("NPAP1 (ENSG00000185823.3_1)\nFDR (Adult)=",round(Ares[which(rownames(Ares)=="ENSG00000185823.3_1"),"padj"], digits=5),
                 "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)=="ENSG00000185823.3_1"),"padj"], digits=5))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# APOE
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/APOE_byFraction_byAge_byLibrary.pdf", width = 7, height = 5)
ggplot(x, aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=ENSG00000130203.9_2)) +
  facet_grid(. ~ Library) + scale_fill_brewer(palette="Dark2") +
  labs(fill="") +
  ylab("log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("APOE (ENSG00000130203.9_2)\nFDR (Adult)=",round(Ares[which(rownames(Ares)=="ENSG00000130203.9_2"),"padj"], digits=5),
                 "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)=="ENSG00000130203.9_2"),"padj"], digits=5))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/APOE_byFraction_byAge.pdf", width = 7, height = 5)
ggplot(x, aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=ENSG00000130203.9_2)) +
  labs(fill="") + scale_fill_brewer(palette="Dark2") +
  ylab("log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("APOE (ENSG00000130203.9_2)\nFDR (Adult)=",round(Ares[which(rownames(Ares)=="ENSG00000130203.9_2"),"padj"], digits=5),
                 "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)=="ENSG00000130203.9_2"),"padj"], digits=5))) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

## Alzheimers
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/AZ_gene_expression_byAge_byFraction_byLibrary.pdf", width = 7, height = 5)
plots = list()
for (i in 1:length(disease[["az"]])){
  print(i)
  plots[[i]] = ggplot(diseaseRPKM[["az"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["az"]][i]))) +
    facet_grid(. ~ Library) + scale_fill_brewer(palette="Dark2") +
    labs(fill="") +
    ylab("log(RPKM+1)") + 
    xlab("") +
    ggtitle(paste0(sym[["az"]][i]," (", disease[["az"]][i], ")\nFDR=",round(Ires.down[which(rownames(Ires.down)==disease[["az"]][i]),"padj"], digits=5))) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(plots[[i]])
}
dev.off()   

pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/AZ_gene_expression_byAge_byFraction.pdf", width = 7, height = 5)
for (i in 1:length(disease[["az"]])){
 print(i)
 plots[[i]] = ggplot(diseaseRPKM[["az"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["az"]][i]))) +
   labs(fill="") + scale_fill_brewer(palette="Dark2") +
   ylab("log(RPKM+1)") +
   xlab("") +
   ggtitle(paste0(sym[["az"]][i]," (", disease[["az"]][i], ")\nFDR=",round(Ires.down[which(rownames(Ires.down)==disease[["az"]][i]),"padj"], digits=5))) +
   theme(title = element_text(size = 20)) +
   theme(text = element_text(size = 20))
 print(plots[[i]])
 }
 dev.off()   
 
 ## Motor neuron disease
 pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/motorNeuron_gene_expression_byAge_byFraction_byLibrary.pdf", width = 7, height = 5)
 for (i in 1:length(disease[["motor"]])){
   print(i)
   plots[[i]] = ggplot(diseaseRPKM[["motor"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["motor"]][i]))) +
     facet_grid(. ~ Library) +
     labs(fill="") + scale_fill_brewer(palette="Dark2") +
     ylab("log(RPKM+1)") + 
     xlab("") +
     ggtitle(paste0(sym[["motor"]][i]," (", disease[["motor"]][i], ")\nFDR (Adult)=",round(Ares[which(rownames(Ares)==disease[["motor"]][i]),"padj"], digits=5),
                    "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)==disease[["motor"]][i]),"padj"], digits=5))) +
     theme(title = element_text(size = 20)) +
     theme(text = element_text(size = 20))
   print(plots[[i]])
 }
 dev.off()   
 
 pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/motorNeuron_gene_expression_byAge_byFraction.pdf", width = 7, height = 5)
 for (i in 1:length(disease[["motor"]])){
   print(i)
   plots[[i]] = ggplot(diseaseRPKM[["motor"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["motor"]][i]))) +
     labs(fill="") + scale_fill_brewer(palette="Dark2") +
     ylab("log(RPKM+1)") +
     xlab("") +
     ggtitle(paste0(sym[["motor"]][i]," (", disease[["motor"]][i], ")\nFDR (Adult)=",round(Ares[which(rownames(Ares)==disease[["motor"]][i]),"padj"], digits=5),
                    "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)==disease[["motor"]][i]),"padj"], digits=5))) +
     theme(title = element_text(size = 20)) +
     theme(text = element_text(size = 20))
   print(plots[[i]])
 }
 dev.off() 

## ALS
 pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/ALS_gene_expression_byAge_byFraction_byLibrary.pdf", width = 7, height = 5)
 for (i in 1:length(disease[["als"]])){
   print(i)
   plots[[i]] = ggplot(diseaseRPKM[["als"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["als"]][i]))) +
     facet_grid(. ~ Library) +
     labs(fill="") + scale_fill_brewer(palette="Dark2") +
     ylab("log(RPKM+1)") + 
     xlab("") +
     ggtitle(paste0(sym[["als"]][i]," (", disease[["als"]][i], ")\nFDR (Adult)=",round(Ares[which(rownames(Ares)==disease[["als"]][i]),"padj"], digits=5),
                    "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)==disease[["als"]][i]),"padj"], digits=5))) +
     theme(title = element_text(size = 20)) +
     theme(text = element_text(size = 20))
   print(plots[[i]])
 }
 dev.off()   
 
 pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/ALS_gene_expression_byAge_byFraction.pdf", width = 7, height = 5)
 for (i in 1:length(disease[["als"]])){
   print(i)
   plots[[i]] = ggplot(diseaseRPKM[["als"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["als"]][i]))) +
     labs(fill="") + scale_fill_brewer(palette="Dark2") +
     ylab("log(RPKM+1)") +
     xlab("") +
     ggtitle(paste0(sym[["als"]][i]," (", disease[["als"]][i], ")\nFDR (Adult)=",round(Ares[which(rownames(Ares)==disease[["als"]][i]),"padj"], digits=5),
                    "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)==disease[["als"]][i]),"padj"], digits=5))) +
     theme(title = element_text(size = 20)) +
     theme(text = element_text(size = 20))
   print(plots[[i]])
 }
 dev.off()
 
 ## lateral sclerosis
 pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/LateralSclerosis_gene_expression_byAge_byFraction_byLibrary.pdf", width = 7, height = 5)
 for (i in 1:length(disease[["ls"]])){
   print(i)
   plots[[i]] = ggplot(diseaseRPKM[["ls"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["ls"]][i]))) +
     facet_grid(. ~ Library) + scale_fill_brewer(palette="Dark2") +
     labs(fill="") +
     ylab("log(RPKM+1)") + 
     xlab("") +
     ggtitle(paste0(sym[["ls"]][i]," (", disease[["ls"]][i], ")\nFDR (Adult)=",round(Ares[which(rownames(Ares)==disease[["ls"]][i]),"padj"], digits=5),
                    "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)==disease[["ls"]][i]),"padj"], digits=5))) +
     theme(title = element_text(size = 20)) +
     theme(text = element_text(size = 20))
   print(plots[[i]])
 }
 dev.off()   
 
 pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/LateralSclerosis_gene_expression_byAge_byFraction.pdf", width = 7, height = 5)
 for (i in 1:length(disease[["ls"]])){
   print(i)
   plots[[i]] = ggplot(diseaseRPKM[["ls"]], aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=get(disease[["ls"]][i]))) +
     labs(fill="") + scale_fill_brewer(palette="Dark2") +
     ylab("log(RPKM+1)") +
     xlab("") +
     ggtitle(paste0(sym[["ls"]][i]," (", disease[["ls"]][i], ")\nFDR (Adult)=",round(Ares[which(rownames(Ares)==disease[["ls"]][i]),"padj"], digits=5),
                    "\nFDR (Prenatal)=",round(Fres.down[which(rownames(Fres.down)==disease[["ls"]][i]),"padj"], digits=5))) +
     theme(title = element_text(size = 20)) +
     theme(text = element_text(size = 20))
   print(plots[[i]])
 }
 dev.off()
