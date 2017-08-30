library(clusterProfiler)
library(GenomicRanges)
library(ggplot2)


load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.downsampled.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23_nodownsamp.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Disease Ontology
plot(compareDO,colorBy="p.adjust",  showCategory = 45, title= "Disease Ontology Enrichment")
DO = as.data.frame(compareDO)
DO.LFC1 = read.csv("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.DO.downsampled.csv")

disease = list(az = as.character(DO[which(DO$Description=="Alzheimer's disease"),"geneID"]),
               als = as.character(DO[which(DO$Description=="amyotrophic lateral sclerosis"),"geneID"]),
               ls = as.character(DO[which(DO$Description=="lateral sclerosis"),"geneID"]),
               motor = as.character(DO[which(DO$Description=="motor neuron disease"),"geneID"]))

disease = lapply(disease, function(x) strsplit(as.character(x), "/", fixed=TRUE))
disease = lapply(disease, function(x) x[[1]])
elementNROWS(disease)

disease = lapply(disease, function(x) geneMap[which(geneMap$EntrezID %in% x),"gencodeID"])

# make rpkm object for ggplot
geneRpkm[1:5,1:5] 
x= t(geneRpkm)
x = data.frame(x)
match(rownames(x), pd$SampleID)
x$Age = pd$Fetal 
x$Fraction = pd$Zone
x$Library = pd$Library

diseaseRPKM = lapply(disease, function(y) x[,which(colnames(x) %in% y | 
                                     colnames(x)=="Age" | colnames(x)== "Fraction"| colnames(x)== "Library")])

sym = lapply(disease, function(a) geneMap[which(geneMap$gencodeID %in% a),"Symbol"])
for (i in 1:length(disease)){
  tmp = diseaseRPKM[[i]]
  tmp[,1:length(disease[[i]])] = log(tmp[,1:length(disease[[i]])]+1)
  diseaseRPKM[[i]] = tmp
}

### plot expression of disease-associated genes
# Prader-Willi
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/NPAP1_prader-willi_byFraction_byAge_byLibrary.pdf", width = 7, height = 5)
ggplot(x, aes(fill = Fraction)) + geom_boxplot(aes(x=Age,y=ENSG00000185823.3_1)) +
  facet_grid(. ~ Library) +
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
  labs(fill="") +
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
  facet_grid(. ~ Library) +
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
  labs(fill="") +
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
    facet_grid(. ~ Library) +
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
   labs(fill="") +
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
     labs(fill="") +
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
     labs(fill="") +
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
     labs(fill="") +
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
     labs(fill="") +
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
     facet_grid(. ~ Library) +
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
