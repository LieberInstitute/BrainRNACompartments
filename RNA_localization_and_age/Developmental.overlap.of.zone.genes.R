library(VennDiagram)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

### Which Fractional DE Genes Overlap? ###
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down))
develop = list(Agepres = data.frame(Agepres), Agerres = data.frame(Agerres), 
               Agepres.down = data.frame(Agepres.down))
develop = Map(cbind, develop, GeneID = lapply(develop, rownames))

# At LFC >= 1
develop = lapply(develop, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange)>=1),])
develop = lapply(develop, function(x) as.character(x$GeneID))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = Map(cbind, SigFracList, Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc")))
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = unlist(DirList, recursive = F)
DirList = lapply(DirList, function(x) as.character(rownames(x)))
DirList = c(DirList,develop)
names(DirList) = c("Adult\nPolyA\nCytosolic", "Adult\nPolyA\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear", 
                   "Adult\nRibozero\nCytosolic", "Adult\nRibozero\nNuclear", "Fetal\nRibozero\nCytosolic", 
                   "Fetal\nRibozero\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear", 
                   "Developmentally\nRegulated (polyA)", "Developmentally\nRegulated (Ribozero)", "Developmentally Regulated\n(polyA Downsampled)")
DirList = lapply(DirList, unique)
elementNROWS(DirList)
Names.Zone.P.Devel.D <- list("Adult:PolyA\n(938)"=DirList[["Adult\nPolyA\nCytosolic"]], 
                              "Fetal:PolyA\n(1)"=DirList[["Fetal\nPolyA\nCytosolic"]], 
                              "Developmentally\nRegulated (10097)"=DirList[["Developmentally Regulated\n(polyA Downsampled)"]])
Names.Zone.P.Devel.UP <- list("Adult:PolyA\n(956)"=DirList[["Adult\nPolyA\nNuclear"]], 
                             "Fetal:PolyA\n(51)"=DirList[["Fetal\nPolyA\nNuclear"]],
                             "Developmentally\nRegulated (10097)"=DirList[["Developmentally Regulated\n(polyA Downsampled)"]])
venn.Zone.P.Devel.UP <- venn.diagram(Names.Zone.P.Devel.UP, 
                                     "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.DEGbyFraction&Age.polya.nuclear.jpeg", 
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.P.Devel.D <- venn.diagram(Names.Zone.P.Devel.D, 
                                    "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.DEGbyFraction&Age.polya.cytosolic.jpeg",
                                    main="Differentially expressed Genes\nUp-regulated in Cytosol (LFC≥1, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)


## Are genes significantly regulated by fraction more or less likely to also be regulated by age?

genes = data.frame(geneID = rownames(Agepres.down), A.LFC = data.frame(Apres)[match(rownames(Agepres.down), rownames(Apres)),"log2FoldChange"],
                   A.padj = data.frame(Apres)[match(rownames(Agepres.down), rownames(Apres)),"padj"],
                   P.LFC = data.frame(Fpres.down)[match(rownames(Agepres.down), rownames(Fpres.down)),"log2FoldChange"],
                   P.padj = data.frame(Fpres.down)[match(rownames(Agepres.down), rownames(Fpres.down)),"padj"],
                   C.LFC = data.frame(Cpres.down)[match(rownames(Agepres.down), rownames(Cpres.down)),"log2FoldChange"],
                   C.padj = data.frame(Cpres.down)[match(rownames(Agepres.down), rownames(Cpres.down)),"padj"],
                   N.LFC = data.frame(Npres)[match(rownames(Agepres.down), rownames(Npres)),"log2FoldChange"],
                   N.padj = data.frame(Npres)[match(rownames(Agepres.down), rownames(Npres)),"padj"])

devel = list(CytosolIncr = genes[which(genes$C.LFC<0 & genes$C.padj<=0.05),], CytosolDecr = genes[which(genes$C.LFC>0 & genes$C.padj<=0.05),],
             NucleusIncr = genes[which(genes$N.LFC<0 & genes$N.padj<=0.05),], NucleusDecr = genes[which(genes$N.LFC>0 & genes$N.padj<=0.05),])
tb = list()
for (i in 1:length(devel)) {
  tmp = devel[[i]]
  tb[[i]] = list(AdultFrac = data.frame(Sig = c(length(unique(tmp[which(tmp$A.padj<=0.05 & tmp$A.LFC<0),"geneID"])),
                                                length(unique(tmp[which(tmp$A.padj<=0.05 & tmp$A.LFC>0),"geneID"]))),
                                        nonSig = c(length(unique(tmp[which(tmp$A.padj>0.05 & tmp$A.LFC<0),"geneID"])),
                                                   length(unique(tmp[which(tmp$A.padj>0.05 & tmp$A.LFC>0),"geneID"]))),
                                        row.names = c("Cytosolic","Nuclear")),
                 PrenatalFrac = data.frame(Sig = c(length(unique(tmp[which(tmp$P.padj<=0.05 & tmp$P.LFC<0),"geneID"])),
                                                   length(unique(tmp[which(tmp$P.padj<=0.05 & tmp$P.LFC>0),"geneID"]))),
                                           nonSig = c(length(unique(tmp[which(tmp$P.padj>0.05 & tmp$P.LFC<0),"geneID"])),
                                                      length(unique(tmp[which(tmp$P.padj>0.05 & tmp$P.LFC>0),"geneID"]))),
                                           row.names = c("Cytosolic","Nuclear")))
}
names(tb) = names(devel)
elementNROWS(tb)
df = rbind(pvalue = unlist(lapply(lapply(tb, function(x) lapply(x, fisher.test)), function(x) lapply(x, function(y) y$p.value))),
           OR = unlist(lapply(lapply(tb, function(x) lapply(x, fisher.test)), function(x) lapply(x, function(y) y$estimate)))) 
Counts = do.call(rbind, unlist(tb, recursive=F))
Counts$DevelopmentGeneSet = rownames(Counts)
Counts$DevelopmentGeneSet = gsub(".Cytosolic", "",Counts$DevelopmentGeneSet)
Counts$DevelopmentGeneSet = gsub(".Nuclear", "",Counts$DevelopmentGeneSet)
Counts$pval = df["pvalue",match(Counts$DevelopmentGeneSet,colnames(df))]
Counts$OR = df["OR",match(Counts$DevelopmentGeneSet,colnames(df))]
write.csv(Counts,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/FracLFCxFracpval_bysigAgeGroup_fisher.csv")

Counts[which(Counts$pval*8<=0.05),]
#                                    Sig nonSig       DevelopmentGeneSet         pval        OR
#CytosolIncr.AdultFrac.Cytosolic    1290   2260    CytosolIncr.AdultFrac 3.478374e-17 1.6741685
#CytosolIncr.AdultFrac.Nuclear       510   1496    CytosolIncr.AdultFrac 3.478374e-17 1.6741685
#CytosolIncr.PrenatalFrac.Cytosolic   27   1564 CytosolIncr.PrenatalFrac 1.904316e-09 0.2832729
#CytosolIncr.PrenatalFrac.Nuclear     82   1345 CytosolIncr.PrenatalFrac 1.904316e-09 0.2832729
#CytosolDecr.PrenatalFrac.Cytosolic   30   2604 CytosolDecr.PrenatalFrac 8.184570e-08 0.3158532
#CytosolDecr.PrenatalFrac.Nuclear     66   1809 CytosolDecr.PrenatalFrac 8.184570e-08 0.3158532
#NucleusIncr.PrenatalFrac.Cytosolic   21   1108 NucleusIncr.PrenatalFrac 2.355003e-10 0.2526443
#NucleusIncr.PrenatalFrac.Nuclear    105   1399 NucleusIncr.PrenatalFrac 2.355003e-10 0.2526443
#NucleusDecr.AdultFrac.Cytosolic     775   1643    NucleusDecr.AdultFrac 6.234787e-30 2.2550791
#NucleusDecr.AdultFrac.Nuclear       348   1664    NucleusDecr.AdultFrac 6.234787e-30 2.2550791
#NucleusDecr.PrenatalFrac.Cytosolic   31   2700 NucleusDecr.PrenatalFrac 2.579350e-06 0.3461096
#NucleusDecr.PrenatalFrac.Nuclear     51   1537 NucleusDecr.PrenatalFrac 2.579350e-06 0.3461096

