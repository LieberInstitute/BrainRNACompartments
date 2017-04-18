library(VennDiagram)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

### Which Fractional DE Genes Overlap? ###
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = unlist(DirList, recursive = F)
names(DirList) = c("Adult\nPolyA\nCytosolic", "Adult\nPolyA\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear", 
                   "Adult\nRibozero\nCytosolic", "Adult\nRibozero\nNuclear", "Fetal\nRibozero\nCytosolic", 
                   "Fetal\nRibozero\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear")
elementNROWS(SigFracList)
Names.Zone <- list("Adult:PolyA\n(1894)"=SigFracList[["Apres"]], "Fetal:PolyA\n(52)"=SigFracList[["Fpres"]], 
                   "Adult:Ribozero\n(1892)"=SigFracList[["Arres"]], "Fetal:Ribozero\n(30)"=SigFracList[["Frres"]])
Names.Zone = lapply(Names.Zone, function(x) rownames(x))
Names.Zone.down <- list("Adult:PolyA\n(1894)"=SigFracList[["Apres"]], "Fetal:PolyA\n(40)"=SigFracList[["Fpres.down"]], 
                   "Adult:Ribozero\n(1892)"=SigFracList[["Arres"]], "Fetal:Ribozero\n(30)"=SigFracList[["Frres"]])
Names.Zone.down = lapply(Names.Zone.down, function(x) rownames(x))

elementNROWS(DirList)
DirList = lapply(DirList, function(x) rownames(x))
Names.ZoneUP <- list("Adult:PolyA\n(956)"=DirList[["Adult\nPolyA\nNuclear"]], "Fetal:PolyA\n(51)"=DirList[[4]], 
                     "Adult:Ribozero\n(868)"=DirList[["Adult\nRibozero\nNuclear"]], "Fetal:Ribozero\n(23)"=DirList[["Fetal\nRibozero\nNuclear"]])
Names.ZoneD <- list("Adult:PolyA\n(938)"=DirList[["Adult\nPolyA\nCytosolic"]], "Fetal:PolyA\n(1)"=DirList[[3]], 
                    "Adult:Ribozero\n(1024)"=DirList[["Adult\nRibozero\nCytosolic"]], "Fetal:Ribozero\n(7)"=DirList[["Fetal\nRibozero\nCytosolic"]])
Names.ZoneUP.down <- list("Adult:PolyA\n(956)"=DirList[["Adult\nPolyA\nNuclear"]], "Fetal:PolyA\n(39)"=DirList[[10]], 
                     "Adult:Ribozero\n(868)"=DirList[["Adult\nRibozero\nNuclear"]], "Fetal:Ribozero\n(23)"=DirList[["Fetal\nRibozero\nNuclear"]])
Names.ZoneD.down <- list("Adult:PolyA\n(938)"=DirList[["Adult\nPolyA\nCytosolic"]], "Fetal:PolyA\n(1)"=DirList[[9]], 
                    "Adult:Ribozero\n(1024)"=DirList[["Adult\nRibozero\nCytosolic"]], "Fetal:Ribozero\n(7)"=DirList[["Fetal\nRibozero\nCytosolic"]])

venn.ZoneUP <- venn.diagram(Names.ZoneUP, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.nuclear.jpeg", 
                            main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                            col = "transparent",
                            fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                            alpha = 0.50,
                            label.col = c("olivedrab4", "white", "darkorchid4", "white","white", "white", "white", "white",
                                          "palevioletred4","white","white", "white", "white", "darkblue", "white"),
                            fontfamily = "Arial",
                            fontface = "bold",
                            cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                            cat.fontfamily = "Arial", margin=0.2)
venn.ZoneD <- venn.diagram(Names.ZoneD, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.cytosolic.jpeg", 
                           main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                           col = "transparent",
                           fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                           alpha = 0.50,
                           label.col = c("olivedrab4", "white", "darkorchid4", "white","white", "white", "white", "white",
                                         "palevioletred4","white","white", "white", "white", "darkblue", "white"),
                           fontfamily = "Arial",
                           fontface = "bold",
                           cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                           cat.fontfamily = "Arial", margin=0.2)
venn.zone <- venn.diagram(Names.Zone, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.zone.jpeg", 
                          main="Differentially expressed Genes\nby Fraction (LFC≥1, FDR<0.05)",
                          col = "transparent",
                          fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                          alpha = 0.50,
                          label.col = c("olivedrab4", "white", "darkorchid4", "white","white", "white", "white", "white",
                                        "palevioletred4","white","white", "white", "white", "darkblue", "white"),
                          fontfamily = "Arial",
                          fontface = "bold",
                          cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                          cat.fontfamily = "Arial", margin=0.2)
venn.ZoneUP.down <- venn.diagram(Names.ZoneUP.down, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.nuclear.downsampled.jpeg", 
                            main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                            col = "transparent",
                            fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                            alpha = 0.50,
                            label.col = c("olivedrab4", "white", "darkorchid4", "white","white", "white", "white", "white",
                                          "palevioletred4","white","white", "white", "white", "darkblue", "white"),
                            fontfamily = "Arial",
                            fontface = "bold",
                            cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                            cat.fontfamily = "Arial", margin=0.2)
venn.ZoneD.down <- venn.diagram(Names.ZoneD.down, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.cytosolic.downsampled.jpeg", 
                           main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                           col = "transparent",
                           fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                           alpha = 0.50,
                           label.col = c("olivedrab4", "white", "darkorchid4", "white","white", "white", "white", "white",
                                         "palevioletred4","white","white", "white", "white", "darkblue", "white"),
                           fontfamily = "Arial",
                           fontface = "bold",
                           cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                           cat.fontfamily = "Arial", margin=0.2)
venn.zone.down <- venn.diagram(Names.Zone.down, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.zone.downsampled.jpeg", 
                          main="Differentially expressed Genes\nby Fraction (LFC≥1, FDR<0.05)",
                          col = "transparent",
                          fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                          alpha = 0.50,
                          label.col = c("olivedrab4", "white", "darkorchid4", "white","white", "white", "white", "white",
                                        "palevioletred4","white","white", "white", "white", "darkblue", "white"),
                          fontfamily = "Arial",
                          fontface = "bold",
                          cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                          cat.fontfamily = "Arial", margin=0.2)