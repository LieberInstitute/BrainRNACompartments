library(VennDiagram)

### Which Fractional DE Genes Overlap? ###
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]]) 
names(DirList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", "Fetal\nRibozero\nDown")
elementLengths(SigFracList)
Names.Zone <- list("Adult:PolyA\n(1648)"=SigFracList[["Apres"]], "Fetal:PolyA\n(75)"=SigFracList[["Fpres"]], 
                   "Adult:Ribozero\n(2869)"=SigFracList[["Arres"]], "Fetal:Ribozero\n(349)"=SigFracList[["Frres"]])
Names.Zone = lapply(Names.Zone, function(x) x$X)

elementLengths(DirList)
DirList = lapply(DirList, function(x) as.character(x$X))
Names.ZoneUP <- list("Adult:PolyA\n(954)"=DirList[["Adult\nPolyA\nUp"]], "Fetal:PolyA\n(74)"=DirList[["Fetal\nPolyA\nUp"]], 
                     "Adult:Ribozero\n(1507)"=DirList[["Adult\nRibozero\nUp"]], "Fetal:Ribozero\n(339)"=DirList[["Fetal\nRibozero\nUp"]])
Names.ZoneD <- list("Adult:PolyA\n(694)"=DirList[["Adult\nPolyA\nDown"]], "Fetal:PolyA\n(1)"=DirList[["Fetal\nPolyA\nDown"]], 
                    "Adult:Ribozero\n(1362)"=DirList[["Adult\nRibozero\nDown"]], "Fetal:Ribozero\n(10)"=DirList[["Fetal\nRibozero\nDown"]])

list.Zone <- calculate.overlap(Names.Zone)
list.ZoneUP <- calculate.overlap(Names.ZoneUP)
list.ZoneD <- calculate.overlap(Names.ZoneD)
venn.ZoneUP <- venn.diagram(Names.ZoneUP, "/Users/amanda/Dropbox/NucVsCytosol/Results/venn.zoneUP.jpeg", main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                            col = "transparent",
                            fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                            alpha = 0.50,
                            label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                          "white", "white", "white", "white", "palevioletred4", "white",
                                          "white", "white", "white", "darkblue", "white"),
                            fontfamily = "Arial",
                            fontface = "bold",
                            cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                            cat.fontfamily = "Arial", margin=0.2)
venn.ZoneD <- venn.diagram(Names.ZoneD, "/Users/amanda/Dropbox/NucVsCytosol/Results/venn.zoneD.jpeg", main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                           col = "transparent",
                           fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                           alpha = 0.50,
                           label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                         "white", "white", "white", "white", "palevioletred4", "white",
                                         "white", "white", "white", "darkblue", "white"),
                           fontfamily = "Arial",
                           fontface = "bold",
                           cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                           cat.fontfamily = "Arial", margin=0.2)
venn.zone <- venn.diagram(Names.Zone, "/Users/amanda/Dropbox/NucVsCytosol/Results/venn.zone.jpeg", main="Differentially expressed Genes\nby Fraction (LFC≥1, FDR<0.05)",
                          col = "transparent",
                          fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                          alpha = 0.50,
                          label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                        "white", "white", "white", "white", "palevioletred4", "white",
                                        "white", "white", "white", "darkblue", "white"),
                          fontfamily = "Arial",
                          fontface = "bold",
                          cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                          cat.fontfamily = "Arial", margin=0.2)

### Characterizing the significantly different genes that don't overlap ###



