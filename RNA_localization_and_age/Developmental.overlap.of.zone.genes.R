library(VennDiagram)

### Which Fractional DE Genes Overlap? ###
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))
develop = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/developmentally-regulated.genes.csv")

# At LFC >= 1
devel = develop[which(develop$padj<=0.05 & abs(develop$log2FoldChange)>=1),]
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigFracList)

Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]],
               Devel = devel) 
names(DirList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", 
                   "Fetal\nRibozero\nDown", "Developmentally\nRegulated")

elementLengths(DirList)
DirList = lapply(DirList, function(x) as.character(x$X))
Names.Zone.P.Devel.UP <- list("Adult:PolyA\n(954)"=DirList[["Adult\nPolyA\nUp"]], 
                              "Fetal:PolyA\n(74)"=DirList[["Fetal\nPolyA\nUp"]], 
                              "Developmentally\nRegulated (12,618)"=DirList[["Developmentally\nRegulated"]])
Names.Zone.P.Devel.D <- list("Adult:PolyA\n(694)"=DirList[["Adult\nPolyA\nDown"]], 
                             "Fetal:PolyA\n(1)"=DirList[["Fetal\nPolyA\nDown"]],
                             "Developmentally\nRegulated (12,618)"=DirList[["Developmentally\nRegulated"]])

venn.Zone.P.Devel.UP <- venn.diagram(Names.Zone.P.Devel.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.P.Devel.UP.jpeg", 
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.P.Devel.D <- venn.diagram(Names.Zone.P.Devel.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.P.Devel.D.jpeg", 
                                    main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)

Names.Zone.R.Devel.UP <- list("Adult:Ribozero\n(1,507)"=DirList[["Adult\nRibozero\nUp"]], 
                              "Fetal:Ribozero\n(339)"=DirList[["Fetal\nRibozero\nUp"]], 
                              "Developmentally\nRegulated (12,618)"=DirList[["Developmentally\nRegulated"]])
Names.Zone.R.Devel.D <- list("Adult:Ribozero\n(1,362)"=DirList[["Adult\nRibozero\nDown"]], 
                             "Fetal:Ribozero\n(10)"=DirList[["Fetal\nRibozero\nDown"]],
                             "Developmentally\nRegulated (12,618)"=DirList[["Developmentally\nRegulated"]])

venn.Zone.R.Devel.UP <- venn.diagram(Names.Zone.R.Devel.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.R.Devel.UP.jpeg", 
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.R.Devel.D <- venn.diagram(Names.Zone.R.Devel.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.R.Devel.D.jpeg", 
                                    main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)

# At lfc >=2

devel = develop[which(develop$padj<=0.05 & abs(develop$log2FoldChange)>=2),]
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=2),])
elementLengths(SigFracList)

Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]],
               Devel = devel) 
names(DirList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", 
                   "Fetal\nRibozero\nDown", "Developmentally\nRegulated")

elementLengths(DirList)
DirList = lapply(DirList, function(x) as.character(x$X))

Names.Zone.P.Devel.UP <- list("Adult:PolyA\n(61)"=DirList[["Adult\nPolyA\nUp"]], 
                              "Fetal:PolyA\n(1)"=DirList[["Fetal\nPolyA\nUp"]], 
                              "Developmentally\nRegulated (6,800)"=DirList[["Developmentally\nRegulated"]])
Names.Zone.P.Devel.D <- list("Adult:PolyA\n(1)"=DirList[["Adult\nPolyA\nDown"]], 
                             "Fetal:PolyA\n(0)"=DirList[["Fetal\nPolyA\nDown"]],
                             "Developmentally\nRegulated (6,800)"=DirList[["Developmentally\nRegulated"]])

venn.Zone.P.Devel.UP <- venn.diagram(Names.Zone.P.Devel.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.P.Devel.UP.lfc2.jpeg", 
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥2, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.P.Devel.D <- venn.diagram(Names.Zone.P.Devel.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.P.Devel.D.lfc2.jpeg", 
                                    main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥2, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)

Names.Zone.R.Devel.UP <- list("Adult:Ribozero\n(174)"=DirList[["Adult\nRibozero\nUp"]], 
                              "Fetal:Ribozero\n(19)"=DirList[["Fetal\nRibozero\nUp"]], 
                              "Developmentally\nRegulated (6,800)"=DirList[["Developmentally\nRegulated"]])
Names.Zone.R.Devel.D <- list("Adult:Ribozero\n(32)"=DirList[["Adult\nRibozero\nDown"]], 
                             "Fetal:Ribozero\n(0)"=DirList[["Fetal\nRibozero\nDown"]],
                             "Developmentally\nRegulated (6,800)"=DirList[["Developmentally\nRegulated"]])

venn.Zone.R.Devel.UP <- venn.diagram(Names.Zone.R.Devel.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.R.Devel.UP.lfc2.jpeg", 
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥2, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.R.Devel.D <- venn.diagram(Names.Zone.R.Devel.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.zone.R.Devel.D.lfc2.jpeg", 
                                    main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥2, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)