library("VennDiagram")

### Which Age DEGs overlap? ###

# Make list of Ensembl IDs of DEGs
AgeList = list(Cpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Cpres.csv"),
                Npres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Npres.csv"),
                Crres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Crres.csv"),
                Nrres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nrres.csv"))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "DownFetal"))
Sign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Cpres.Up = DirList[["Cpres"]][["UpFetal"]], Cpres.Down = DirList[["Cpres"]][["DownFetal"]],
               Npres.Up = DirList[["Npres"]][["UpFetal"]], Npres.Down = DirList[["Npres"]][["DownFetal"]],
               Crres.Up = DirList[["Crres"]][["UpFetal"]], Crres.Down = DirList[["Crres"]][["DownFetal"]],
               Nrres.Up = DirList[["Nrres"]][["UpFetal"]], Nrres.Down = DirList[["Nrres"]][["DownFetal"]]) 
names(DirList) = c("Cytosol\nPolyA\nUp", "Cytosol\nPolyA\nDown", "Nucleus\nPolyA\nUp", "Nucleus\nPolyA\nDown", 
                   "Cytosol\nRibozero\nUp", "Cytosol\nRibozero\nDown", "Nucleus\nRibozero\nUp", "Nucleus\nRibozero\nDown")
elementLengths(SigAgeList)
Names.Age <- list("Cytosol:PolyA\n(9273)"=SigAgeList[["Cpres"]], "Nucleus:PolyA\n(7870)"=SigAgeList[["Npres"]], 
                  "Cytosol:Ribozero\n(8270)"=SigAgeList[["Crres"]], "Nucleus:Ribozero\n(6221)"=SigAgeList[["Nrres"]])
Names.Age = lapply(Names.Age, function(x) x$X)

elementLengths(DirList)
DirList = lapply(DirList, function(x) as.character(x$X))
Names.AgeUP <- list("Cytosol:PolyA\n(4238)"=DirList[["Cytosol\nPolyA\nUp"]], "Nucleus:PolyA\n(3427)"=DirList[["Nucleus\nPolyA\nUp"]], 
                    "Cytosol:Ribozero\n(3182)"=DirList[["Cytosol\nRibozero\nUp"]], "Nucleus:Ribozero\n(2942)"=DirList[["Nucleus\nRibozero\nUp"]])
Names.AgeUP = lapply(Names.AgeUP, function(x) x$X)
Names.AgeD <- list("Cytosol:PolyA\n(5035)"=DirList[["Cytosol\nPolyA\nDown"]], "Nucleus:PolyA\n(4443)"=DirList[["Nucleus\nPolyA\nDown"]], 
                   "Cytosol:Ribozero\n(5088)"=DirList[["Cytosol\nRibozero\nDown"]], "Nucleus:Ribozero\n(3279)"=DirList[["Nucleus\nRibozero\nDown"]])
Names.AgeD = lapply(Names.AgeD, function(x) x$X)

### Which Age DEGs overlap? ###
list.age <- calculate.overlap(Names.Age)
list.ageUP <- calculate.overlap(Names.AgeUP)
list.ageD <- calculate.overlap(Names.AgeD)
venn.age <- venn.diagram(Names.Age, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.age.jpeg", 
                         main="Differentially expressed Genes\nOver Brain Development\n(LFC≥1, FDR<0.05)",
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
venn.ageUP <- venn.diagram(Names.AgeUP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.ageUP.jpeg", 
                           main="Genes With Significantly Decreasing Expression\nOver Brain Development\n(LFC≥1, FDR<0.05)",
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
venn.ageD <- venn.diagram(Names.AgeD, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.ageD.jpeg",
                          main="Genes With Significantly Increasing Expression\nOver Brain Development\n(LFC≥1, FDR<0.05)",
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
