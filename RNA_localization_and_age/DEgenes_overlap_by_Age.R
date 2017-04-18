library("VennDiagram")

### Which Age DEGs overlap? ###

# Make list of IDs of DEGs
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres), 
               Crres = data.frame(Crres), Nrres = data.frame(Nrres),
               Cpres.down = data.frame(Cpres.down))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "DownFetal"))
Sign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = unlist(DirList, recursive = F) 
names(DirList) = c("Cytosol\nPolyA\nIncreasing", "Cytosol\nPolyA\nDecreasing", "Nucleus\nPolyA\nIncreasing", "Nucleus\nPolyA\nDecreasing", 
                   "Cytosol\nRibozero\nIncreasing", "Cytosol\nRibozero\nDecreasing", "Nucleus\nRibozero\nIncreasing", 
                   "Nucleus\nRibozero\nDecreasing", "Cytosol\nPolyA\nIncreasing", "Cytosol\nPolyA\nDecreasing")
elementNROWS(SigAgeList)
Names.Age <- list("Cytosol:PolyA\n(9021)"=SigAgeList[["Cpres"]], "Nucleus:PolyA\n(7660)"=SigAgeList[["Npres"]], 
                  "Cytosol:Ribozero\n(7163)"=SigAgeList[["Crres"]], "Nucleus:Ribozero\n(5436)"=SigAgeList[["Nrres"]])
Names.Age = lapply(Names.Age, function(x) rownames(x))
Names.Age.down <- list("Cytosol:PolyA\n(8356)"=SigAgeList[["Cpres.down"]], "Nucleus:PolyA\n(7870)"=SigAgeList[["Npres"]], 
                  "Cytosol:Ribozero\n(7163)"=SigAgeList[["Crres"]], "Nucleus:Ribozero\n(5436)"=SigAgeList[["Nrres"]])
Names.Age.down = lapply(Names.Age.down, function(x) rownames(x))

elementNROWS(DirList)
DirList = lapply(DirList, function(x) rownames(x))
Names.AgeUP <- list("Cytosol:PolyA\n(4170)"=DirList[[2]], "Nucleus:PolyA\n(3417)"=DirList[["Nucleus\nPolyA\nDecreasing"]], 
                    "Cytosol:Ribozero\n(3043)"=DirList[["Cytosol\nRibozero\nDecreasing"]], "Nucleus:Ribozero\n(2600)"=DirList[["Nucleus\nRibozero\nDecreasing"]])
Names.AgeD <- list("Cytosol:PolyA\n(4851)"=DirList[[1]], "Nucleus:PolyA\n(4243)"=DirList[["Nucleus\nPolyA\nIncreasing"]], 
                   "Cytosol:Ribozero\n(4120)"=DirList[["Cytosol\nRibozero\nIncreasing"]], "Nucleus:Ribozero\n(2836)"=DirList[["Nucleus\nRibozero\nIncreasing"]])
Names.AgeUP.down <- list("Cytosol:PolyA\n(3808)"=DirList[[10]], "Nucleus:PolyA\n(3417)"=DirList[["Nucleus\nPolyA\nDecreasing"]], 
                    "Cytosol:Ribozero\n(3043)"=DirList[["Cytosol\nRibozero\nDecreasing"]], "Nucleus:Ribozero\n(2600)"=DirList[["Nucleus\nRibozero\nDecreasing"]])
Names.AgeD.down <- list("Cytosol:PolyA\n(4548)"=DirList[[9]], "Nucleus:PolyA\n(4243)"=DirList[["Nucleus\nPolyA\nIncreasing"]], 
                   "Cytosol:Ribozero\n(4120)"=DirList[["Cytosol\nRibozero\nIncreasing"]], "Nucleus:Ribozero\n(2836)"=DirList[["Nucleus\nRibozero\nIncreasing"]])

### Which Age DEGs overlap? ###
venn.age <- venn.diagram(Names.Age, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.age.jpeg", 
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
venn.ageUP <- venn.diagram(Names.AgeUP, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.ageDecreasing.jpeg", 
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
venn.ageD <- venn.diagram(Names.AgeD, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.ageIncreasing.jpeg",
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
venn.age.down <- venn.diagram(Names.Age.down, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.age.downsampled.jpeg", 
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
venn.ageUP.down <- venn.diagram(Names.AgeUP.down, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.ageDecreasing.downsampled.jpeg", 
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
venn.ageD.down <- venn.diagram(Names.AgeD.down, "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.ageIncreasing.downsampled.jpeg",
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