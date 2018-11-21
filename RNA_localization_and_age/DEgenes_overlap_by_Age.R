library("VennDiagram")

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

### Which Age DEGs overlap? ###

# Make list of IDs of DEGs
AgeList = list(Cpres = data.frame(Cpres.down), Npres = data.frame(Npres), 
               Crres = data.frame(Crres), Nrres = data.frame(Nrres))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "DownFetal"))
Sign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = unlist(DirList, recursive = F) 
names(DirList) = c("Cytosol\nPolyA\nIncreasing", "Cytosol\nPolyA\nDecreasing", "Nucleus\nPolyA\nIncreasing", "Nucleus\nPolyA\nDecreasing", 
                   "Cytosol\nRibozero\nIncreasing", "Cytosol\nRibozero\nDecreasing", "Nucleus\nRibozero\nIncreasing", "Nucleus\nRibozero\nDecreasing")
elementNROWS(SigAgeList)
names = c("Cytosol:PolyA\n(","Nucleus:PolyA\n(","Cytosol:Ribozero\n(","Nucleus:Ribozero\n(")
names(SigAgeList) <- paste0(names, elementNROWS(SigAgeList), ")")
Names.Age.down = lapply(SigAgeList, function(x) rownames(x))

elementNROWS(DirList)
DirList = lapply(DirList, function(x) rownames(x))
dir = list(up = grep("Decreasing", names(DirList)), down = grep("Increasing", names(DirList)))
Names.AgeUP.down = DirList[dir[["up"]]]
names(Names.AgeUP.down) = gsub("Decreasing", "", names(Names.AgeUP.down))
names(Names.AgeUP.down) = paste0(names(Names.AgeUP.down), "(", elementNROWS(Names.AgeUP.down), ")")
Names.AgeD.down = DirList[dir[["down"]]]
names(Names.AgeD.down) = gsub("Increasing", "", names(Names.AgeD.down))
names(Names.AgeD.down) = paste0(names(Names.AgeD.down), "(", elementNROWS(Names.AgeD.down), ")")
names(Names.Age.down) = gsub("Cytosol", "Cytoplasm", names(Names.Age.down))


### Which Age DEGs overlap? ###
## Only using downsampled bams
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


## Redo in black and white

names(Names.Age.down) = c("Cytoplasm\nPolyA\n(8356)","Nucleus\nPolyA\n(7660)",
                          "Cytoplasm\nRibozero\n(7163)","Nucleus\nRibozero\n(5436)")
venn.diagram(Names.Age.down, "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/venn.age.downsampled.jpeg", 
             main="Differentially expressed Genes\nOver Brain Development\n(LFC≥1, FDR<0.05)",
             fill = c("#1b9e77","#d95f02","#1b9e77", "#d95f02"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.3,
             cat.fontfamily = "Arial", margin=0.2)


# How many of the total of all age genes in all four groups not broken down by direction overlap?

round(3435/elementNROWS(SigAgeList)*100,1)
#Cytosol:PolyA\n(8356)    Nucleus:PolyA\n(7660) Cytosol:Ribozero\n(7163) Nucleus:Ribozero\n(5436) 
#41.1                     44.8                     48.0                     63.2