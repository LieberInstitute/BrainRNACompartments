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
devel = do.call(c, develop[1:2])
devel.down = do.call(c,develop[2:3])
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = unlist(DirList, recursive = F)
DirList = lapply(DirList, function(x) as.character(rownames(x)))
DirList = c(DirList,list(devel, devel.down))
names(DirList) = c("Adult\nPolyA\nCytosolic", "Adult\nPolyA\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear", 
                   "Adult\nRibozero\nCytosolic", "Adult\nRibozero\nNuclear", "Fetal\nRibozero\nCytosolic", 
                   "Fetal\nRibozero\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear", 
                   "Developmentally\nRegulated", "Developmentally\nRegulated")
DirList = lapply(DirList, unique)
elementNROWS(DirList)
Names.Zone.P.Devel.D <- list("Adult:PolyA\n(938)"=DirList[[1]], 
                              "Fetal:PolyA\n(1)"=DirList[[3]], 
                              "Developmentally\nRegulated (12991)"=DirList[[11]])
Names.Zone.P.Devel.UP <- list("Adult:PolyA\n(694)"=DirList[[2]], 
                             "Fetal:PolyA\n(51)"=DirList[[4]],
                             "Developmentally\nRegulated (12991)"=DirList[[11]])
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
Names.Zone.P.Devel.D.down <- list("Adult:PolyA\n(938)"=DirList[[1]], 
                                  "Fetal:PolyA\n(1)"=DirList[[9]], 
                                  "Developmentally\nRegulated (12667)"=DirList[[12]])
Names.Zone.P.Devel.UP.down <- list("Adult:PolyA\n(694)"=DirList[[2]], 
                                   "Fetal:PolyA\n(39)"=DirList[[10]],
                                   "Developmentally\nRegulated (12667)"=DirList[[12]])
venn.Zone.P.Devel.UP.down <- venn.diagram(Names.Zone.P.Devel.UP.down, 
                                     "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.DEGbyFraction&Age.polya.nuclear.downsampled.jpeg", 
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.P.Devel.D.down <- venn.diagram(Names.Zone.P.Devel.D.down, 
                                    "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.DEGbyFraction&Age.polya.cytosolic.downsampled.jpeg",
                                    main="Differentially expressed Genes\nUp-regulated in Cytosol (LFC≥1, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)
Names.Zone.R.Devel.UP <- list("Adult:Ribozero\n(868)"=DirList[["Adult\nRibozero\nNuclear"]], 
                              "Fetal:Ribozero\n(23)"=DirList[["Fetal\nRibozero\nNuclear"]], 
                              "Developmentally\nRegulated (12991)"=DirList[[11]])
Names.Zone.R.Devel.D <- list("Adult:Ribozero\n(1024)"=DirList[["Adult\nRibozero\nCytosolic"]], 
                             "Fetal:Ribozero\n(7)"=DirList[["Fetal\nRibozero\nCytosolic"]],
                             "Developmentally\nRegulated (12991)"=DirList[[11]])
venn.Zone.R.Devel.UP <- venn.diagram(Names.Zone.R.Devel.UP, 
                                     "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.DEGbyFraction&Age.ribo.nuclear.jpeg",
                                     main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                     col = "transparent",
                                     fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                     alpha = 0.50,
                                     fontfamily = "Arial",
                                     fontface = "bold",
                                     cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                     cat.fontfamily = "Arial", margin=0.2)
venn.Zone.R.Devel.D <- venn.diagram(Names.Zone.R.Devel.D, 
                                    "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.DEGbyFraction&Age.ribo.cytosolic.jpeg", 
                                    main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                    col = "transparent",
                                    fill = c("lightpink2","cornflowerblue", "olivedrab2"),
                                    alpha = 0.50,
                                    fontfamily = "Arial",
                                    fontface = "bold",
                                    cat.col = c("palevioletred4", "darkblue", "olivedrab4"),
                                    cat.fontfamily = "Arial", margin=0.2)