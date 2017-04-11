load("./Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/rawCounts_nucleusVsCytosol_n24.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_Amanda_polya_n12.rda")
exonCountsP = exonCounts
geneCountsP = geneCounts
jCountsP = jCounts
jMapP = jMap
metricsP = metrics
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_Amanda_ribozero_n12.rda")
exonCountsR = exonCounts
geneCountsR = geneCounts
jCountsR = jCounts
jMapR = jMap
metricsR = metrics
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_Amanda_downsamp_n2.rda")
exonCountsD = exonCounts
geneCountsD = geneCounts
jCountsD = jCounts
jMapD = jMap
metricsD = metrics

# Combine files
metrics = rbind(metricsP,metricsR,metricsD)
geneCounts = cbind(geneCountsP, geneCountsR)
exonCounts = cbind(exonCountsP, exonCountsR)
geneCounts.down = cbind(geneCounts[,c(1:6,8,10:24)],geneCountsD)
exonCounts.down = cbind(exonCounts[,c(1:6,8,10:24)],exonCountsD)
colnames(exonCounts.down) = colnames(geneCounts.down) = gsub("downsamp", "polyA", colnames(geneCounts.down))
jMap = c(jMapP,jMapR,jMapD)
jMap = sort(unique(jMap))

rm(exonCountsP, exonCountsR,exonCountsD,geneCountsP, geneCountsR,geneCountsD,metricsP,metricsR,metricsD,
   jCounts,jMapP,jMapR,jMapD)

# filter 
geneCounts = geneCounts[which(rowSums(geneCounts) > 0),which(colnames(geneCounts)!="Br1113N1_RiboZero")]
exonCounts = exonCounts[,which(colnames(exonCounts)!="Br1113N1_RiboZero")]
pd = pd[which(pd$SampleID!="Br1113N1_RiboZero"),]
geneCounts.down = geneCounts.down[which(rowSums(geneCounts.down) > 0),which(colnames(geneCounts.down)!="Br1113N1_RiboZero")]
exonCounts.down = exonCounts.down[,which(colnames(exonCounts.down)!="Br1113N1_RiboZero")]
pd$Label = factor(paste(pd$Fetal, pd$Zone,pd$Library, sep="\n"), 
                  levels = c("Adult\nCytosol\npolyA", "Fetal\nCytosol\npolyA", "Adult\nNucleus\npolyA",
                             "Fetal\nNucleus\npolyA", "Adult\nCytosol\nRiboZero", "Fetal\nCytosol\nRiboZero",
                             "Adult\nNucleus\nRiboZero", "Fetal\nNucleus\nRiboZero"))
pd$WorkingID = c("Adult1_Cytosol_polyA", "Adult1_Nucleus_polyA", "Adult2_Cytosol_polyA", "Adult2_Nucleus_polyA",
                "Adult3_Cytosol_polyA", "Adult3_Nucleus_polyA", "Fetal1_Cytosol_polyA", "Fetal1_Nucleus_polyA",
                "Fetal2_Cytosol_polyA", "Fetal2_Nucleus_polyA", "Fetal3_Cytosol_polyA", "Fetal3_Nucleus_polyA",
                "Adult1_Cytosol_RiboZero", "Adult2_Cytosol_RiboZero", "Adult2_Nucleus_RiboZero",
                "Adult3_Cytosol_RiboZero", "Adult3_Nucleus_RiboZero", "Fetal1_Cytosol_RiboZero", "Fetal1_Nucleus_RiboZero",
                "Fetal2_Cytosol_RiboZero", "Fetal2_Nucleus_RiboZero", "Fetal3_Cytosol_RiboZero", "Fetal3_Nucleus_RiboZero")
pd = pd[order(pd$SampleID),]
rownames(pd) = pd$SampleID
geneCounts = geneCounts[,order(colnames(geneCounts))]
geneCounts.down = geneCounts.down[,order(colnames(geneCounts.down))]
exonCounts = exonCounts[,order(colnames(exonCounts))]
exonCounts.down = exonCounts.down[,order(colnames(exonCounts.down))]
pd[which(pd$Library=="polyA"),"leftReads"] = gsub("/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/FASTQ/PolyA_Separation/",
                                                  "/dcl01/lieber/ajaffe/Brain/CytosolNucleus/polyA/FASTQ/", 
                                                  pd[which(pd$Library=="polyA"),"leftReads"])
pd[which(pd$Library=="polyA"),"rightReads"] = gsub("/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/FASTQ/PolyA_Separation/",
                                                  "/dcl01/lieber/ajaffe/Brain/CytosolNucleus/polyA/FASTQ/", 
                                                  pd[which(pd$Library=="polyA"),"rightReads"])
pd[which(pd$Library=="RiboZero"),"leftReads"] = gsub("/dcs01/lieber/ajaffe/Brain/DLPFC_Ribozero/FASTQ/CytosolNucleus/",
                                                     "/dcl01/lieber/ajaffe/Brain/CytosolNucleus/RiboZero/FASTQ/", 
                                                     pd[which(pd$Library=="RiboZero"),"leftReads"])
pd[which(pd$Library=="RiboZero"),"rightReads"] = gsub("/dcs01/lieber/ajaffe/Brain/DLPFC_Ribozero/FASTQ/CytosolNucleus/",
                                                     "/dcl01/lieber/ajaffe/Brain/CytosolNucleus/RiboZero/FASTQ/", 
                                                     pd[which(pd$Library=="RiboZero"),"rightReads"])
pd[which(pd$Library=="polyA"),"bamFile"] = gsub("/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/FASTQ/PolyA_Separation/TopHat/",
                                                "/dcl01/lieber/ajaffe/Amanda/NucVsCyt/pipeline_output/polyA/HISAT2_out/", 
                                                pd[which(pd$Library=="polyA"),"bamFile"])
pd[which(pd$Library=="polyA"),"bamFile"] = gsub("/accepted_hits.bam",
                                                "_polyA_accepted_hits.sorted.bam", 
                                                pd[which(pd$Library=="polyA"),"bamFile"])
pd[which(pd$Library=="RiboZero"),"bamFile"] = gsub("/dcs01/lieber/ajaffe/Brain/DLPFC_Ribozero/FASTQ/CytosolNucleus/TopHat/",
                                                "/dcl01/lieber/ajaffe/Amanda/NucVsCyt/pipeline_output/RiboZero/HISAT2_out/", 
                                                pd[which(pd$Library=="RiboZero"),"bamFile"])
pd[which(pd$Library=="RiboZero"),"bamFile"] = gsub("/accepted_hits.bam",
                                                "_RiboZero_accepted_hits.sorted.bam", 
                                                pd[which(pd$Library=="RiboZero"),"bamFile"])
pd = pd[,which(colnames(pd)!="totalMapped" & colnames(pd)!="mitoMapped")]

save(geneCounts, geneCounts.down, geneMap, metrics, pd, 
     exonCounts, exonCounts.down, exonMap, jMap, jCountsP,jCountsR,jCountsD,
     file="./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Load and Process RPKM
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_Amanda_polya_n12.rda")
exonRpkmP = exonRpkm
geneRpkmP = geneRpkm
jRpkmP = jRpkm
jMapP = jMap
metricsP = metrics
txNumReadsP = txNumReads
txTpmP = txTpm
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_Amanda_ribozero_n12.rda")
exonRpkmR = exonRpkm
geneRpkmR = geneRpkm
jRpkmR = jRpkm
jMapR = jMap
metricsR = metrics
txNumReadsR = txNumReads
txTpmR = txTpm
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_Amanda_downsamp_n2.rda")
exonRpkmD = exonRpkm
geneRpkmD = geneRpkm
jRpkmD = jRpkm
jMapD = jMap
metricsD = metrics
txNumReadsD = txNumReads
txTpmD = txTpm

# Combine files
metrics = rbind(metricsP,metricsR,metricsD)
geneRpkm = cbind(geneRpkmP, geneRpkmR)
exonRpkm = cbind(exonRpkmP, exonRpkmR)
txNumReads = cbind(txNumReadsP, txNumReadsR)
txTpm = cbind(txTpmP,txTpmR)
geneRpkm.down = cbind(geneRpkm[,c(1:6,8,10:24)],geneRpkmD)
exonRpkm.down = cbind(exonRpkm[,c(1:6,8,10:24)],exonRpkmD)
txNumReads.down = cbind(txNumReads[,c(1:6,8,10:24)],txNumReadsD)
txTpm.down = cbind(txTpm[,c(1:6,8,10:24)],txTpmD)
colnames(exonRpkm.down) = colnames(geneRpkm.down) = colnames(txTpm.down) = colnames(txNumReads.down) = gsub("downsamp", "polyA", colnames(geneRpkm.down))
jMap = c(jMapP,jMapR,jMapD)
jMap = sort(unique(jMap))

rm(exonRpkmP, exonRpkmR,exonRpkmD,geneRpkmP, geneRpkmR,geneRpkmD,metricsP,metricsR,metricsD,txNumReadsP,txNumReadsR,txNumReadsD,
   txTpmP,txTpmR,txTpmD,jRpkm,jMapP,jMapR,jMapD)

# filter 
geneRpkm = geneRpkm[which(rowSums(geneRpkm) > 0),which(colnames(geneRpkm)!="Br1113N1_RiboZero")]
exonRpkm = exonRpkm[,which(colnames(exonRpkm)!="Br1113N1_RiboZero")]
geneRpkm.down = geneRpkm.down[which(rowSums(geneRpkm.down) > 0),which(colnames(geneRpkm.down)!="Br1113N1_RiboZero")]
exonRpkm.down = exonRpkm.down[,which(colnames(exonRpkm.down)!="Br1113N1_RiboZero")]
geneRpkm = geneRpkm[,order(colnames(geneRpkm))]
geneRpkm.down = geneRpkm.down[,order(colnames(geneRpkm.down))]
exonRpkm = exonRpkm[,order(colnames(exonRpkm))]
exonRpkm.down = exonRpkm.down[,order(colnames(exonRpkm.down))]

save(geneRpkm, geneRpkm.down, geneMap, metrics, pd, txMap, txNumReads, txNumReads.down, txTpm, txTpm.down,
     exonRpkm, exonRpkm.down, exonMap, jMap, jRpkmP,jRpkmR,jRpkmD,
     file="./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")
