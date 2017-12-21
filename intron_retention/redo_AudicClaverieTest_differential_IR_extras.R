library(GenomicRanges)



paste0("mkdir ",pd$SampleID)
paste0("scp ./Dropbox/sorted_figures/IRfinder/PolyA/",pd[pd$Fetal=="Adult","SampleID"],
       "/IRFinder-IR-nondir-cleanIntrons-adultShared.txt aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/",
       pd[pd$Fetal=="Adult","SampleID"],"/")

paste0("scp ./Dropbox/sorted_figures/IRfinder/PolyA/",pd[pd$Fetal=="Prenatal","SampleID"],
       "/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/",
       pd[pd$Fetal=="Prenatal","SampleID"],"/")

paste0("scp ./Dropbox/sorted_figures/IRfinder/PolyA/",pd[pd$Zone=="Cytosol","SampleID"],
       "/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/",
       pd[pd$Zone=="Cytosol","SampleID"],"/")

paste0("scp ./Dropbox/sorted_figures/IRfinder/PolyA/",pd[pd$Zone=="Nucleus","SampleID"],
       "/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/",
       pd[pd$Zone=="Nucleus","SampleID"],"/")


## Filter pooled samples to match the results files that were copied from Dropbox from IRFinder_GLM.r

pooled = list(APC_Pooled = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/IRfinder_redo/Filtered/PolyA/APC_Pooled/IRFinder-IR-nondir.txt", header=T),
              FPC_Pooled = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/IRfinder_redo/PolyA/FPC_Pooled/IRFinder-IR-nondir.txt", header=T),
              APN_Pooled = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/IRfinder_redo/Filtered/PolyA/APN_Pooled/IRFinder-IR-nondir.txt", header=T),
              FPN_Pooled = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/IRfinder_redo/Filtered/PolyA/FPN_Pooled/IRFinder-IR-nondir.txt", header=T))
pooled = Map(cbind, pooled, intronID = lapply(pooled, function(x) paste0(x$GeneIntronDetails,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)))

filtered = list(adult = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/Br2046C/IRFinder-IR-nondir-cleanIntrons-adultShared.txt", header=T),
                prenatal = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/Br5341N1/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt", header=T),
                cyt = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/Br2046C/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt", header=T),
                nuc = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/Br5341N1/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt", header=T))
filtered = Map(cbind, filtered, intronID = lapply(filtered, function(x) paste0(x$GeneIntronDetails,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)))

path = "/dcl01/lieber/ajaffe/Amanda/NucVsCyt/11.29.17/"
write.table(pooled$APC_Pooled[pooled$APC_Pooled$intronID %in% filtered$adult$intronID,], 
            file=paste0(path,"APC_Pooled/IRFinder-IR-nondir-cleanIntrons-adultShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$APC_Pooled[pooled$APC_Pooled$intronID %in% filtered$cyt$intronID,], 
            file=paste0(path,"APC_Pooled/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$FPC_Pooled[pooled$FPC_Pooled$intronID %in% filtered$prenatal$intronID,], 
            file=paste0(path,"FPC_Pooled/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$FPC_Pooled[pooled$FPC_Pooled$intronID %in% filtered$cyt$intronID,], 
            file=paste0(path,"FPC_Pooled/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$APN_Pooled[pooled$APN_Pooled$intronID %in% filtered$adult$intronID,], 
            file=paste0(path,"APN_Pooled/IRFinder-IR-nondir-cleanIntrons-adultShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$APN_Pooled[pooled$APN_Pooled$intronID %in% filtered$nuc$intronID,], 
            file=paste0(path,"APN_Pooled/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$FPN_Pooled[pooled$FPN_Pooled$intronID %in% filtered$prenatal$intronID,], 
            file=paste0(path,"FPN_Pooled/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt"), row.names = F, quote = F, sep = "\t")
write.table(pooled$FPN_Pooled[pooled$FPN_Pooled$intronID %in% filtered$nuc$intronID,], 
            file=paste0(path,"FPN_Pooled/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt"), row.names = F, quote = F, sep = "\t")





