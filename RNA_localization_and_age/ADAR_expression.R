load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")


# Check ADAR expression

adars = geneMap[grep("ADAR", geneMap$Symbol),]

Apres[rownames(adars),] # none significant

Fpres[rownames(adars),] # none significant

ADAR = rbind(cbind(as.data.frame(Apres[rownames(adars),1:6]), 
                   Symbol = adars[match(rownames(Apres[rownames(adars),]),adars$gencodeID),"Symbol"],
                   gencodeID = adars[match(rownames(Apres[rownames(adars),]),adars$gencodeID),"gencodeID"],
                   Comparison = "By Fraction in Adult"),
             cbind(as.data.frame(Fpres[rownames(adars),]), 
                   Symbol = adars[match(rownames(Fpres[rownames(adars),]),adars$gencodeID),"Symbol"],
                   gencodeID = adars[match(rownames(Fpres[rownames(adars),]),adars$gencodeID),"gencodeID"],
                   Comparison = "By Fraction in Prenatal"),
             cbind(as.data.frame(Cpres[rownames(adars),]), 
                   Symbol = adars[match(rownames(Cpres[rownames(adars),]),adars$gencodeID),"Symbol"],
                   gencodeID = adars[match(rownames(Cpres[rownames(adars),]),adars$gencodeID),"gencodeID"],
                   Comparison = "By Age in Cytosol"),
             cbind(as.data.frame(Npres[rownames(adars),]), 
                   Symbol = adars[match(rownames(Npres[rownames(adars),]),adars$gencodeID),"Symbol"],
                   gencodeID = adars[match(rownames(Npres[rownames(adars),]),adars$gencodeID),"gencodeID"],
                   Comparison = "By Age in Nucleus"))

write.csv(ADAR, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/ADAR_gene_expression.csv", quote = F, row.names = F)

ADAR[which(ADAR$padj<=0.05),]
# ADAR is not DE by fraction or age, 
# but ADARB2, ADARB2-AS1, ADARB1 are greater in adult than prenatal in cytosol, 
# and ADARB2 and ADARB1 are greater in adult than prenatal in nucleus.


