library(GenomicRanges)
library(RColorBrewer)
library(jaffelab)

load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("./Dropbox/sorted_figures/github_controlled/disease/data/gwasResults_lifted.rda")


## Enrichment in GWAS genes and gene sets

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = ss(pgc$Position..hg19.,":")
tmp = ss(pgc$Position..hg19.,":",2)
pgc$start = as.numeric(ss(tmp, "-"))
pgc$end = as.numeric(ss(tmp, "-",2))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
pgcGR$Dx = "SCZ"

gwas = c(as.list(split(gwasLift, gwasLift$Dx)), list(SCZ = pgcGR))

geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
genes = lapply(gwas, function(x) findOverlaps(geneMapGR, x))
genes = lapply(genes, function(x) geneMapGR[queryHits(x)])
gwasGenes = lapply(genes, function(x) unique(as.character(na.omit(x$EntrezID))))

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')


## Enrichment in genes differentially expressed by fraction

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"EntrezID"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$EntrezGene.ID) %in% geneuniverse ), ] # drop genes that are not present in the test set
aej_sets_expressed$EntrezGene.ID = as.character(aej_sets_expressed$EntrezGene.ID)

splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
splitSets = c(gwasGenes, lapply(splitSets, function(x) as.character(na.omit(unique(x$EntrezGene.ID)))))
splitSets = splitSets[which(names(splitSets)!="SCZ PGC GWAS")]

inGroup = lapply(sig, function(x) as.character(na.omit(unique(x$EntrezID))))
outGroup = lapply(sig, function(x) geneuniverse[!(geneuniverse %in% as.character(na.omit(x$EntrezID)))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x),sum(!(inG %in% x)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x), sum(!(outG %in% x)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)
enrich = lapply(enrich, data.frame)
enrich = do.call(rbind, Map(cbind, Comparison = as.list(names(enrich)), lapply(enrich, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")

write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/github_controlled/disease/data/Birnbaum_geneSet_GWAS_enrichment_FractionDEGs.csv")
enrich = read.csv("./Dropbox/sorted_figures/github_controlled/disease/data/Birnbaum_geneSet_GWAS_enrichment_FractionDEGs.csv")

enrich[enrich$FDR<=0.05,]







