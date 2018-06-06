library(GenomicRanges)
library(data.table)

load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/introns_object.rda")

### How many genes are represented by how many retained introns?
lapply(introns, head)
introns = c(introns, list(Cytoplasmic = rbind(introns[["Adult:Cytoplasm-Increased"]],introns[["Prenatal:Cytoplasm-Increased"]]),
                          Nuclear = rbind(introns[["Adult:Nucleus-Increased"]],introns[["Prenatal:Nucleus-Increased"]]),
                          Adult = rbind(introns[["Cytoplasm:Adult-Increased"]],introns[["Nucleus:Adult-Increased"]]),
                          Prenatal = rbind(introns[["Cytoplasm:Prenatal-Increased"]],introns[["Nucleus:Prenatal-Increased"]])))
genes = data.frame(numgenes = unlist(lapply(introns, function(x) length(unique(x$ensID)))),
                   numintrons = unlist(lapply(introns, function(x) length(unique(x$intronID)))))
genes$group = c("All", "Cytoplasmic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytoplasmic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                "Adult Retention\nIn Cytoplasm","Prenatal Retention\nIn Cytoplasm","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus",
                "Pooled Cytoplasmic\nRetention","Pooled Nuclear\nRetention","Pooled Adult\nRetention","Pooled Prenatal\nRetention")
genes$retention = c("NA","Cytoplasmic","Nuclear","Cytoplasmic","Nuclear","Adult","Prenatal","Adult","Prenatal","Cytoplasmic","Nuclear","Adult","Prenatal")

## by Fraction: Cytoplasmic (is there a relationship between the number of introns found dIR and the number of genes?)

tables = list()
for (i in 2:nrow(genes)) {
  tables[[(i-1)]] = fisher.test(data.frame(c(genes[i,1],genes[i,2]),c(genes[1,1]-genes[i,1],genes[1,2]-genes[i,2])))
}
names(tables) = rownames(genes)[2:nrow(genes)]
tables = c(tables, list("Frac:Nuc.VS.Cyt" = fisher.test(data.frame(c(genes["Cytoplasmic",1],genes["Cytoplasmic",2]),c(genes["Nuclear",1],genes["Nuclear",2]))),
                        "Adult:Nuc.VS.Cyt" = fisher.test(data.frame(c(genes["Adult:Cytoplasm-Increased",1],genes["Adult:Cytoplasm-Increased",2]),
                                                                    c(genes["Adult:Nucleus-Increased",1],genes["Adult:Nucleus-Increased",2]))),
                        "Prenatal:Nuc.VS.Cyt" = fisher.test(data.frame(c(genes["Prenatal:Cytoplasm-Increased",1],genes["Prenatal:Cytoplasm-Increased",2]),
                                                                       c(genes["Prenatal:Nucleus-Increased",1],genes["Prenatal:Nucleus-Increased",2]))),
                        "Age:Adult.VS.Pren" = fisher.test(data.frame(c(genes["Adult",1],genes["Adult",2]),c(genes["Prenatal",1],genes["Prenatal",2]))),
                        "Cyt:Adult.VS.Pren" = fisher.test(data.frame(c(genes["Cytoplasm:Adult-Increased",1],genes["Cytoplasm:Adult-Increased",2]),
                                                                     c(genes["Cytoplasm:Prenatal-Increased",1],genes["Cytoplasm:Prenatal-Increased",2]))),
                        "Nuc:Adult.VS.Pren" = fisher.test(data.frame(c(genes["Nucleus:Adult-Increased",1],genes["Nucleus:Adult-Increased",2]),
                                                                     c(genes["Nucleus:Prenatal-Increased",1],genes["Nucleus:Prenatal-Increased",2])))))

df = data.frame(Comparison = names(tables), pval = unlist(lapply(tables, function(x) x$p.value)), OddsRatio = unlist(lapply(tables, function(x) x$estimate)), row.names = NULL)
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ratio_retainedIntrons_to_uniqueGenes.csv")
df[df$FDR<=0.05,]
#                     Comparison         pval OddsRatio          FDR
#2       Adult:Nucleus-Increased 6.563496e-15  8.317378 5.907146e-14
#4    Prenatal:Nucleus-Increased 2.335464e-03  7.813589 7.006393e-03
#6  Cytoplasm:Prenatal-Increased 8.126824e-05  8.334628 2.925657e-04
#8    Nucleus:Prenatal-Increased 3.715396e-08  7.901259 1.671928e-07
#10                      Nuclear 9.188406e-17  8.429571 1.653913e-15
#11                        Adult 1.432685e-02  9.373242 3.684047e-02
#12                     Prenatal 5.033587e-10  8.160341 3.020152e-09

## How many introns per gene represented?
pergeneIR = lapply(introns, function(x) data.frame(numIR=data.table(x)[,length(unique(intronID)), by="ensID"]))
stat = data.frame(mean = unlist(lapply(pergeneIR, function(x) mean(x$numIR.V1))), sd = unlist(lapply(pergeneIR,function(x) sd(x$numIR.V1))))
morethan1 = lapply(pergeneIR, function(x) x[which(x$numIR.V1>1),])
elementNROWS(morethan1)
elementNROWS(pergeneIR)

## number of genes with more than one retained intron
tables = list()
for (i in 2:length(pergeneIR)) {
  tables[[(i-1)]] = fisher.test(data.frame(c(elementNROWS(morethan1)[i],elementNROWS(pergeneIR)[i]-elementNROWS(morethan1)[i]),
                                           c(elementNROWS(morethan1)[1]-elementNROWS(morethan1)[i],
                                             (elementNROWS(pergeneIR)[1]-elementNROWS(morethan1)[1])-(elementNROWS(pergeneIR)[i]-elementNROWS(morethan1)[i]))))
}
names(tables) = names(pergeneIR)[2:length(pergeneIR)]
tables = c(tables, list("Adult:Nuc.VS.Cyt" = fisher.test(data.frame(c(elementNROWS(morethan1)["Adult:Cytoplasm-Increased"],
                                                                      elementNROWS(pergeneIR)["Adult:Cytoplasm-Increased"]-elementNROWS(morethan1)["Adult:Cytoplasm-Increased"]),
                                                                    c(elementNROWS(morethan1)["Adult:Nucleus-Increased"],
                                                                      elementNROWS(pergeneIR)["Adult:Nucleus-Increased"]-elementNROWS(morethan1)["Adult:Nucleus-Increased"]))),
                        "Prenatal:Nuc.VS.Cyt" = fisher.test(data.frame(c(elementNROWS(morethan1)["Prenatal:Cytoplasm-Increased"],
                                                                         elementNROWS(pergeneIR)["Prenatal:Cytoplasm-Increased"]-elementNROWS(morethan1)["Prenatal:Cytoplasm-Increased"]),
                                                                       c(elementNROWS(morethan1)["Prenatal:Nucleus-Increased"],
                                                                         elementNROWS(pergeneIR)["Prenatal:Nucleus-Increased"]-elementNROWS(morethan1)["Prenatal:Nucleus-Increased"]))),
                        "Cyt:Adult.VS.Pren" = fisher.test(data.frame(c(elementNROWS(morethan1)["Cytoplasm:Adult-Increased"],
                                                                       elementNROWS(pergeneIR)["Cytoplasm:Adult-Increased"]-elementNROWS(morethan1)["Cytoplasm:Adult-Increased"]),
                                                                     c(elementNROWS(morethan1)["Cytoplasm:Prenatal-Increased"],
                                                                       elementNROWS(pergeneIR)["Cytoplasm:Prenatal-Increased"]-elementNROWS(morethan1)["Cytoplasm:Prenatal-Increased"]))),
                        "Nuc:Adult.VS.Pren" = fisher.test(data.frame(c(elementNROWS(morethan1)["Nucleus:Adult-Increased"],
                                                                       elementNROWS(pergeneIR)["Nucleus:Adult-Increased"]-elementNROWS(morethan1)["Nucleus:Adult-Increased"]),
                                                                     c(elementNROWS(morethan1)["Nucleus:Prenatal-Increased"],
                                                                       elementNROWS(pergeneIR)["Nucleus:Prenatal-Increased"]-elementNROWS(morethan1)["Nucleus:Prenatal-Increased"])))))
df = data.frame(Comparison = names(tables), pval = unlist(lapply(tables, function(x) x$p.value)), OddsRatio = unlist(lapply(tables, function(x) x$estimate)), row.names = NULL)
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ratio_retainedIntrons_per_uniqueGenes.csv")
df[df$FDR<=0.05,]
#                    Comparison         pval  OddsRatio          FDR
#2      Adult:Nucleus-Increased 8.742442e-22 0.01913714 1.049093e-20
#4   Prenatal:Nucleus-Increased 8.168124e-04 0.03273576 2.450437e-03
#6 Cytoplasm:Prenatal-Increased 2.000103e-06 0.01863621 8.000410e-06
#7      Nucleus:Adult-Increased 1.343243e-02 0.00000000 3.223782e-02
#8   Nucleus:Prenatal-Increased 2.608347e-10 0.03004254 1.565008e-09

x = do.call(rbind, pergeneIR[which(elementNROWS(pergeneIR)>0 & elementNROWS(pergeneIR)<1000)])
length(unique(x[which(x$numIR.V1==1),"numIR.ensID"]))/length(unique(x$numIR.ensID))
# 0.9038462