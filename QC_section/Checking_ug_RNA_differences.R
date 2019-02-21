library(ggplot2)



df = openxlsx::read.xlsx("./Dropbox/sorted_figures/github_controlled/QC_section/data/suppTable5_separationPheno.xlsx", startRow=2,colNames=T)
br = unique(df$BrNum)
ratio = data.frame(BrNum = br, Age = c(rep.int("Adult", 3), rep.int("Prenatal", 3)))

for (i in 1:length(br)) {
  ratio[i,"Nuclear.Yield"] = df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Nucleus"),"Yield(ug/mg)"]
  
  ratio[i,"Cytoplasm.Yield"] = df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Cytosol"),"Yield(ug/mg)"]
  
  ratio[i,"Nuclear.Total"] = df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Nucleus"),"ugTotal"]
  
  ratio[i,"Cytoplasm.Total"] = df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Cytosol"),"ugTotal"]
  
  ratio[i,"ratio"] = df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Nucleus"),"Yield(ug/mg)"] / 
    (df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Nucleus"),"Yield(ug/mg)"]+df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Cytosol"),"Yield(ug/mg)"]) * 100
  
  ratio[i,"total.Perc"] = df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Nucleus"),"ugTotal"] / 
    (df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Nucleus"),"ugTotal"]+df[which(df$BrNum==ratio[i,"BrNum"] & df$Zone=="Cytosol"),"ugTotal"]) * 100
  
}
ratio = reshape2::melt(ratio, measure.vars=c("total.Perc","Nuclear.Total","Cytoplasm.Total"))
ratio$variable = gsub("total.Perc", "Percent Total (ug)", ratio$variable)
ratio$variable = gsub("Nuclear.Total", "Total Nuclear (ug)", ratio$variable)
ratio$variable = gsub("Cytoplasm.Total", "Total Cytoplasm (ug)", ratio$variable)
ratio$Fraction = ifelse(ratio$variable=="Total Nuclear (ug)", "Nucleus","Cytoplasm")
ratio$Group = paste0(ratio$Age, "\n", ratio$Fraction)
ratio$Group = factor(ratio$Group, levels = c("Adult\nCytoplasm", "Prenatal\nCytoplasm", "Adult\nNucleus", "Prenatal\nNucleus"))

pdf("./Dropbox/sorted_figures/github_controlled/QC_section/data/nuclear.RNA.proportion.pdf",width=3.5,height=3.25)
ggplot(ratio[which(ratio$variable=="Percent Total (ug)"),], 
       aes(x = Age, y = value, fill=Age)) + geom_boxplot(outlier.shape=NA) + geom_jitter() +
  xlab("") + ylab("") + scale_fill_brewer(palette="Set1") + ylim(0,15) + ylab("Percent") +
  ggtitle("Percent Total RNA\nin Nucleus") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
ggplot(ratio[which(ratio$variable=="Total Nuclear (ug)"),], 
       aes(x = Age, y = value, fill=Age)) + geom_boxplot(outlier.shape=NA) + geom_jitter() +
  xlab("") + ylab("") + scale_fill_brewer(palette="Set1") + ylab("ug") +
  ggtitle("Total RNA\nin Nucleus") + ylim(0,1.2) +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
ggplot(ratio[which(ratio$variable=="Total Cytoplasm (ug)"),], 
       aes(x = Age, y = value, fill=Age)) + geom_boxplot(outlier.shape=NA) + geom_jitter() +
  xlab("") + ylab("") + scale_fill_brewer(palette="Set1") + ylab("ug") + ylim(0,13) +
  ggtitle("Total RNA\nin Cytoplasm") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
dev.off()

pdf("./Dropbox/sorted_figures/github_controlled/QC_section/data/RNA.proportion.pdf",width=5.5,height=3.25)
ggplot(ratio[which(ratio$variable!="Percent Total (ug)"),], 
       aes(x = Group, y = value, fill=Group)) + geom_boxplot(outlier.shape=NA) + geom_jitter() +
  xlab("") + ylab("") + ylim(0,15) + ylab("ug") +
  ggtitle("Total RNA") + theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
dev.off()

res = list(Percent.Total = t.test(ratio[which(ratio$variable=="Percent Total (ug)" & ratio$Age=="Adult"),"value"], 
                                  ratio[which(ratio$variable=="Percent Total (ug)" & ratio$Age=="Prenatal"),"value"]),
           Nuclear.ug = t.test(ratio[which(ratio$variable=="Total Nuclear (ug)" & ratio$Age=="Adult"),"value"], 
                               ratio[which(ratio$variable=="Total Nuclear (ug)" & ratio$Age=="Prenatal"),"value"]),
           Cytoplasm.ug = t.test(ratio[which(ratio$variable=="Total Cytoplasm (ug)" & ratio$Age=="Adult"),"value"], 
                                 ratio[which(ratio$variable=="Total Cytoplasm (ug)" & ratio$Age=="Prenatal"),"value"]))

res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, mean.Adult = x$estimate[1], 
                                                            mean.Prenatal = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#           Group       Tstat mean.NucGroup mean.Other       pval        FDR
#t  Percent.Total  4.89079659    11.5913746  8.1255509 0.03024013 0.08599446
#t1    Nuclear.ug  0.09723657     0.8953333  0.8896667 0.92858536 0.92858536
#t2  Cytoplasm.ug -3.00823482     6.8458333 10.2145000 0.05732964 0.08599446


