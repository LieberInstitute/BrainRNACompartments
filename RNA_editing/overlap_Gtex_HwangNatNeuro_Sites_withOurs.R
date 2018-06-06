library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


## load public editing site lists
gtex_all = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/GTEX_Supp1_allSites.xlsx')
gtex_all = cbind(gtex_all, seqnames = unlist(strsplit(gtex_all$sites,"_",fixed=T))[grep("chr",unlist(strsplit(gtex_all$sites,"_",fixed=T)))],
                 start = unlist(strsplit(gtex_all$sites,"_",fixed=T))[-grep("chr",unlist(strsplit(gtex_all$sites,"_",fixed=T)))], site = gtex_all$sites)
gtex_all = makeGRangesFromDataFrame(gtex_all, end.field = "start", keep.extra.columns = T)

gtex_ts = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/GTEX_Supp3_tissueSpecificSites.xlsx')
gtex_ts = makeGRangesFromDataFrame(gtex_ts[,1:4], seqnames.field="Chromosome",start.field="Position", end.field="Position",keep.extra.columns = T)

increasing = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/Hwang_NatNeuro_SuppTables.xlsx',sheet =3)
increasing = increasing[-1,]
incr = makeGRangesFromDataFrame(increasing, seqnames.field = "Chromosome",start.field = "Coordinate",end.field = "Coordinate", strand.field = "Strand")

editing_anno$start = editing_anno$end
AG = editing_anno[collapsedconversion=="A:G / T:C",,]
editing_annogr = makeGRangesFromDataFrame(editing_anno)
AGgr = makeGRangesFromDataFrame(AG)
reduced = unique(granges(editing_annogr)) # 25045
AGred = unique(granges(AGgr)) # 18907


## How many of our sites have been identified already?
red_all = findOverlaps(AGred, unique(granges(gtex_all)))
length(unique(queryHits(red_all))) # 13054 A to I sites are identified in the Gtex samples out of 18907 
red_ts = findOverlaps(AGred, unique(granges(gtex_ts)))
length(unique(queryHits(red_ts))) # 66 are identified as tissue-specific in GTEX 
ov_inc = findOverlaps(AGred, incr)
length(unique(queryHits(ov_inc))) # 576 of the 742 increasing editing sites identified in Taeyoung's paper are in our set


## Annotate our sites according to the GTEX data

editing_anno$id = paste0(editing_anno$seqnames, ":",editing_anno$start, "-", editing_anno$end)
ov_all = findOverlaps(makeGRangesFromDataFrame(editing_anno), gtex_all)
editing_anno = rbind(cbind(editing_anno[queryHits(ov_all),,], new = "In GTEX"),cbind(editing_anno[-unique(queryHits(ov_all)),,], new = "Novel"))
ov_ts = findOverlaps(makeGRangesFromDataFrame(editing_anno), gtex_ts)
editing_anno = rbind(cbind(editing_anno[queryHits(ov_ts),,], tissuespecific = "Tissue-Specific"),
                     cbind(editing_anno[-unique(queryHits(ov_ts)),,], tissuespecific = "Multi-Tissue"))
ov_inc = findOverlaps(makeGRangesFromDataFrame(editing_anno), incr)
editing_anno = rbind(cbind(editing_anno[queryHits(ov_inc),,], Increasing = "Increasing"),
                     cbind(editing_anno[-unique(queryHits(ov_inc)),,], Increasing = "Not Increasing"))


editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by=c("Fraction","new")]
#   Fraction     new    V1
#1:  Nucleus In GTEX 10593
#2:  Nucleus   Novel  4484
#3:  Cytosol   Novel  3133
#4:  Cytosol In GTEX  6965
(4484-3133)/3133*100 # 43.12161
editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by=c("Age","new")]
#        Age     new   V1
#1:    Adult In GTEX 9856
#2: Prenatal In GTEX 5952
#3:    Adult   Novel 3002
#4: Prenatal   Novel 3416
(3416-3002)/3002*100 # 13.79081
editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by="Fraction"]
#   Fraction    V1
#1:  Cytosol 10098
#2:  Nucleus 15077
(15077-10098)/10098*100 # 49.30679

gtex = list()
mcols(gtex_all) = mcols(gtex_all)[,!colnames(mcols(gtex_all)) %in% c("site", "sites")]
for (i in 1:53) {
  gtex[[i]] = granges(gtex_all)
  mcols(gtex[[i]])$present = mcols(gtex_all)[,i]
}
names(gtex) = colnames(mcols(gtex_all))
gtex = lapply(gtex, function(x) x[!is.na(x$present)])
gtex = lapply(gtex, function(x) x[x$present>0])

ov = lapply(gtex, function(x) findOverlaps(reduced, unique(granges(x))))
min(elementNROWS(lapply(ov, function(x) unique(subjectHits(x))))/elementNROWS(gtex)*100) #5.015265
df = data.frame(tissue = names(gtex),
                gtexHits = elementNROWS(lapply(ov, function(x) unique(subjectHits(x)))), 
                percentGTEX = elementNROWS(lapply(ov, function(x) unique(subjectHits(x))))/elementNROWS(gtex)*100,
                percentOURS = elementNROWS(lapply(ov, function(x) unique(queryHits(x))))/length(reduced)*100,row.names = NULL)
df[grep("Brain", df$tissue),"Brain"] = "Brain"
df[-grep("Brain", df$tissue),"Brain"] = "Other"
write.csv(df, quote=F,file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/GTEX_editingSites_inOurData_byTissue.csv")

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/GTEX_editingSites_inOurData_byTissue.pdf",width=18,height=8)
df$tissue = factor(df$tissue, levels = as.character(df[order(df$percentOURS, decreasing = T),"tissue"]))
ggplot(df, aes(x = tissue, y = percentOURS, fill=Brain)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  labs(fill="") + ylab("Percent") + xlab("") +
  scale_fill_manual(values=c("blue1","gray30")) +
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) +
  ggtitle("Percent of Our Editing Sites Found in Each GTEX Tissue") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
df$tissue = factor(df$tissue, levels = as.character(df[order(df$percentGTEX, decreasing = T),"tissue"]))
ggplot(df, aes(x = tissue, y = percentGTEX, fill=Brain)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  labs(fill="") + ylab("Percent") + xlab("") +
  scale_fill_manual(values=c("blue1","gray30")) +
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) +
  ggtitle("Percent of GTEX Editing Sites Found in Each Tissue in Our Study") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Plot the distribution

editing_anno$Fraction = gsub("Cytosol","Cytoplasm", editing_anno$Fraction)
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/GTEX_HwangNatNeuro_editingSites_inOurData.pdf",width=8,height=5)
x = editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by=c("Fraction","Age","new")]
x$perc = c(round(x$V1[1]/sum(x$V1[c(1,5)])*100,1),round(x$V1[c(2)]/sum(x$V1[c(2,6)])*100,1),round(x$V1[c(3)]/sum(x$V1[c(3,7)])*100,1),
           round(x$V1[c(4)]/sum(x$V1[c(4,8)])*100,1),round(x$V1[5]/sum(x$V1[c(1,5)])*100,1),round(x$V1[c(6)]/sum(x$V1[c(2,6)])*100,1),
           round(x$V1[c(7)]/sum(x$V1[c(3,7)])*100,1),round(x$V1[c(8)]/sum(x$V1[c(4,8)])*100,1))
ggplot(x, aes(x = Fraction, y = V1, fill = new)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  facet_grid(. ~ Age) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to GTEX") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by=c("Fraction","Age","tissuespecific")]
x$perc = c(round(x$V1[1]/sum(x$V1[c(1,3)])*100,1),round(x$V1[c(2)]/sum(x$V1[c(2,4)])*100,1),round(x$V1[c(3)]/sum(x$V1[c(1,3)])*100,1),
           round(x$V1[4]/sum(x$V1[c(2,4)])*100,1),round(x$V1[c(5)]/sum(x$V1[c(5,7)])*100,1),round(x$V1[c(6)]/sum(x$V1[c(6,8)])*100,1),
           round(x$V1[c(7)]/sum(x$V1[c(5,7)])*100,1),round(x$V1[c(8)]/sum(x$V1[c(6,8)])*100,1))
ggplot(x[tissuespecific=="Tissue-Specific",,], aes(x = Age, y = V1, fill = Fraction)) + 
  geom_bar(stat = "identity",position=position_dodge(width=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") + scale_fill_brewer(palette="Dark2") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Tissue-Specific Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by=c("Fraction","Age","Increasing")]
x$perc = c(round(x$V1[1:2]/sum(x$V1[1:2])*100,1),round(x$V1[3:4]/sum(x$V1[3:4])*100,1),round(x$V1[5:6]/sum(x$V1[5:6])*100,1),round(x$V1[7:8]/sum(x$V1[7:8])*100,1))
ggplot(x[Increasing=="Increasing",,], aes(x = Age, y = V1, fill = Fraction)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") + scale_fill_brewer(palette="Dark2") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to\nDevelopmentally-Increasing Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

## Isolate the unique sites present in all samples in each group

unique_bySamp_all = list(cytosolOnly = unique_bySamp$cytosolOnly[-grep("no", unique_bySamp$cytosolOnly$cytosolAll),],
                         nucleusOnly = unique_bySamp$nucleusOnly[-grep("no", unique_bySamp$nucleusOnly$nucleusAll),], 
                         adultOnly = unique_bySamp$adultOnly[-grep("no", unique_bySamp$adultOnly$adultAll),], 
                         prenatalOnly = unique_bySamp$prenatalOnly[-grep("no", unique_bySamp$prenatalOnly$prenatalAll),], 
                         ANnotAC = unique_bySamp$ANnotAC[-grep("no", unique_bySamp$ANnotAC$allAN),], 
                         ACnotAN = unique_bySamp$ACnotAN[-grep("no", unique_bySamp$ACnotAN$allAC),], 
                         ANnotPN = unique_bySamp$ANnotPN[-grep("no", unique_bySamp$ANnotPN$allAN),], 
                         PNnotAN = unique_bySamp$PNnotAN[-grep("no", unique_bySamp$PNnotAN$allPN),],
                         ACnotPC = unique_bySamp$ACnotPC[-grep("no", unique_bySamp$ACnotPC$allAC),], 
                         PCnotAC = unique_bySamp$PCnotAC[-grep("no", unique_bySamp$PCnotAC$allPC),], 
                         PCnotPN = unique_bySamp$PCnotPN[-grep("no", unique_bySamp$PCnotPN$allPC),], 
                         PNnotPC = unique_bySamp$PNnotPC[-grep("no", unique_bySamp$PNnotPC$allPN),])
unique_all = lapply(unique_bySamp_all[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = Map(cbind, unique_all, ensID = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]), name = as.list(names(unique_all)))
elementNROWS(lapply(unique_all, function(x) unique(x$id)))
unique_all = do.call(rbind,unique_all)
unique_all$name = factor(unique_all$name, levels= c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))
unique_all = as.data.frame(unique_all)[,colnames(unique_all) %in% c("seqnames","start","end","strand","ref","alt","depth","valdepth","ref.count","alt.count","sampleID",
                                                                    "editingID","Fraction","Age","rate","overlappingGene","annotation","nearestSymbol","nearestID",
                                                                    "distToGene","EntrezID","ensID","Type","name")]
unique_all$name = gsub("adultOnly", "Adult Only", unique_all$name)
unique_all$name = gsub("prenatalOnly", "Prenatal Only",unique_all$name)
unique_all$name = gsub("ANnotAC", "AN not AC",unique_all$name)
unique_all$name = gsub("ACnotAN", "AC not AN", unique_all$name)
unique_all$name = gsub("ANnotPN", "AN not PN", unique_all$name)
unique_all$name = gsub("PNnotAN", "PN not AN", unique_all$name)
unique_all$name = gsub("ACnotPC", "AC not PC", unique_all$name)
unique_all$name = gsub("PCnotAC", "PC not AC", unique_all$name)
unique_all$name = gsub("PCnotPN", "PC not PN", unique_all$name)
unique_all$name = gsub("PNnotPC", "PN not PC", unique_all$name)
write.csv(unique_all, quote = F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/table4_editingSites_uniqueAll3.csv") 

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/GTEX_HwangNatNeuro_editingSites_inOurData_uniqueAll3.pdf",width=18,height=8)
x = unique_all[,length(unique(id)),by=c("name","new")]
for (i in 1:length(unique(x$name))){ 
  x$perc[grep(unique(x$name)[i],x$name)] = round(x$V1[grep(unique(x$name)[i],x$name)]/sum(x$V1[grep(unique(x$name)[i],x$name)])*100,1)}
ggplot(x, aes(x = new, y = V1, fill = new)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  facet_grid(. ~ name) +
  labs(fill="") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to GTEX") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
x = unique_all[,length(unique(id)),by=c("name","tissuespecific")]
for (i in 1:length(unique(x$name))){ 
  x$perc[grep(unique(x$name)[i],x$name)] = round(x$V1[grep(unique(x$name)[i],x$name)]/sum(x$V1[grep(unique(x$name)[i],x$name)])*100,1)}
ggplot(x, aes(x = tissuespecific, y = V1, fill = tissuespecific)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  facet_grid(. ~ name) +
  labs(fill="") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to GTEX") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
x = unique_all[,length(unique(id)),by=c("name","Increasing")]
for (i in 1:length(unique(x$name))){ 
  x$perc[grep(unique(x$name)[i],x$name)] = round(x$V1[grep(unique(x$name)[i],x$name)]/sum(x$V1[grep(unique(x$name)[i],x$name)])*100,1)}
ggplot(x, aes(x = Increasing, y = V1, fill = Increasing)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  facet_grid(. ~ name) +
  labs(fill="") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to\nDevelopmentally-Increasing Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())
dev.off()


## Which Tissue-specific sites in the previous datasets are included in ours?
gtex_ts = data.table(as.data.frame(unique(gtex_ts[subjectHits(ov_ts)])))
gtex_ts$ID = paste0(gtex_ts$seqnames, ":",gtex_ts$start, "-", gtex_ts$end)
gtex_ts = gtex_ts[,length(unique(ID)),by="Specifically.edited.tissue"]
gtex_ts$prop = round(gtex_ts$V1/sum(gtex_ts$V1)*100,1)

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/GTEX_TissueSpecific_Annotation.pdf")
ggplot(gtex_ts, aes(x = Specifically.edited.tissue, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(prop,"%")), vjust = -.5) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Tissue-Specific GTEX Editing Sites\nIdentified in Our Data") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Are more of our unique-in-all-samples sites developmentally increasing in Taeyoung's paper than expected by chance?

al = as.data.frame(editing_anno[collapsedconversion=="A:G / T:C",length(unique(id)),by="Increasing"])
un = as.data.frame(unique_all[collapsedconversion=="A:G / T:C",length(unique(id)),by=c("name","Increasing")])
names = as.character(unique(unique_all$name))
inc = list()
for (i in 1:length(names)){ 
  if (nrow(un[un$name==names[i] & un$Increasing=="Increasing",])>0) {
    inc[[i]] = data.frame(inanno = c(un[un$name==names[i] & un$Increasing=="Increasing","V1"], 
                                     un[un$name==names[i] & un$Increasing=="Not Increasing","V1"]),
                          notanno = c(al[al$Increasing=="Increasing","V1"]-un[un$name==names[i] & un$Increasing=="Increasing","V1"], 
                                      al[al$Increasing=="Not Increasing","V1"]-un[un$name==names[i] & un$Increasing=="Not Increasing","V1"]),
                          row.names = c("Increasing","Not"))
  } else {
    inc[[i]] = data.frame(inanno = c(0, un[un$name==names[i] & un$Increasing=="Not Increasing","V1"]),
                          notanno = c(al[al$Increasing=="Increasing","V1"], al[al$Increasing=="Not Increasing","V1"]-un[un$name==names[i] & un$Increasing=="Not Increasing","V1"]),
                          row.names = c("Increasing","Not"))
  }
}
names(inc) = names
fisher = lapply(inc, fisher.test)
df = do.call(rbind, Map(cbind, Annotation = as.list(names(fisher)), lapply(fisher, function(x) data.frame(pval=x$p.value, 
                                                                                                          OddsRatio=x$estimate, row.names = NULL))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(do.call(rbind, Map(cbind, Annotation = as.list(names(inc)), lapply(inc, function(x) data.frame(Dir=rownames(x), x, row.names = NULL)))),quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/HwangNatNeuro_editingSites_contingencytables_inOurData.csv")
write.csv(df, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/HwangNatNeuro_editingSites_enrichment_inOurData.csv", quote = F)
df[df$FDR<=0.05,]
#  Annotation         pval OddsRatio          FDR
#1  adultOnly 7.423829e-26 13.455781 2.474610e-25
#5    ANnotPN 9.568469e-53  8.752737 9.568469e-52
#6    PNnotAN 3.289963e-04  0.000000 8.224907e-04
#7    ACnotPC 7.602388e-45 12.348638 3.801194e-44


## Are the increasing sites more edited in adult than prenatal?

ttests = list(all = t.test(x = unique(editing_anno[collapsedconversion=="A:G / T:C" & Increasing=="Increasing" & Age=="Adult",list(rate),by=c("id","sampleID")])$rate, 
                           y = unique(editing_anno[collapsedconversion=="A:G / T:C" & Increasing=="Increasing" & Age=="Prenatal",list(rate),by=c("id","sampleID")])$rate),
              nucleus = t.test(x = unique(editing_anno[collapsedconversion=="A:G / T:C" & Increasing=="Increasing" & Group=="Adult:Nucleus",list(rate),by=c("id","sampleID")])$rate, 
                               y = unique(editing_anno[collapsedconversion=="A:G / T:C" & Increasing=="Increasing" & Group=="Prenatal:Nucleus",list(rate),by=c("id","sampleID")])$rate),
              cytosol = t.test(x = unique(editing_anno[collapsedconversion=="A:G / T:C" & Increasing=="Increasing" & Group=="Adult:Cytosol",list(rate),by=c("id","sampleID")])$rate, 
                               y = unique(editing_anno[collapsedconversion=="A:G / T:C" & Increasing=="Increasing" & Group=="Prenatal:Cytosol",list(rate),by=c("id","sampleID")])$rate))
df = do.call(rbind, Map(cbind, In = as.list(names(ttests)), 
                        lapply(ttests, function(x) data.frame(Tstat = x$statistic,pval = x$p.value, confInt1 = x$conf.int[1], 
                                                              confInt2 = x$conf.int[2], estMeans1 = x$estimate[1], estMeans2 = x$estimate[2]))))
df$FDR = p.adjust(df$pval, method="fdr")
write.csv(df,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/HwangNatNeuro_editingRate_ttest_byAge_inIncreasingSites.csv",quote=F)
df[df$FDR<=0.05,]
#        In   Tstat      pval  confInt1    confInt2 estMeans1 estMeans2        FDR
#t2 cytosol -2.4877 0.0140796 -0.127862 -0.01460146 0.5335726 0.6048043 0.04223879