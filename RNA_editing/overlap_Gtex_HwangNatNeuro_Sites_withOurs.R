library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


## load public editing site lists
gtex_all = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/GTEX_Supp1_allSites.xlsx')
gtex_all = cbind(gtex_all[,grep("Brain",colnames(gtex_all))], 
          seqnames = unlist(strsplit(gtex_all$sites,"_",fixed=T))[grep("chr",unlist(strsplit(gtex_all$sites,"_",fixed=T)))],
          start = unlist(strsplit(gtex_all$sites,"_",fixed=T))[-grep("chr",unlist(strsplit(gtex_all$sites,"_",fixed=T)))], site = gtex_all$sites)
gtex_all = makeGRangesFromDataFrame(gtex_all, end.field = "start", keep.extra.columns = T)

gtex_ts = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/GTEX_Supp3_tissueSpecificSites.xlsx')
gtex_ts = makeGRangesFromDataFrame(gtex_ts[,1:4], seqnames.field="Chromosome",start.field="Position", end.field="Position",keep.extra.columns = T)

increasing = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/Hwang_NatNeuro_SuppTables.xlsx',sheet =3)
increasing = increasing[-1,]
incr = makeGRangesFromDataFrame(increasing, seqnames.field = "Chromosome",start.field = "Coordinate",end.field = "Coordinate", strand.field = "Strand")

editing_annogr = editing_anno[collapsedconversion=="A:G / T:C",,]
editing_annogr$start = editing_annogr$end
editing_annogr = makeGRangesFromDataFrame(editing_annogr)
reduced = reduce(editing_annogr)


## How many of our sites have been identified already?
red_all = findOverlaps(reduced, reduce(gtex_all))
length(unique(queryHits(red_all))) # 12245 A to I sites are identified in the Gtex samples out of 17876 
red_ts = findOverlaps(reduced, reduce(gtex_ts))
length(unique(queryHits(red_ts))) # 66 are identified as tissue-specific in GTEX 
ov_inc = findOverlaps(reduced, incr)
length(unique(queryHits(ov_inc))) # 554 of the 742 increasing editing sites identified in Taeyoung's paper are in our set


## Annotate our sites according to the GTEX data
ov_all = findOverlaps(editing_annogr, gtex_all)
ov_ts = findOverlaps(editing_annogr, gtex_ts)
ov_inc = findOverlaps(editing_annogr, incr)

editing_anno$rnum = 1:nrow(editing_anno)
editing_anno$new = ifelse(editing_anno$rnum %in% queryHits(ov_all), "In GTEX","Novel")
editing_anno$tissuespecific = ifelse(editing_anno$rnum %in% queryHits(ov_ts), "Tissue-Specific","Multi-Tissue")
editing_anno$Increasing = ifelse(editing_anno$rnum %in% queryHits(ov_inc), "Increasing","Not Increasing")


## Plot the distribution

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/GTEX_HwangNatNeuro_editingSites_inOurData.pdf",width=8,height=5)
x = editing_anno[,length(unique(editingID)),by=c("Fraction","Age","new")]
x$perc = c(round(x$V1[1:2]/sum(x$V1[1:2])*100,1),round(x$V1[3:4]/sum(x$V1[3:4])*100,1),round(x$V1[5:6]/sum(x$V1[5:6])*100,1),round(x$V1[7:8]/sum(x$V1[7:8])*100,1))
ggplot(x, aes(x = Fraction, y = V1, fill = new)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  facet_grid(. ~ Age) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to GTEX") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = editing_anno[,length(unique(editingID)),by=c("Fraction","Age","tissuespecific")]
x$perc = c(round(x$V1[1:2]/sum(x$V1[1:2])*100,1),round(x$V1[3:4]/sum(x$V1[3:4])*100,1),round(x$V1[5:6]/sum(x$V1[5:6])*100,1),round(x$V1[7:8]/sum(x$V1[7:8])*100,1))
ggplot(x[tissuespecific=="Tissue-Specific",,], aes(x = Age, y = V1, fill = Fraction)) + 
  geom_bar(stat = "identity",position=position_dodge(width=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Our Editing Sites Compared to GTEX") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = editing_anno[,length(unique(editingID)),by=c("Fraction","Age","Increasing")]
x$perc = c(round(x$V1[1:2]/sum(x$V1[1:2])*100,1),round(x$V1[3:4]/sum(x$V1[3:4])*100,1),round(x$V1[5:6]/sum(x$V1[5:6])*100,1),round(x$V1[7:8]/sum(x$V1[7:8])*100,1))
ggplot(x[Increasing=="Increasing",,], aes(x = Age, y = V1, fill = Fraction)) + geom_bar(stat = "identity",position=position_dodge(width=1)) +
  geom_text(aes(label = paste0(perc,"%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
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
unique_all = do.call(rbind,unique_all)
unique_all$name = factor(unique_all$name, levels= c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/GTEX_HwangNatNeuro_editingSites_inOurData_uniqueAll3.pdf",width=18,height=8)
x = unique_all[,length(unique(editingID)),by=c("name","new")]
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
x = unique_all[,length(unique(editingID)),by=c("name","tissuespecific")]
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
x = unique_all[,length(unique(editingID)),by=c("name","Increasing")]
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

al = as.data.frame(editing_anno[,length(unique(editingID)),by="Increasing"])
un = as.data.frame(unique_all[,length(unique(editingID)),by=c("name","Increasing")])
names = unique(unique_all$name)
inc = list()
for (i in 1:length(names)){
  inc[[i]] = data.frame(inanno = c(un[un$name==names[i] & un$Increasing=="Increasing","V1"], 
                                   un[un$name==names[i] & un$Increasing=="Not Increasing","V1"]),
                        notanno = c(al[al$Increasing=="Increasing","V1"]-un[un$name==names[i] & un$Increasing=="Increasing","V1"], 
                                    al[al$Increasing=="Not Increasing","V1"]-un[un$name==names[i] & un$Increasing=="Not Increasing","V1"]),
                        row.names = c("Increasing","Not"))
}
names(inc) = names
fisher = lapply(inc, fisher.test)
write.csv(data.frame(inc),quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/HwangNatNeuro_editingSites_contingencytables_inOurData.csv")
write.csv(rbind(pval = data.frame(lapply(fisher, function(x) x$p.value)),OR = data.frame(lapply(fisher, function(x) round(x$estimate,3)))),
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/HwangNatNeuro_editingSites_enrichment_inOurData.csv",
          quote = F)


## Are the increasing sites more edited in adult than prenatal?

ttests = list(all = t.test(x = editing_anno[Increasing=="Increasing" & Age=="Adult",list(rate),], 
                           y = editing_anno[Increasing=="Increasing" & Age=="Prenatal",list(rate),]),
              nucleus = t.test(x = editing_anno[Increasing=="Increasing" & Group=="Adult:Nucleus",list(rate),], 
                               y = editing_anno[Increasing=="Increasing" & Group=="Prenatal:Nucleus",list(rate),]),
              cytosol = t.test(x = editing_anno[Increasing=="Increasing" & Group=="Adult:Cytosol",list(rate),], 
                               y = editing_anno[Increasing=="Increasing" & Group=="Prenatal:Cytosol",list(rate),]))

write.csv(rbind(Tstat = data.frame(lapply(ttests, function(x) x$statistic)),pval = data.frame(lapply(ttests, function(x) x$p.value)),
                confInt = data.frame(lapply(ttests, function(x) x$conf.int)),estMeans = data.frame(lapply(ttests, function(x) x$estimate))),
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/HwangNatNeuro_editingRate_ttest_byAge_inIncreasingSites.csv",
          quote=F)