library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


### For RNA Binding protein analysis, create DNA strings for PWM analysis

editing_anno$start = editing_anno$end
keep = editing_anno[overlappingGene!="NA" & collapsedconversion=="A:G / T:C",,]
keep_gr = makeGRangesFromDataFrame(keep)
keep_gr = unique(keep_gr)
keep_df = data.frame(keep_gr)
keep_df$start = keep_df$start - 10
keep_df$end = keep_df$end + 10
keep_gr = makeGRangesFromDataFrame(keep_df)
seq = getSeq(Hsapiens, keep_gr)
seq = as.list(as.character(seq, use.names=TRUE))
revcomp = unlist(lapply(seq, function(x) c2s(rev(comp(s2c(x), forceToLower = F)))))
seq = unlist(seq)
write.table(paste0(">", keep_df$seqnames, ":",keep_df$start,"-",keep_df$end,":", keep_df$strand,
                   "\n", revcomp), col.names = F, quote = F,row.names = F,
            file = "./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/editing_site_sequence_revcomp.txt")
write.table(paste0(">", keep_df$seqnames, ":",keep_df$start,"-",keep_df$end,":", keep_df$strand,
                   "\n", seq), col.names = F, quote = F,row.names = F,
            file = "./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/editing_site_sequence.txt")

names(revcomp) = names(seq) = paste0(keep_df$seqnames, ":",keep_df$start,"-",keep_df$end,":", keep_df$strand)

sep = seprc = list()
mut = mutrc = vector(mode = "character", length = length(seq))
for (i in 1:length(seq)) {
  sep[[i]] = s2c(seq[[i]])
  seprc[[i]] = s2c(revcomp[i])
  mutrc[i] = c2s(c(seprc[[i]][1:10], "C", seprc[[i]][12:21]))
  mut[i] = c2s(c(sep[[i]][1:10], "G", sep[[i]][12:21]))
}
names(mut) = names(mutrc) = names(seq)
sequence = list(refseq = DNAStringSet(seq), refseqrc = DNAStringSet(revcomp),
                altseq = DNAStringSet(mut), altseqrc = DNAStringSet(mutrc))


## load RBP database from ATtRACT

rbpdb = read.table("./Dropbox/BrainRNACompartments/ATtRACT/ATtRACT_db.txt", header=T, sep = "\t")
rbpdb = rbpdb[rbpdb$Organism=="Homo_sapiens",]
ids = as.character(unique(rbpdb$Gene_name))


## load motif databases

pwm = read.table("./Dropbox/BrainRNACompartments/ATtRACT/pwm.txt", header = FALSE, 
                 sep = "\t", col.names = paste0("V",seq_len(4)), fill = TRUE)
colnames(pwm) = c("A","C","G","T")
ast = grep(">", pwm$A, fixed = T)
sep = list()
for (i in 1:length(ast)) { sep[[i]] = pwm[c(ast[i]:(ast[i+1]-1)),] }
sep = c(sep, list(pwm[12798:12802,]))
names(sep) = lapply(sep, function(x) gsub(">","",as.character(x[1,1])))
sep = lapply(sep, function(x) x[2:nrow(x),])
for (i in 1:length(sep)) {
  tmp = sep[[i]]
  tmp$A = as.numeric(as.character(tmp$A))
  tmp$C = as.numeric(as.character(tmp$C))
  tmp$G = as.numeric(as.character(tmp$G))
  tmp$'T' = as.numeric(as.character(tmp$'T'))
  sep[[i]] = t(tmp)
  colnames(sep[[i]]) = NULL
}
sep = sep[names(sep) %in% rbpdb$Matrix_id]
head(sep)
#pos = lapply(sep, function(x) round(x*10000))
#pos = lapply(pos, function(m) apply (m, c (1, 2), function (x) { (as.integer(x)) }))
#bg = makeBackground(pos, organism = "hg19", type = "logn", algorithm = "human",bg.len=50)
pwms = toPWM(sep)


## Calculate enrichment for all sequences together

res_all = lapply(sequence, function(x) motifEnrichment(x, pwms))
save(res_all, file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/ATtRACT_res_all.rda")

res_ind = list()
for (i in 1:length(sequence[[1]])) {
  res_ind[[i]] = lapply(sequence, function(x) motifEnrichment(x[i], pwms))
}
save(res_ind, compress = "gzip",
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/ATtRACT_res_ind.rda")


report = lapply(res_all, function(x) groupReport(x), by.top.motifs=TRUE)
lapply(report, head)
plot(report[1:10], fontsize=7, id.fontsize=5)


## for a single site


report = sequenceReport(res, 1)
report
# plot the motif with P-value < 0.05
plot(report[report$p.value < 0.05], fontsize=7, id.fontsize=6)










txdb = loadDb("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
transcripts = transcripts(txdb)
ex = exonsBy(txdb, use.names=T)
ov = findOverlaps(transcripts, keep_gr)
tx = cbind(data.frame(transcripts[queryHits(ov)]), data.frame(keep_gr[subjectHits(ov)]))
ex = ex[names(ex) %in% tx$tx_name]
seqs = extractTranscriptSeqs(Hsapiens, ex)

sernacoord = paste0(keep$seqnames, ":", keep$start-10,"-", keep$end+10, ":", keep$strand)
rnacoord = unique(rnacoord)
length(rnacoord) # 18907 editing sites
