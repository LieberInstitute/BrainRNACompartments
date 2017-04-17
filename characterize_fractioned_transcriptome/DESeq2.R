library('DESeq2')

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

#Parse pd into groups
Cytosol <- pd[which(pd$Zone=="Cytosol"),]
Nucleus <- pd[which(pd$Zone=="Nucleus"),]
Adult <- pd[which(pd$Fetal=="Adult"),]
Fetal <- pd[which(pd$Fetal=="Prenatal"),]
PolyA <- pd[which(pd$Library=="polyA"),]
Ribozero <- pd[which(pd$Library=="RiboZero"),]

Adult.Ribo <-pd[which(pd$Fetal=="Adult" & pd$Library=="RiboZero"),]
Fetal.Ribo <-pd[which(pd$Fetal=="Prenatal" & pd$Library=="RiboZero"),]
Adult.polyA <-pd[which(pd$Fetal=="Adult" & pd$Library=="polyA"),]
Fetal.polyA <-pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),]
Cyt.Ribo <-  pd[which(pd$Zone=="Cytosol" & pd$Library=="RiboZero"),]
Nuc.Ribo <-pd[which(pd$Zone=="Nucleus" & pd$Library=="RiboZero"),]
Nuc.polyA <-pd[which(pd$Zone=="Nucleus"  & pd$Library=="polyA"),]
Cyt.polyA <-pd[which(pd$Zone=="Cytosol"  & pd$Library=="polyA"),]
Adult.cyt <-pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Adult"),]
Adult.nuc <-pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Adult"),]
Fetal.cyt <-pd[which(pd$Zone=="Cytosol" & pd$Fetal=="Prenatal"),]
Fetal.nuc <-pd[which(pd$Zone=="Nucleus" & pd$Fetal=="Prenatal"),]

# Parse count matrices
polya.counts <- geneCounts[,which(colnames(geneCounts)%in%PolyA$SampleID)]
ribozero.counts <- geneCounts[,which(colnames(geneCounts)%in%Ribozero$SampleID)]
Adult.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult$SampleID)]
Fetal.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal$SampleID)]
Cytosol.counts <- geneCounts[,which(colnames(geneCounts)%in%Cytosol$SampleID)]
Nucleus.counts <- geneCounts[,which(colnames(geneCounts)%in%Nucleus$SampleID)]

Adult.Ribo.counts <-geneCounts[,which(colnames(geneCounts)%in%Adult.Ribo$SampleID)]
Fetal.Ribo.counts <-geneCounts[,which(colnames(geneCounts)%in%Fetal.Ribo$SampleID)]
Adult.polyA.counts<-geneCounts[,which(colnames(geneCounts)%in%Adult.polyA$SampleID)]
Fetal.polyA.counts<-geneCounts[,which(colnames(geneCounts)%in%Fetal.polyA$SampleID)]
Cyt.Ribo.counts <-  geneCounts[,which(colnames(geneCounts)%in%Cyt.Ribo$SampleID)]
Nuc.Ribo.counts <-  geneCounts[,which(colnames(geneCounts)%in%Nuc.Ribo$SampleID)]
Nuc.polyA.counts <- geneCounts[,which(colnames(geneCounts)%in%Nuc.polyA$SampleID)]
Cyt.polyA.counts <- geneCounts[,which(colnames(geneCounts)%in%Cyt.polyA$SampleID)]
Adult.cyt.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult.cyt$SampleID)]
Adult.nuc.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult.nuc$SampleID)]
Fetal.cyt.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal.cyt$SampleID)]
Fetal.nuc.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal.nuc$SampleID)]

polya.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%PolyA$SampleID)]
ribozero.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Ribozero$SampleID)]
Adult.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Adult$SampleID)]
Fetal.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Fetal$SampleID)]
Cytosol.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Cytosol$SampleID)]
Nucleus.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Nucleus$SampleID)]

Adult.Ribo.counts.down <-geneCounts.down[,which(colnames(geneCounts.down)%in%Adult.Ribo$SampleID)]
Fetal.Ribo.counts.down <-geneCounts.down[,which(colnames(geneCounts.down)%in%Fetal.Ribo$SampleID)]
Adult.polyA.counts.down<-geneCounts.down[,which(colnames(geneCounts.down)%in%Adult.polyA$SampleID)]
Fetal.polyA.counts.down<-geneCounts.down[,which(colnames(geneCounts.down)%in%Fetal.polyA$SampleID)]
Cyt.Ribo.counts.down <-  geneCounts.down[,which(colnames(geneCounts.down)%in%Cyt.Ribo$SampleID)]
Nuc.Ribo.counts.down <-  geneCounts.down[,which(colnames(geneCounts.down)%in%Nuc.Ribo$SampleID)]
Nuc.polyA.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Nuc.polyA$SampleID)]
Cyt.polyA.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Cyt.polyA$SampleID)]
Adult.cyt.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Adult.cyt$SampleID)]
Adult.nuc.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Adult.nuc$SampleID)]
Fetal.cyt.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Fetal.cyt$SampleID)]
Fetal.nuc.counts.down <- geneCounts.down[,which(colnames(geneCounts.down)%in%Fetal.nuc$SampleID)]

## Using counts from original BAMs
#How many genes are differentially expressed by fraction?
dds <- DESeqDataSetFromMatrix(countData = geneCounts, colData = pd, design = ~ Library + Zone + Fetal)
dds <- DESeq(dds)
res <- results(dds)

Irdds <- DESeqDataSetFromMatrix(countData = ribozero.counts, colData = Ribozero, design = ~ Fetal + Zone + Fetal:Zone)
Irdds <- DESeq(Irdds)
Irres <- results(Irdds)
Ipdds <- DESeqDataSetFromMatrix(countData = polya.counts, colData = PolyA, design = ~ Fetal + Zone + Fetal:Zone)
Ipdds <- DESeq(Ipdds)
Ipres <- results(Ipdds)

Zrdds <- DESeqDataSetFromMatrix(countData = ribozero.counts, colData = Ribozero, design = ~ Fetal + Zone)
Zrdds <- DESeq(Zrdds)
Zrres <- results(Zrdds)
Zpdds <- DESeqDataSetFromMatrix(countData = polya.counts, colData = PolyA, design = ~ Fetal + Zone)
Zpdds <- DESeq(Zpdds)
Zpres <- results(Zpdds)

Ap.dds <- DESeqDataSetFromMatrix(countData = Adult.polyA.counts, colData = Adult.polyA, design = ~ Zone)
Ap.dds <- DESeq(Ap.dds)
Apres <- results(Ap.dds)
Fp.dds <- DESeqDataSetFromMatrix(countData = Fetal.polyA.counts, colData = Fetal.polyA, design = ~ Zone)
Fp.dds <- DESeq(Fp.dds)
Fpres <- results(Fp.dds)
Ar.dds <- DESeqDataSetFromMatrix(countData = Adult.Ribo.counts, colData = Adult.Ribo, design = ~ Zone)
Ar.dds <- DESeq(Ar.dds)
Arres <- results(Ar.dds)
Fr.dds <- DESeqDataSetFromMatrix(countData = Fetal.Ribo.counts, colData = Fetal.Ribo, design = ~ Zone)
Fr.dds <- DESeq(Fr.dds)
Frres <- results(Fr.dds)

#How many genes are differentially expressed by Age?
Ager.dds <- DESeqDataSetFromMatrix(countData = ribozero.counts, colData = Ribozero, design = ~ Zone + Fetal)
Ager.dds <- DESeq(Ager.dds)
Agerres <- results(Ager.dds)
resultsNames(Ager.dds)
Agep.dds <- DESeqDataSetFromMatrix(countData = polya.counts, colData = PolyA, design = ~ Zone + Fetal)
Agep.dds <- DESeq(Agep.dds)
Agepres <- results(Agep.dds)

Cp.dds <- DESeqDataSetFromMatrix(countData = Cyt.polyA.counts, colData = Cyt.polyA, design = ~ Fetal)
Cp.dds <- DESeq(Cp.dds)
Cpres <- results(Cp.dds)
Np.dds <- DESeqDataSetFromMatrix(countData = Nuc.polyA.counts, colData = Nuc.polyA, design = ~ Fetal)
Np.dds <- DESeq(Np.dds)
Npres <- results(Np.dds)
Cr.dds <- DESeqDataSetFromMatrix(countData = Cyt.Ribo.counts, colData = Cyt.Ribo, design = ~ Fetal)
Cr.dds <- DESeq(Cr.dds)
Crres <- results(Cr.dds)
Nr.dds <- DESeqDataSetFromMatrix(countData = Nuc.Ribo.counts, colData = Nuc.Ribo, design = ~ Fetal)
Nr.dds <- DESeq(Nr.dds)
Nrres <- results(Nr.dds)

#Gene expression by library
Ln.dds <- DESeqDataSetFromMatrix(countData = Nucleus.counts, colData = Nucleus, design = ~ Fetal + Library)
Ln.dds <- DESeq(Ln.dds)
Lnres <- results(Ln.dds)

# Combining both libraries
Idds = DESeqDataSetFromMatrix(countData = geneCounts, colData = pd, design = ~ Library + Fetal + Zone + Fetal:Zone)
Idds = DESeq(Idds)
Ires = results(Idds)

Adds = DESeqDataSetFromMatrix(countData = Adult.counts, colData = Adult, design = ~ Library + Zone)
Adds = DESeq(Adds)
Ares = results(Adds)

Fdds = DESeqDataSetFromMatrix(countData = Fetal.counts, colData = Fetal, design = ~ Library + Zone)
Fdds = DESeq(Fdds)
Fres = results(Fdds)

Cdds = DESeqDataSetFromMatrix(countData = Cytosol.counts, colData = Cytosol, design = ~ Library + Fetal)
Cdds = DESeq(Cdds)
Cres = results(Cdds)

Ndds = DESeqDataSetFromMatrix(countData = Nucleus.counts, colData = Nucleus, design = ~ Library + Fetal)
Ndds = DESeq(Ndds)
Nres = results(Ndds)

Fdds.down = DESeqDataSetFromMatrix(countData = Fetal.counts.down, colData = Fetal, design = ~ Library + Zone)
Fdds.down = DESeq(Fdds.down)
Fres.down = results(Fdds.down)

Cdds.down = DESeqDataSetFromMatrix(countData = Cytosol.counts.down, colData = Cytosol, design = ~ Library + Fetal)
Cdds.down = DESeq(Cdds.down)
Cres.down = results(Cdds.down)


## Using counts from downsampled BAMs
#How many genes are differentially expressed by fraction?
dds.down <- DESeqDataSetFromMatrix(countData = geneCounts.down, colData = pd, design = ~ Library + Zone + Fetal)
dds.down <- DESeq(dds.down)
res.down <- results(dds.down)

Ipdds.down <- DESeqDataSetFromMatrix(countData = polya.counts.down, colData = PolyA, design = ~ Fetal + Zone + Fetal:Zone)
Ipdds.down <- DESeq(Ipdds.down)
Ipres.down <- results(Ipdds.down)

Zpdds.down <- DESeqDataSetFromMatrix(countData = polya.counts.down, colData = PolyA, design = ~ Fetal + Zone)
Zpdds.down <- DESeq(Zpdds.down)
Zpres.down <- results(Zpdds.down)

Fp.dds.down <- DESeqDataSetFromMatrix(countData = Fetal.polyA.counts.down, colData = Fetal.polyA, design = ~ Zone)
Fp.dds.down <- DESeq(Fp.dds.down)
Fpres.down <- results(Fp.dds.down)

#How many genes are differentially expressed by Age?
Agep.dds.down <- DESeqDataSetFromMatrix(countData = polya.counts.down, colData = PolyA, design = ~ Zone + Fetal)
Agep.dds.down <- DESeq(Agep.dds.down)
Agepres.down <- results(Agep.dds.down)

Cp.dds.down <- DESeqDataSetFromMatrix(countData = Cyt.polyA.counts.down, colData = Cyt.polyA, design = ~ Fetal)
Cp.dds.down <- DESeq(Cp.dds.down)
Cpres.down <- results(Cp.dds.down)

#Gene expression by library
Ln.dds.down <- DESeqDataSetFromMatrix(countData = Nucleus.counts.down, colData = Nucleus, design = ~ Fetal + Library)
Ln.dds.down <- DESeq(Ln.dds.down)
Lnres.down <- results(Ln.dds.down)


save(res.down,Ipres.down,Zpres.down,Fpres.down,Agepres.down,Cpres.down,Lnres.down,
     res,Irres,Ipres,Zrres,Zpres,Apres,Fpres,Arres,Frres,
     Agerres,Agepres,Cpres,Npres,Crres,Nrres,Lnres,
     Ires,Ares,Fres,Cres,Nres,Fres.down,Cres.down, geneMap,
     file="./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")