# read in Dougherty table
dopd = read.table("./Dropbox/sorted_figures/new/dougherty.pd.txt", header=T)
doCounts = read.csv("./Dropbox/sorted_figures/new/GSE73391_Counts.csv", header=T)
rownames(doCounts) = doCounts$X
doCounts = doCounts[,-1]
colnames(doCounts) = gsub("X920_", "", colnames(doCounts))
pd = data.frame(ID = colnames(doCounts), Fraction = c(rep.int("Nucleus", 24), "Cytosol","Cytosol"),
                Cell.Type = c(rep.int("Neuron",4),rep.int("Astrocyte",4),rep.int("Oligodendrocyte",4),
                              rep.int("Neuron",4),rep.int("Astrocyte",4),rep.int("Oligodendrocyte",4),"Presort","Presort"),
                Sex = c(rep.int("M",12), rep.int("F",12), "NA","NA"), Sorted = gsub("\\_.*", "", colnames(doCounts)))
pd = pd[order(pd$Sorted),]
pd[1:14,"Cell.Type"] = "Presort"
pd = pd[,1:4]
pd = pd[order(pd$ID),]
doCounts = doCounts[,order(colnames(doCounts))]
doCounts = as.matrix(doCounts)
save(doCounts, pd, file="./Dropbox/sorted_figures/new/Dougherty_counts.rda")
