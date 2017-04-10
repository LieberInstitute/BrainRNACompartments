# Is IR increased in the nucleus?
# Main effect: ALL filtered introns
Measured in polyA, the nucleus has more IR (p-value < 2.2e-16)
#  mean of x  mean of y 
#0.01578728 0.01004347

# in adult polyA
#t = 51.234, df = 187730, p-value < 2.2e-16
#  mean of x   mean of y 
#0.018589335 0.007632962 

# in fetal polyA
#t = 14.141, df = 350770, p-value < 2.2e-16
#  mean of x  mean of y 
#0.01389595 0.01142055


# Main effect: introns >=75% retained
Is a greater percent of introns retained at 75% or more in the nucleus compared to cytosol? No.
# in all polyA
#t = 0.11417, df = 7.6003, p-value = 0.456
#  mean of x mean of y 
#0.1046333 0.1009900

How about 50%? Almost.
#t = 1.7437, df = 9.9664, p-value = 0.05596
#  mean of x mean of y 
#0.3478683 0.2138433

How about in adult polyA, at 75%? Yes.
#t = 2.644, df = 3.1771, p-value = 0.03644
#  mean of x  mean of y 
#0.12810667 0.05880333

in fetal polyA, at 75%? No.
#t = -1.3617, df = 2.0908, p-value = 0.8493
#  mean of x mean of y 
#0.0811600 0.1431767
in fetal polyA, at 50%? No.
#t = -0.40062, df = 2.7842, p-value = 0.6413
#  mean of x mean of y 
#0.2604233 0.2967300 

## Is IR increased in adults?
# Main effect: All filtered introns
Measured in polyA, the Adult has more IR than fetal
#t = 5.6143, df = 464860, p-value = 9.874e-09
#  mean of x  mean of y 
#0.01341493 0.01262374

# in Nuclear polyA
Measured in nuclear polyA, the Adult has more IR than fetal
#t = 21.301, df = 214160, p-value < 2.2e-16
#  mean of x  mean of y 
#0.01858933 0.01389595 

# in Cytosol polyA
Measured in cytosol polyA, the Adult does not have more IR than fetal
#t = -22.705, df = 275350, p-value = 1
#  mean of x   mean of y 
#0.007632962 0.011420549


# Main effect: introns >=75% retained
# in all polyA
Adult does not have more 75% retained introns than fetal
#t = -0.59634, df = 9.4844, p-value = 0.7175
#  mean of x mean of y 
#0.0934550 0.1121683
Adult does not have more 50% retained introns than fetal
#t = 0.05194, df = 7.6364, p-value = 0.48
#  mean of x mean of y 
#0.2831350 0.2785767

# in Nuclear polyA
Adult almost has more 75% retained introns than fetal in the nucleus
#t = 1.9761, df = 2.3527, p-value = 0.08367
#  mean of x mean of y 
#0.1281067 0.0811600
Adult almost has more 50% retained introns than fetal in the nucleus
#t = 1.9576, df = 2.8095, p-value = 0.07569
#  mean of x mean of y 
#0.4353133 0.2604233

# in cytosol polyA
Adult does not have more 75% retained introns than fetal in the cytosol
#t = -1.8001, df = 2.3304, p-value = 0.9023
#  mean of x  mean of y 
#0.05880333 0.14317667
Adult does not have more 50% retained introns than fetal in the cytosol
#t = -1.9755, df = 2.1297, p-value = 0.9105
#  mean of x mean of y 
#0.1309567 0.2967300

## Differences by intron
These are the genes/introns whose IR was compared by the variables below, that pass filtering step.
Most differentially retained introns are retained in the nucleus, and are preferentially 
in fetal cytosol compared to adult cytosol but adult nucleus compared to fetal nucleus.
Few genes with changing IR by the measured variable had a mean IR ratio >=50, but those that 
did were almost all preferentially in the nucleus.
#alternative hypothesis: true odds ratio is not equal to 1

#data:  polyA adult
#p-value = 1

#data:  polyA fetal
#p-value = 2.19e-05

#data:  polyA Cytosol
#p-value = 0.001621

#data:  polyA nucleus
#p-value = 0.002444

elementLengths(lapply(IRclean, function(x) x[which(x$A.IRratio>=0.5 | x$B.IRratio>=0.5),]))
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 3                 4                10                 2                 1
elementLengths(lapply(lapply(IRclean, function(x) x[which(x$A.IRratio>=0.5 | x$B.IRratio>=0.5),]), function(x) 
  x[which(x$Sign=="MoreIRInNuc.Fetal"),]))
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 2                 4                 3                 2                 1

# Are IR genes significantly regulated by fraction?

These tests were on the intron with the highest IR ratio in genes 
differentially expressed by fraction or age. 
Genes differentially expressed by fraction in adult polyA 
contain introns that are significantly more retained in nuclus than cytosol
#LFC>1: t = -38.132, df = 5636.5, p-value < 2.2e-16
#  mean of x  mean of y 
#0.02157524 0.13782572

#NO LFC Cutoff: t = -50.035, df = 17027, p-value < 2.2e-16
#  mean of x  mean of y 
#0.02447499 0.09868088
t.test(FPC$IRratio, FPN$IRratio, alternative = "less")
#data:  FPC$IRratio and FPN$IRratio
#t = -33.876, df = 2177, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.1452936
#sample estimates:
#  mean of x   mean of y 
#0.004986751 0.157698449 

comparing 2 variables: significantly enriched IR in nucleus, and
significantly differentially expressed by fraction at the gene levela
dult = data.frame(c(69,4348), c(70,12568))
fetal = data.frame(c(69,4348), c(71,12567))
#data:  adult
#p-value = 2.074e-09
#data:  fetal
#p-value = 2.687e-09

in the cytosol, IR ratio is higher in fetal than adult
t.test(CPA$IRratio, CPF$IRratio)
#t = -7.7772, df = 39783, p-value = 7.593e-15
#  mean of x mean of y 
#0.0398155 0.0480012 
in the nucleus, IR ratio is higher in adult than fetal
#t = 6.7818, df = 33380, p-value = 1.206e-11
#  mean of x  mean of y 
#0.04770549 0.03977136

there is no interaction between whether a gene is significantly age-related or not 
and whether IR increases in the adult or fetal
#data:  cyt
#p-value = 0.1922
#alternative hypothesis: true odds ratio is not equal to 1

#data:  nuc
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1

##Working in the opposite direction
compare DEG significance to IR significance
#data:  adult polyA
#p-value < 2.2e-16

#data: fetal polyA
#p-value < 2.2e-16

# genes >50% retained in at least one sample
Compare IR percent >50 and DEG significance
#data:  both_retained
#p-value = 5.193e-10
both_exported = data.frame(c(0,906), c(79,9229))
#data:  both_exported
#p-value = 0.001091
Fet_retained =  data.frame(c(6,139), c(73,9996))
#data:  Fet_retained
#p-value = 0.0008756
Ad_retained =  data.frame(c(28,1641), c(51,8494))
#data:  Ad_retained
#p-value = 3.928e-05
ret_Ad_exp_Fet = data.frame(c(0,11), c(79,10124))
#data:  ret_Ad_exp_Fet
#p-value = 1
ret_Fet_exp_Ad = data.frame(c(0,4), c(79,10131))
#data:  ret_Fet_exp_Ad
#p-value = 1
interacting = data.frame(c(13,916), c(66,9219))
#data:  interacting
#p-value = 0.02981

# Intron retention between retained and exported genes
#data:  ret$IRratio and exp$IRratio
#t = 13.297, df = 521.58, p-value < 2.2e-16
#  mean of x mean of y 
#0.1166033 0.0271476

# plot gene that changes LFC by fraction and IR ratio
# plot gene that changes age:fraction that also changes IR ratio
