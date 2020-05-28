#***********************************************************************************************************************
# Demo example for the use of the functions
# R version: 3.6.3
# Operating system: MacOS Catalina(10.15.4)
# Processor:2.7 GHz Dual-Core Intel Core i5
# Memory:8 GB 1867 MHz DDR3
#***********************************************************************************************************************


#***********************************************************************************************************************
# Demo code for Get_Summary_Statistics_for_Comorbidity() and Get_Summary_Statistics_for_Only_Single_Trait
#***********************************************************************************************************************
library(data.table)

sumstats = fread("Harmonized_data_Toy_example.txt")

# The first few lines of the sumstats look like this: 
# SNP effect_allele.exposure other_allele.exposure effect_allele.outcome
# 1:       rs10                      A                     C                     A
# 2:  rs1000000                      A                     G                     A
# 3: rs10000010                      C                     T                     C
# 4: rs10000012                      C                     G                     C
# 5: rs10000013                      C                     A                     C
# 6: rs10000017                      T                     C                     T
# other_allele.outcome beta.exposure beta.outcome eaf.exposure eaf.outcome remove palindromic
# 1:                    C       0.04276       0.0470      0.07591       0.033  FALSE       FALSE
# 2:                    G       0.00319      -0.0130      0.21756       0.373  FALSE       FALSE
# 3:                    T       0.00122       0.0023      0.49824       0.425  FALSE       FALSE
# 4:                    G       0.00253       0.0220      0.87103       0.808  FALSE        TRUE
# 5:                    A       0.00758       0.0340      0.22350       0.167  FALSE       FALSE
# 6:                    C       0.00506      -0.0120      0.20607       0.223  FALSE       FALSE
# ambiguous id.outcome se.outcome outcome mr_keep.outcome pval.outcome pval_origin.outcome
# 1:     FALSE     HdthAE      0.044 Obesity            TRUE    0.1427193            inferred
# 2:     FALSE     HdthAE      0.015 Obesity            TRUE    0.1930623            inferred
# 3:     FALSE     HdthAE      0.012 Obesity            TRUE    0.4240017            inferred
# 4:     FALSE     HdthAE      0.018 Obesity            TRUE    0.1108118            inferred
# 5:     FALSE     HdthAE      0.015 Obesity            TRUE    0.0117053            inferred
# 6:     FALSE     HdthAE      0.015 Obesity            TRUE    0.2118554            inferred
# se.exposure pval.exposure exposure mr_keep.exposure pval_origin.exposure id.exposure action
# 1:     0.02052      0.037195      CAD             TRUE             reported      WouUDZ      2
# 2:     0.00979      0.744846      CAD             TRUE             reported      WouUDZ      2
# 3:     0.00815      0.881128      CAD             TRUE             reported      WouUDZ      2
# 4:     0.01249      0.839695      CAD             TRUE             reported      WouUDZ      2
# 5:     0.00992      0.444560      CAD             TRUE             reported      WouUDZ      2
# 6:     0.01032      0.624046      CAD             TRUE             reported      WouUDZ      2
# mr_keep samplesize.outcome
# 1:    TRUE                 NA
# 2:    TRUE                 NA
# 3:    TRUE                 NA
# 4:    TRUE                 NA
# 5:    TRUE                 NA
# 6:    TRUE                 NA

Size_TraitA = 148715
Size_TraitB = 98697
TraitA_preval = 0.405
TraitB_preval = 0.3
Comorbid_preval = 0.1422289
LDSR_intercept = 0.0172

Size_TraitA = 148715
Size_TraitB = 98697
TraitA_preval = 0.405
TraitB_preval = 0.3
Comorbid_preval = 0.1422289
LDSR_intercept = 0.0172

t1=proc.time()
comor_result = Get_Summary_Statistics_for_Comorbidity(sumstats,Size_TraitA,Size_TraitB,TraitA_preval,TraitB_preval,Comorbid_preval,LDSR_intercept)
# The first few lines of the comor_result look like this: 
# SNP A1 A2 trait1_A1F trait2_A1F trait1_beta trait1_se trait1_pval trait2_beta
# rs10019946  G  A    0.84500      0.858     0.04318   0.01207    0.000348       0.054
# rs1000722  G  A    0.69613      0.650     0.01677   0.00892    0.060070       0.052
# rs10017412  G  A    0.50589      0.517     0.01081   0.00810    0.182166       0.048
# rs10015696  C  A    0.77792      0.718     0.02668   0.00975    0.006229       0.046
# rs10022347  G  A    0.53411      0.551     0.01035   0.00817    0.205215       0.046
# rs10009336  C  T    0.82956      0.825     0.02808   0.01087    0.009775       0.052
# trait2_se  trait2_pval   avg_eaf Comorb_beta   Comorb_se  Comorb_pval
# 0.018 1.349898e-03 0.8508361  0.06224200 0.014120332 1.043440e-05
# 0.014 1.018892e-04 0.6754208  0.04590803 0.011043834 3.226060e-05
# 0.012 3.167124e-05 0.5108776  0.03971895 0.009631789 3.727773e-05
# 0.015 1.082300e-03 0.7510200  0.04742704 0.011830525 6.100987e-05
# 0.012 6.320923e-05 0.5416925  0.03805698 0.009650545 8.029916e-05
# 0.017 1.111040e-03 0.8275129  0.05229117 0.013305230 8.490289e-05
fwrite(comor_result,"Toy_example_Comorbidity_results.txt",sep="\t")
proc.time()-t1
# user  system elapsed 
# 18.036   0.400  19.722 

t2=proc.time()
single_result = Get_Summary_Statistics_for_Only_Single_Trait(sumstats,Size_TraitA,Size_TraitB,TraitA_preval,TraitB_preval,Comorbid_preval,LDSR_intercept)
# The first few lines of the single_result look like this: 
# SNP A1 A2 trait1_A1F trait1_beta trait1_se trait1_pval trait2_A1F trait2_beta
# rs10013305  G  A    0.71713     0.02723   0.00910    0.002787      0.776      -0.042
# rs1000940  A  G    0.67750     0.01265   0.00872    0.146685      0.775      -0.058
# rs10015965  A  G    0.64499     0.02966   0.00868    0.000636      0.724      -0.024
# rs10011097  T  G    0.71706     0.02764   0.00960    0.004013      0.783      -0.036
# rs10018280  T  G    0.71448     0.02601   0.00909    0.004213      0.775      -0.035
# rs10014719  G  A    0.71906     0.02699   0.00977    0.005734      0.783      -0.036
# trait2_se  trait2_pval   avg_eaf        single_beta          single_se
# 0.017 6.744552e-03 0.7435586 0.0671505619511518 0.0172605339535787
# 0.013 4.068667e-06 0.7212708 0.0583041023470487 0.0150250829354872
# 0.013 3.243494e-02 0.6804601 0.0578195577080177 0.0149536179194003
# 0.014 5.063995e-03 0.7466625 0.0635075025138351  0.016445098365256
# 0.014 6.209665e-03 0.7416493 0.0605301085037842 0.0158725122460065
# 0.014 5.063995e-03 0.7477647 0.0626057615264842 0.0166339433964121
# single_pval
# 0.000100074797039374
# 0.000104262828414881
# 0.000110366287248796
# 0.000112559689128026
# 0.000137002634960711
# 0.000167394006908089
proc.time()-t2  
# user  system elapsed 
# 25.606   0.463  27.336 
fwrite(comor_result,"Toy_example_Only_Single_Trait_results.txt",sep="\t")

#***********************************************************************************************************************
# Demo code for transformation of coefficients derived from linear regression to that of logistic regression
#***********************************************************************************************************************
library(dplyr) 
library(data.table)

sumstats = fread("Quantitive_trait_toy_example.txt")
# The first few lines of the sumstats look like this: 
# SNP_hg18       SNP_hg19        SNP A1 A2   BETA     SE         N         P
# chr2:21302112  chr2:21448607   rs428696  C  T 0.1005 0.0057 126558.92 1.972e-65
# chr5:74664131  chr5:74628375  rs3843481  T  A 0.0682 0.0038 162085.00 7.521e-65
# chr5:74639235  chr5:74603479  rs2335418  A  G 0.0654 0.0037 173001.00 8.583e-65
# chr1:109628776 chr1:109827253   rs672569  G  A 0.1431 0.0082  89876.00 2.083e-64
# chr5:74913559  chr5:74877803  rs5744672  C  T 0.0654 0.0037 172857.00 4.528e-64
# chr2:21162627  chr2:21309122 rs28562532  C  T 0.1340 0.0076  83137.01 6.892e-64
# Freq.A1.1000G.EUR
# 0.7916
# 0.4103
# 0.4578
# 0.8430
# 0.4169
# 0.8193

mean_trait = 133.6
sd_trait = 39.0
thres = 190
beta_cont = sumstats$BETA
se_cont = sumstats$SE
eaf = sumstats$Freq.A1.1000G.EUR
rank_normalized_y = TRUE

beta_esti = Binary_GWAS(beta_cont = beta_cont, 
                       mean_trait = mean_trait,  
                       sd_trait = sd_trait, 
                       thres_raw = thres, 
                       eaf = eaf, 
                       rank_normalized_y = rank_normalized_y) 
se_esti = Binary_GWAS_SE( beta_cont = beta_cont, 
                         se_cont = se_cont, 
                         mean_trait = mean_trait,  
                         sd_trait = sd_trait, 
                         thres_raw = thres,
                         eaf = eaf, 
                         rank_normalized_y = rank_normalized_y) 
zscore_esti = beta_esti/se_esti 
pval_esti = pnorm(-abs(zscore_esti), lower.tail=TRUE)
result_binary = cbind(sumstats[,-c(6,7,9)],beta_esti,se_esti,pval_esti)
colnames_result_binary = c("SNP_hg18","SNP_hg19","SNP","A1","A2","N","Freq.A1","BETA","SE","P")
colnames(result_binary) = colnames_result_binary
# The first few lines of the result_binary look like this: 
# SNP_hg18       SNP_hg19        SNP A1 A2         N Freq.A1      BETA          SE
# chr2:21302112  chr2:21448607   rs428696  C  T 126558.92  0.7916 0.2107437 0.012291245
# chr5:74664131  chr5:74628375  rs3843481  T  A 162085.00  0.4103 0.1402076 0.007868969
# chr5:74639235  chr5:74603479  rs2335418  A  G 173001.00  0.4578 0.1346151 0.007679854
# chr1:109628776 chr1:109827253   rs672569  G  A  89876.00  0.8430 0.3044443 0.018230000
# chr5:74913559  chr5:74877803  rs5744672  C  T 172857.00  0.4169 0.1344402 0.007659646
# chr2:21162627  chr2:21309122 rs28562532  C  T  83137.01  0.8193 0.2840442 0.016771617
# P
# 3.376011e-66
# 2.571352e-71
# 4.353517e-69
# 6.533568e-63
# 2.883484e-69
# 1.220762e-64
fwrite(as.data.frame(result_binary),"Binary_summary_statistics_for_Quantitive_trait_toy_example.txt")
