#***********************************************************************************************************************
# Demo example for the use of the two function
# Get_Summary_Statistics_for_Comorbidity() and Get_Summary_Statistics_for_Only_Single_Trait()
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
fwrite(comor_result,"Toy_example_Comorbidity_results.txt",sep="\t")
proc.time()-t1
# user  system elapsed 
# 18.036   0.400  19.722 

t2=proc.time()
single_result = Get_Summary_Statistics_for_Only_Single_Trait(sumstats,Size_TraitA,Size_TraitB,TraitA_preval,TraitB_preval,Comorbid_preval,LDSR_intercept)
proc.time()-t2  
# user  system elapsed 
# 25.606   0.463  27.336 
fwrite(comor_result,"Toy_example_Only_Single_Trait_results.txt",sep="\t")
