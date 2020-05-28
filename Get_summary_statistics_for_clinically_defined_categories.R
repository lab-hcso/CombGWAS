#**********************************************************************************************
# Function: transform beta derived from linear regression to that from logistic regression
# beta_cont: beta derived from linear regression
# mean_trait: mean value of the studied quantitative trait
# sd_trait: standard deviation of the studied quantitative trait
# thres_raw: cut-off for the clincially-defined categories
# eaf: effect allele frequency
# rank_normalized_y: normalized the studied quantitative trait
# retured value: beta from logistic regression based on clinically-defined categories
#**********************************************************************************************
Binary_GWAS <- function(beta_cont, mean_trait,  sd_trait, thres_raw,  eaf , rank_normalized_y ) {  ##note that se of the coef is not required
  
  if (rank_normalized_y==TRUE) {
    thres = (thres_raw - mean_trait) / sd_trait  ##standarize the threshold on a normal scale if the trait is already standardized
    mean_trait <-  0 
    sd_trait <- 1
  } else {thres = thres_raw} 
  p <- eaf 
  mean_X = 2*p*(1-p)  + 2* p^2  #X is the genotype originally coded as 0,1,2
  #****************************************************************************************
  # Reference for the calculation of sd_X:#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2668004/#app1
  # Suppose y is the posterior probability vector of a genotype(coded as 0,1,2)
  # According to the reference:
  # Var(X) = E(X*X|y)-E(X|y)*E(X|y) 
  # E(X*X|y) = y(1)+4y(2) ;  E(X|y)*E(X|y) = [y(1)+2y(2)]*[y(1)+2y(2)]
  # Var(X) = 2*p*(1-p) + 4*p*p - (2*p)*(2*p)
  #****************************************************************************************
  sd_X = sqrt( 2*p*(1-p) )  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2668004/#app1
  
  ##b_standardized = beta_cont * sd_X / sd_trait 
  beta1 = beta_cont 
  
  ##variance explained by SNP X =  beta1^2 * var(X)
  var.given.SNP =   sd_trait^2 - beta1^2 * sd_X^2  #total_var(y) - beta1^2 * var(X) #the residual variance of y given the genotype, usu. this is very close to the total var. of y, as individual SNPs have very small var explained
  beta0 = mean_trait - beta1*mean_X
  
  Etrait.aa = beta0
  Pr.trait.aa = pnorm( thres, mean = Etrait.aa, sd = sqrt(var.given.SNP) , lower.tail=FALSE) 
  
  Etrait.Aa = beta0 + beta1*1
  Pr.trait.Aa = pnorm( thres, mean = Etrait.Aa, sd = sqrt(var.given.SNP) , lower.tail=FALSE) 
  
  Etrait.AA = beta0 + beta1*2
  Pr.trait.AA = pnorm( thres, mean = Etrait.AA, sd = sqrt(var.given.SNP) , lower.tail=FALSE) 
  
  ##computing odds ratio 
  odds_Aa = Pr.trait.Aa/ (1- Pr.trait.Aa) 
  odds_aa = Pr.trait.aa/ (1- Pr.trait.aa) 
  OR_Aa_to_aa = odds_Aa / odds_aa 
  
  odds_AA = Pr.trait.AA/ (1- Pr.trait.AA) 
  OR_AA_to_aa = odds_AA / odds_aa 
  
  PAa = 2*p*(1-p) 
  PAA = p^2 
  
  beta_binary = PAa/(PAa+PAA) * log(OR_Aa_to_aa) + PAA/(PAa+PAA) * log(OR_AA_to_aa)/2 
  return (beta_binary)
}

#**********************************************************************************************
# Function: transform se derived from linear regression to that from logistic regression
# beta_cont: beta derived from linear regression
# se_cont: se derived from linear regression
# mean_trait: mean value of the studied quantitative trait
# sd_trait: standard deviation of the studied quantitative trait
# thres_raw: cut-off for the clincially-defined categories
# eaf: effect allele frequency
# rank_normalized_y: normalized the studied quantitative trait
# retured value: se from logistic regression based on clinically-defined categories
#**********************************************************************************************
Binary_GWAS_SE <- function(beta_cont, se_cont, mean_trait,  sd_trait, thres_raw,  eaf , rank_normalized_y = TRUE) {
  library(numDeriv)
  deriv = grad(func=Binary_GWAS, 
               x = beta_cont, 
               mean_trait = mean_trait,  sd_trait= sd_trait, thres_raw=thres_raw, 
               eaf = eaf, 
               rank_normalized_y = TRUE,
               method="Richardson") 
  
  delta_var = deriv^2 * (se_cont)^2 
  delta_se = sqrt( delta_var ) 
  return(delta_se) 
}
