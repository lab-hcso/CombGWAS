
#*****************************************************************************************************************
# Function: Infer summary statistics of comorbid disorder from corresponding traits, eg, SCZ and BIP(Bipolar) 
# Parameters:
# harmon_data: harmonized summary statistics for two traits derived from TwosampleMR
# Size_TraitA: Sample size of trait A
# Size_TraitB: Sample size of trait B
# TraitA_preval: Disease prevalence of trait A
# TraitB_preval: Disease prevalence of trait B
# Comorbid_preval: Disease prevalence of comorbid trait
# LDSR_intercept : LD score regression intercept  
# Notes: Summary statistics are based on the shared set of SNPs for two traits
#*****************************************************************************************************************
GWISForComorbidity <- function(harmon_data,Size_TraitA,Size_TraitB,TraitA_preval,TraitB_preval,Comorbid_preval,LDSR_intercept){
  
  library(locfdr)
  library(sfsmisc)
  library(msm)
  library(rootSolve)
  library(numDeriv)
  
  #**************************************************************************************************
  # k.1st: disease prevalence of 1st trait
  # k.2nd: disease prevalence of 2nd trait
  # LDSR_intercept: LD score regression intercept
  # beta.logistic.1st: logistic regression coefficient for the 1st trait
  # beta.logistic.2nd: logistic regression coefficient for the 2nd trait
  # se.beta.logistic.1st: logistic regression standard error for the 1st trait
  # se.beta.logistic.2nd: logistic regression standard error for the 2nd trait
  # p1: effect allele frequency for the 1st trait
  # p2: effect allele frequency for the 2nd trait
  #**************************************************************************************************
  k.1st = TraitA_preval
  k.2nd = TraitB_preval
  comorbid.preval = Comorbid_preval
  LDSR_intercept <- LDSR_intercept
  # N is the estimated sample size of the comorbidity
  N = (Size_TraitA + Size_TraitB)/2
  totalSNPs <- nrow(harmon_data)
  
  # Cal weighted-average A1_Freq
  cal_weighted_avg_eaf<-function(harmon_data,Size_TraitA,Size_TraitB){
    denom  <- sqrt(Size_TraitA) + sqrt(Size_TraitB)
    avg_eaf<- sqrt(Size_TraitA)/denom * harmon_data$eaf.exposure + 
      sqrt(Size_TraitB)/denom * harmon_data$eaf.outcome
    cbind(harmon_data, avg_eaf)
  }
  harmon_data <- cal_weighted_avg_eaf(harmon_data,Size_TraitA,Size_TraitB)
  
  # Result Variables
  beta1.comor.logistic        <- rep(NA,totalSNPs)
  new.se.beta1.comor.logistic <- rep(NA,totalSNPs)
  zscore_logistic             <- rep(NA,totalSNPs)
  pval_logistic_2tailed       <- rep(NA,totalSNPs)
  pval_logistic               <- rep(NA,totalSNPs)
  
  for(i in 1:totalSNPs)tryCatch({
    
    beta.logistic.1st = harmon_data$beta.exposure[i]
    beta.logistic.2nd = harmon_data$beta.outcome[i]
    se.beta.logistic.1st = harmon_data$se.exposure[i]
    se.beta.logistic.2nd = harmon_data$se.outcome[i]
    p1 = harmon_data$eaf.exposure[i]
    p2 = harmon_data$eaf.outcome[i]
    p  = harmon_data$avg_eaf[i]
    
    var_1st = k.1st*(1-k.1st)
    var_2nd = k.2nd*(1-k.2nd)
    cov_1st_2nd = comorbid.preval - k.1st*k.2nd
    
    #*********************************************
    # directly solve quadratic equation
    #*********************************************
    # Constructing Quadratic Formula from https://rpubs.com/kikihatzistavrou/80124
    
    ## here we modify the code such that only beta in the correct range is allowed
    # the linear corrected beta should be closer to zero than the original one 
    
    # Constructing delta
    delta<-function(a,b,c){
      b^2-4*a*c
    }
    
    quad_equation_result <- function(a,b,c){
      if(delta(a,b,c) > 0){ # first case D>0
        x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
        x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
        result = c(x_1,x_2)
      }
      else if(delta(a,b,c) == 0){ # second case D=0
        result =  -b/(2*a)
      }
      else {result <- NA} # third case D<0
      
      return(result)
    }
    
    #*******************************************
    # convert logistic beta to linear beta
    #*******************************************
    funct_equation <- function(beta.logistic, k, p)  {
      a = exp(beta.logistic)
      quad_a = p*(1-p) - a*p*(1-p)
      quad_b = a*k*(1-p)+ a*p -a*p*k + k*p + (1-p)*(1-k) 
      quad_c = -(   a*k-a*k^2 -k*(1-k)   )
      quad_solution = quad_equation_result(quad_a, quad_b, quad_c) 
      b1 = quad_solution[abs(quad_solution) < abs(beta.logistic)]
      return(b1)
    }
    
    ## b0 = y_mean - b1* x_mean
    ##b0 is the beta0 in linear regression
    b0.func <- function(b1,k,p) {
      k - b1*(2*p)
    }
    
    #*****************
    # 1st trait 
    #*****************
    b1_1st = funct_equation(beta.logistic = beta.logistic.1st, 
                            k = k.1st, 
                            p = p1)
    b0_1st = b0.func(b1 =b1_1st , k= k.1st, p=p1)

    funct_equation_grad <- function(x)  {
      beta.logistic <- x[1]
      k <- x[2]
      p <- x[3]
      a = exp(beta.logistic)
      quad_a = p*(1-p) - a*p*(1-p)
      quad_b = a*k*(1-p)+ a*p -a*p*k + k*p + (1-p)*(1-k) 
      quad_c = -(   a*k-a*k^2 -k*(1-k)   )
      quad_solution = quad_equation_result(quad_a, quad_b, quad_c) 
      b1 = quad_solution[abs(quad_solution) < abs(beta.logistic)]
      return(b1)
    }
    
    gradient = grad(funct_equation_grad, 
                    c(beta.logistic.1st, k.1st, p1),
                    method="Richardson")  ##substitute the 2nd number with the observed logistic regresion beta1
    
    cov_matrix_1st = matrix (c(se.beta.logistic.1st^2, 0,          0,
                               0, k.1st*(1-k.1st)/Size_TraitA,  0,
                               0,  0       , p1*(1-p1)/(2*Size_TraitA)), nrow=3, byrow=TRUE) 
    
    var.b1.linear.1st = matrix(gradient, nrow=1) %*% cov_matrix_1st %*% matrix(gradient, nrow=3)
    ## to derive variance of beta0
    # https://stats.stackexchange.com/questions/64195/how-do-i-calculate-the-variance-of-the-ols-estimator-beta-0-conditional-on
    #var(beta0) = sum(xi^2)/n * var(beta1)
    var.b0.linear.1st = (2*p1*(1-p1) + 4*p1^2 )*var.b1.linear.1st
    
    
    #*****************
    # 2nd trait 
    #*****************
    b1_2nd = funct_equation(beta.logistic = beta.logistic.2nd, 
                            k = k.2nd, 
                            p = p2)
    b0_2nd = b0.func(b1 =b1_2nd , k= k.2nd, p=p2)
    
    gradient = grad(funct_equation_grad, 
                    x = c(beta.logistic.2nd, k.2nd, p2),
                    method="Richardson")  ##substitute the 2nd number with the observed logistic regresion beta1
    
    cov_matrix_2nd = matrix (c(se.beta.logistic.2nd^2, 0,          0,
                               0, k.2nd*(1-k.2nd)/Size_TraitB,  0,
                               0,  0       , p2*(1-p2)/(2*Size_TraitB)), 
                             nrow=3, byrow=TRUE) 
    
    var.b1.linear.2nd = matrix(gradient, nrow=1) %*% cov_matrix_2nd %*% matrix(gradient, nrow=3)  
    
    
    
    ## to derive variance of beta0
    # https://stats.stackexchange.com/questions/64195/how-do-i-calculate-the-variance-of-the-ols-estimator-beta-0-conditional-on
    #var(beta0) = sum(xi^2)/n * var(beta1)
    var.b0.linear.2nd = (2*p2*(1-p2) + 4*p2^2 )*var.b1.linear.2nd
    
    
    #*************************************************************************
    #  now apply the method of GWIS 2ND ORDER expansion (beta_1st_trait * beta_2nd_trait to test for comorbidity)
    #  the 2nd order formula is the same as the 1st order formula 
    #*************************************************************************
    
    gamma0 = b0_1st * b0_2nd + cov_1st_2nd        
    
    
    AF1 = 2*p*(1-p) / (2*p*(1-p) + p^2)
    AF2 = ((p^2) ) / (2*p*(1-p) + p^2)
    
    beta1.comor = 
      AF1* (  (b0_1st + b1_1st)  *(b0_2nd + b1_2nd) + cov_1st_2nd  - gamma0 ) + 
      AF2/2*( (b0_1st + 2*b1_1st)*(b0_2nd + 2*b1_2nd) + cov_1st_2nd - gamma0)
    #**************************************
    # compute new SE by delta method
    #**************************************
    ## note that here we consider var and covar of all four b0_1st (x1), b1_1st (x2), b0_2nd (x3) and b1_2nd (x4)
    ## covariance are only present between x2 and x4
    
    # http://www.stat.columbia.edu/~fwood/Teaching/w4315/Fall2009/homework_1_solution.pdf
    # https://www.unige.ch/sciences/astro/files/5413/8971/4090/2_Segransan_StatClassUnige.pdf
    # https://stats.stackexchange.com/questions/64195/how-do-i-calculate-the-variance-of-the-ols-estimator-beta-0-conditional-on
    # https://stats.stackexchange.com/questions/38721/covariance-of-a-variable-and-a-linear-combination-of-other-variables
    # 
    # https://stats.stackexchange.com/questions/127316/clarification-the-covariance-of-intercept-and-slope-in-simple-linear-regression
    #cov of b0 and b1 in a regression : cov(b0,b1) = -x_mean * var(b1)
    
    ###need to run LD score regression 1st
    cov.b1_1st.b1_2nd = sqrt(var.b1.linear.1st) * LDSR_intercept *sqrt(var.b1.linear.2nd)
    
    ##for the SAME trait
    cov.b0_1st.b1_1st  = -(2*p1)* var.b1.linear.1st #? -(2*p1) #cov of b0 and b1 in a regression : cov(b0,b1) = -x_mean * var(b1)
    cov.b0_2nd.b1_2nd  = -(2*p2)* var.b1.linear.2nd #? -(2*p2)
    
    cov.b0_1st.b1_2nd = -(2*p)*cov.b1_1st.b1_2nd #? -(2*sqrt(p1*p2))
    cov.b0_2nd.b1_1st = -(2*p)*cov.b1_1st.b1_2nd #? -(2*sqrt(p1*p2))
    
    cov.b0_1st.b0_2nd = -(2*p)^2*cov.b1_1st.b1_2nd #? -(2*sqrt(p1*p2))
    
    
    covmat = matrix( c(var.b0.linear.1st, cov.b0_1st.b1_1st, cov.b0_1st.b0_2nd, cov.b0_1st.b1_2nd,
                       cov.b0_1st.b1_1st, var.b1.linear.1st, cov.b0_2nd.b1_1st, cov.b1_1st.b1_2nd,
                       cov.b0_1st.b0_2nd, cov.b0_2nd.b1_1st, var.b0.linear.2nd, cov.b0_2nd.b1_2nd,
                       cov.b0_1st.b1_2nd, cov.b1_1st.b1_2nd, cov.b0_2nd.b1_2nd, var.b1.linear.2nd),
                     nrow = 4, ncol=4 , byrow=FALSE
    )
    
    #CORRECT VERSION
    formula1 = sprintf("~ %f* (  (x1 + x2)  *(x3 + x4)  - x1*x3  ) + 
                       %f/2*(  (x1 + 2*x2)*(x3 + 2*x4)- x1*x3  )", AF1, AF2) 
    new.se.beta1.comor <- deltamethod( as.formula(formula1) , 
                                       mean= c(b0_1st, b1_1st, b0_2nd, b1_2nd) ,
                                       cov = covmat 
    )
    #******************************************************
    # now convert the comorbid beta from linear back to logistic scale
    #**********************************************
    
    converion_to_logisticBeta <- function(betaLinear, k , p) {
      numer = (k + betaLinear*(1-p)) * (1-k + betaLinear*p)
      denom = (k - betaLinear*p) * (1-k-betaLinear*(1-p))
      log(numer/denom)
    }
    
    beta1.comor.logistic[i] = converion_to_logisticBeta(beta1.comor, k = comorbid.preval, p = p )
    
    formula_new = sprintf("~ log( (x2 + x1*(1-x3)) * (1-x2 + x1*x3) / (x2 - x1*x3) / (1-x2-x1*(1-x3)) )") 
    cov_matrix_comor = matrix (c(new.se.beta1.comor^2, 0,          0,
                                 0, comorbid.preval*(1-comorbid.preval)/N,  0,
                                 0,  0       , p*(1-p)/(2*N) ), 
                               nrow=3, byrow=TRUE)  #? N = (N1+N2)/2
    
    new.se.beta1.comor.logistic[i] = 
      deltamethod( as.formula(formula_new), 
                   mean= c(beta1.comor ,comorbid.preval, p) , 
                   cov = cov_matrix_comor
      ) 
    
    zscore_linear = beta1.comor/new.se.beta1.comor
    
    ##note we need one-tailed test (we assume the direction of association should be the same for each trait)
    if (beta1.comor>0) {
      pval_linear =    pnorm( -abs(zscore_linear) )
    } 
    else {
      pval_linear = 1- pnorm( -abs(zscore_linear) ) 
    }

    zscore_logistic[i] = beta1.comor.logistic[i]/new.se.beta1.comor.logistic[i]
    pval_logistic_2tailed[i] = pnorm( -abs(zscore_logistic[i]) )*2
    
    if (beta1.comor>0) {
      pval_logistic[i] =  pnorm( -abs(zscore_logistic[i]) )
    } 
    else {
      pval_logistic[i] = 1- pnorm( -abs(zscore_logistic[i]) ) 
    }
  },error=function(e){
    cat('\ni=',i,'\n')
    print(harmon_data[i,])
  }) 
  
  
  rmse <- function(xpred, xreal){
    sqrt(  mean( (xpred-xreal)^2 )  ) 
  }
  
  Results <-data.frame(         SNP=harmon_data$SNP,
                                A1 = harmon_data$effect_allele.exposure,
                                A2 = harmon_data$other_allele.exposure,
                                trait1_A1F  = harmon_data$eaf.exposure,
                                trait2_A1F  = harmon_data$eaf.outcome,
                                trait1_beta = harmon_data$beta.exposure,
                                trait1_se   = harmon_data$se.exposure,
                                trait1_pval = harmon_data$pval.exposure,
                                trait2_beta = harmon_data$beta.outcome,
                                trait2_se   = harmon_data$se.outcome,
                                trait2_pval = harmon_data$pval.outcome,
                                avg_eaf     = harmon_data$avg_eaf,
                                Comorb_beta = beta1.comor.logistic, 
                                Comorb_se   = new.se.beta1.comor.logistic,
                                Comorb_pval = pval_logistic_2tailed
  )
  setorder(Results,'Comorb_pval')
  return (Results)
}
