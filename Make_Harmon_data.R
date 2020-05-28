library(TwoSampleMR)
library(data.table)

#********************************************************************************************
# Function: harmonize GWAS summary statistics for two different traits
# D1: GWAS file name for the 1st trait excluding the file type name
# D2: GWAS file name for the 2nd trait excluding the file type name
# harmon_data_file: filepath for the harmonized summary statistics of two traits
#********************************************************************************************
make_harmon_data <- function(D1,D2,harmon_data_file){

	#library(TwoSampleMR)
	D1_df <- fread(paste0(D1,'.txt'))
	D2_df <- fread(paste0(D2,'.txt'))
	D1_df$Phenotype = D1
	D2_df$Phenotype = D2
	exposure    <- format_data(D1_df,type='exposure',effect_allele_col='A1',other_allele_col='A2',eaf_col='A1_Freq')
	outcome     <- format_data(D2_df,type='outcome', effect_allele_col='A1',other_allele_col='A2',eaf_col='A1_Freq')
	harmon_data <- harmonise_data(exposure_dat=exposure, outcome_dat=outcome, action=2)
	harmon_data <-subset(harmon_data, mr_keep == TRUE & beta.outcome !=0 & beta.exposure != 0 & !is.na(eaf.outcome) & !is.na(eaf.exposure)); #nrow(harmon_data)
	fwrite(harmon_data,harmon_data_file,sep='\t')
}
