################################################################################
##### Example Scripts for EpigenTL Functions                               #####
##### Updated February 12, 2024                                            #####
##### By: Kristen McGreevy                                                 ##### 
################################################################################



################################ Contents ######################################

## Loading Files, Functions, and Packages
## Saliva.2.Blood.DNAmBiomarkers Function
## Example 1. Calculating C Algorithms in Saliva DNAm 
## Example 2. Calculating C+S Algorithms in Saliva DNAm 
## Epigen.TL.Lasso Function
## Data Setup for EpigenTL.Lasso Function
## Example 3. Oracle 1df EpigenTL Lasso 
## Example 4. Oracle EpigenTL Lasso 
## Example 5. Estimate A0 EpigenTL Lasso 




################################################################################
################## Loading Files, Functions, and Packages ######################

### In order to run these functions, you must have glmnet and dplyr packages. 
library(glmnet)
library(dplyr)


### load the files from where you downloaded the files to.
load("~/ExampleFiles.RData")
# these load 5 items into your global environment:
# X_matrix_Github, Y_matrix_GitHub, n_vector_GitHub, 
# ExSample_SalivaCpGs, ExSample_SalivaDNAm


source("~/EpigenTL_SourceFunctions.R")
# this will load helper functions, background functions, and the main functions
# into your global environment. 
# the only functions you need to interact with are 
# Saliva.2.Blood.DNAmBiomarkers and
# Epigen.TL.Lasso


################################################################################
################## Saliva.2.Blood.DNAmBiomarkers Function ######################

# Saliva.2.Blood.DNAmBiomarkers function requires 2 or 3 inputs depending on 
# if you are calculating the C or C+S algorithms. 

## X: is the matrix or dataframe of saliva DNA methylation beta values with 
## samples in rows and methylation sites in columns. To see the list of CpG 
## loci needed for computation, please call colnames to CS_Algorithms_GitHub or 
## C_Algorithms_GitHub. 
## method: should be either "C+S" or "C" to indicate which set of algorithms 
## you are interested in calculating. If "C+S", you must also provide SalivaDNAmBiom. 
## SalivaDNAmBiom: should be the matrix or dataframe of Saliva DNAm biomarkers 
## (ie DNAm biomarkers directly calculated with saliva methylation values) if
## using the C+S algorithms. Otherwise, this should be NULL. Default is NULL. 


## This function outputs a dataframe with predicted DNAm Biomarkers in each 
## column with rows in the same order as the originally supplied X. 



################################################################################
############ Example 1. Calculating C Algorithms in Saliva DNAm ################

# ExSample_SalivaCpGs is loaded in your global environment. It has saliva
# DNA methylation values for 109 samples. The CpG sites are in columns 

## Function Call
dim(ExSample_SalivaCpGs) # 109 samples and 6662 CpG loci (Only 1307 needed though)
CAlgo_SalivaPred <- Saliva.2.Blood.DNAmBiomarkers(X = ExSample_SalivaCpGs, 
                                              method = "C", SalivaDNAmBiom = NULL) 

# looking at some of the predictions, which have the biomarker name and _C_Pred. 
head(CAlgo_SalivaPred)
dim(CAlgo_SalivaPred) # 19 biomarkers for all 109 samples




################################################################################
########### Example 2. Calculating C+S Algorithms in Saliva DNAm ###############

## To Calculate C+S Algorithms, you also need the matrix of Saliva DNAm Biomarkers
## which are the calculated DNAm biomarkers using saliva DNAm values

colnames(ExSample_SalivaDNAm)

# Function Call
CSAlgo_SalivaPred <- Saliva.2.Blood.DNAmBiomarkers(X = ExSample_SalivaCpGs, 
                                              method = "C+S", SalivaDNAmBiom = ExSample_SalivaDNAm) 

# looking at some of the predictions, which have the biomarker name and _CS_Pred. 
head(CSAlgo_SalivaPred)
dim(CSAlgo_SalivaPred) # 9 biomarkers for all 109 samples




################################################################################
######################### Epigen.TL.Lasso Function #############################


## This performs Transfer Learning Lasso allowing either the Auxiliary Dataset 
## information to be known (Oracle or Oracle 1df) or estimated in the process 
## (Estimate A0). 

## X_matrix: is the matrix of covariates (nxp) concatenated between the Target 
## and all Auxiliary Datasets, in that order. 
## n should therefore be n_0 + n_1 + ... + n_k.
## Please make sure that X_matrix is a matrix and not a dataframe to allow for 
## matrix multiplication to be carried out in the function. 
## Y_vector: is the outcome vector to develop the TL Lasso model to estimate. It should 
## align with your X_matrix with Target and Auxiliary outcomes concatenated and be a nx1 vector.
## N_vector: is a vector of number of observations in the Target and each Auxiliary
## datasets in the order they are concatenated. 
## AuxInformation: is either "Estimate A0", "Oracle", or "Oracle 1df" to signify the 
## informative auxiliary datasets need estimated  or they should be treated
## as known. If known, they can either be treated as equally informative and combined 
## into 1 dataset (Oracle), or they can be treated individually (Oracle 1df).
## RhatCount: is either "n0/3" or a value between 1, ..., p to specify how many
## marginal correlations to consider when calculating information. This should be set to NULL
## if AuxInformation is "Oracle" or "Oracle 1df".
## LambdaType: is either "Constant" or "CV" to indicate if the lambda calculated 
## from the target dataset should be used and adjusted based on the aux dataset size ("Constant")
## or whether the optimal lambda should be calculated via cross validation ("CV") for
## each set of auxiliary information. Default is "CV".
## seedstart: is the seed to set in the calculation for reproducibility. If not specified,
## it is set to 123.


## This function will return a data.frame with TL coefficients in the columns. 
## The first column, "Variable" labels the intercept and column from your X matrix
## that it corresponds to. Columns 2:7 or 2:4 have the final coefficients with slight 
## variations in their calculation. "min" and "1se" correspond to the lambda that 
## either minimizes CV error or is at most 1se above it. "allcoef", "halfcoef",
## and "lambcoef" correspond to the coefficient thresholding used before 
## combining coefficients across the auxiliary and target dataset. If "Constant"
## LambdaType was specified, only the minimum lambda from the target data is used
## and so "1se" coefficients are not presented.






################################################################################
################## Data Setup for EpigenTL.Lasso Function ######################

# To demonstrate how this function is called, I provide sample dataframes. 
# This example will show how we can use saliva DNAm to predict blood DNAm 
# biomarkers by incorporating information from multiple tissues and datasets. 

# The set of covariates has both Target and Aux data: X_matrix_Github
# This HAS to be organized where the Target data is provided in the first rows, 
# and Aux data are provided afterwards. 
# X_matrix_Github has 109 saliva samples, 27 buccal (GSE111165), 
# 19 buccal (GSE214901), and 136 adipose (TwinsUK) in that order.
# There are 1307 CpGs in the columns. 

dim(X_matrix_Github)
head(colnames(X_matrix_Github))

# The set of outcomes are stored in a matrix called Y_matrix_GitHub. To call
# this function, you must only provide one outcome at a time, so we will 
# use just 1 column in it for demonstration. 
# This is in the same order as the X_matrix, which corresponds to the outcomes
# for 109 blood, 27 blood, 19 brain, and 136 skin samples. 
# There are 45 DNAm Biomarkers in the columns. 
dim(Y_matrix_GitHub)
head(colnames(Y_matrix_GitHub))


# The counts in each sample correspond to values in n_vector_GitHub. 
# This tells the function which samples are Target or Aux. 
n_vector_GitHub





################################################################################
################### Example 3. Oracle 1df EpigenTL Lasso #######################

# Take 1 outcome of interest from the Y matrix to make the vector
Y_OutcomeVector <- Y_matrix_GitHub[, "DNAmAge"]


## This uses CV Lambda, meaning the optimal lambda is calculated within each Aux Data
DNAmAge_EpigenTL <- Epigen.TL.Lasso(X_matrix = X_matrix_Github, 
                                    Y_vector = Y_OutcomeVector, 
                                    N_vector = n_vector_GitHub, 
                                    AuxInformation = "Oracle 1df", 
                                    RhatCount = NULL, seedstart = 253)

# Coefficients are stored in DNAmAge_EpigenTL

# Coefficients that use the min Lambda and Don't have Coefficient Thresholding:
sum(DNAmAge_EpigenTL$beta_min_allcoef != 0) # 161 coefficients 

# If we instead take Coefficients Thresholded at Lambda level, fewer coefficients are retained:
sum(DNAmAge_EpigenTL$beta_min_lambcoef != 0) # 88 coefficients

# look at the first few CpGs that were chosen and their coefficient value. 
head(DNAmAge_EpigenTL[DNAmAge_EpigenTL$beta_min_allcoef != 0, c("Variable", "beta_min_allcoef")])




################################################################################
#################### Example 4. Oracle EpigenTL Lasso  #########################

# Take 1 outcome of interest from the Y matrix to make the vector
Y_OutcomeVector <- Y_matrix_GitHub[, "DNAmPAI1"]

# Take all Aux Data as 1 informative dataset and use Constant Lambda 
# which uses the target data to calculate the lambda and updates this value. 
DNAmPAI1_EpigenTL_ConstLamb <- Epigen.TL.Lasso(X_matrix = X_matrix_Github, 
                                    Y_vector = Y_OutcomeVector, 
                                    N_vector = n_vector_GitHub, 
                                    AuxInformation = "Oracle", 
                                    RhatCount = NULL, LambdaType = "Constant", 
                                    seedstart = 253)

sum(DNAmPAI1_EpigenTL_ConstLamb$beta_min_allcoef != 0) # 81 

# look at the first few CpGs that were chosen and their coefficient value. 
head(DNAmPAI1_EpigenTL_ConstLamb[DNAmPAI1_EpigenTL_ConstLamb$beta_min_allcoef != 0, 
                                 c("Variable", "beta_min_allcoef")])



# Take all Aux Data as 1 informative dataset and use CV Lambda 
# which calculates the optimal lambda for Lasso within target and aux data. 
DNAmPAI1_EpigenTL_CVLamb <- Epigen.TL.Lasso(X_matrix = X_matrix_Github, 
                                    Y_vector = Y_OutcomeVector, 
                                    N_vector = n_vector_GitHub, 
                                    AuxInformation = "Oracle", 
                                    RhatCount = NULL, LambdaType = "CV", 
                                    seedstart = 253)

sum(DNAmPAI1_EpigenTL_CVLamb$beta_min_allcoef != 0) # 73 coefficients instead. 




################################################################################
################## Example 5. Estimate A0 EpigenTL Lasso  ######################


# Instead of treating all the Aux data as informative, I will let the algorithm
# figure out which ones are helpful for predicting the biomarkers across tissues. 
# Now I have to supply RhatCount, which is how many columns in X to use to 
# calculate informativeness

# Take 1 outcome of interest from the Y matrix to make the vector
Y_OutcomeVector <- Y_matrix_GitHub[, "DNAmFitAge"]

DNAmFitAge_EpigenTL_EstA0 <- Epigen.TL.Lasso(X_matrix = X_matrix_Github, 
                                             Y_vector = Y_OutcomeVector, 
                                             N_vector = n_vector_GitHub, 
                                             AuxInformation = "Estimate A0", 
                                             RhatCount = 500, LambdaType = "CV", 
                                             seedstart = 269)
# Coefficients are stored in DNAmFitAge_EpigenTL_EstA0


# Coefficients that use the min Lambda and Don't have Coefficient Thresholding:
sum(DNAmFitAge_EpigenTL_EstA0$beta_min_allcoef != 0) # 188 non-zero coeff

# Coefficients that use the min Lambda and have Coefficient Thresholding at Lambda:
sum(DNAmFitAge_EpigenTL_EstA0$beta_min_lambcoef != 0) # 40 non-zero coef 


# look at the first few CpGs that were chosen and their coefficient value. 
head(DNAmFitAge_EpigenTL_EstA0[DNAmFitAge_EpigenTL_EstA0$beta_min_lambcoef != 0, 
                                 c("Variable", "beta_min_lambcoef")])


