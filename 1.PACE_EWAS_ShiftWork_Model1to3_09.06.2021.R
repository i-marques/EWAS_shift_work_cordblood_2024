############################################################## INSTRUCTIONS ##########################################################


# "##### PLEASE CHANGE" - this will give you instructions when there is the need to adapt this script


############################################################## LOAD LIBRARIES #######################################################

library(foreign)
library(data.table)# to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) #Huber?s estimation of the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(plyr) # to process results
library(matrixStats) #for rowSums, rowIQRs etc.
library(tableone)
library(limma)
library(stats)


############################################################## STUDY PARAMETERS #######################################################

####### PLEASE CHANGE cohort and analysis.date only 

cohort <- "CohortNAME" #change to name of your cohort/study
analysis.date <- "20210521" # change to the date at which you perform the analyses


## Names for models, you do not need to change anything here##
basic.model <- "basic.model" 
main.model <- "main.model" 
main.cells <- "main.cells"


############################################################## MODELS ###################################################################


############################################################## SET MODEL PARAMETERS #######################################################


####### Please make sure these names are identical to the names in your .dat phenotype file! #######


sample.id <-c("sample.id") ### Case identifier
trait <-c("mat.shiftwork") ### The exposure variables of interest --> (in our case maternal shift work)
cell.names <- c("bcell", "mono", "cd4t", "cd8t", "gran", "nk", "nRBC") ### Salas reference set for cell type correction
batch <-c("batch") ### Adjust for batch effects by including the most important covariate(s), such as plate, in the models
child.sex <- c("child.sex")
main.covariates <- c("mat.age", "mat.edu", "mat.smoking") ### Will be included in each model apart from the basic model


############################################################## PHENOTYPE FILE #######################################################

####### PLEASE CHANGE - Load phenotype data

## This should be a .dat file (tab delimited) containing information on: sample IDs (1st column), exposure (maternal shift work), cell counts (Salas reference, 7 columns), batch, maternal covariates, and child sex.
## There should be no extra columns


dataSW <- read.table("....dat", header = T, stringsAsFactors=FALSE) ### between quotes, specify the name of your phenotype file

names(dataSW)


## This part is to make sure that categorical covariates are treated as factors and continuous covariates as numeric.
dataSW$mat.shiftwork = as.factor(dataSW$mat.shiftwork)
dataSW$bcell = as.numeric(dataSW$bcell)
dataSW$mono = as.numeric(dataSW$mono)
dataSW$cd4t = as.numeric(dataSW$cd4t)
dataSW$cd8t = as.numeric(dataSW$cd8t)
dataSW$gran = as.numeric(dataSW$gran)
dataSW$nk = as.numeric(dataSW$nk)
dataSW$nRBC = as.numeric(dataSW$nRBC)
dataSW$batch = as.factor(dataSW$batch) 
dataSW$mat.age = as.numeric(dataSW$mat.age)
dataSW$mat.edu = as.factor(dataSW$mat.edu)
dataSW$mat.smoking = as.factor(dataSW$mat.smoking)
dataSW$child.sex = as.factor(dataSW$child.sex)


## Summarizes your data before excluding any NA's
summary(dataSW)



## Check if all variables are present in phenotype file
for(i in 1:length(c(sample.id, trait, cell.names, batch, child.sex, main.covariates))) {
  print(ifelse(c(sample.id, trait, cell.names, batch, child.sex, main.covariates)[i] %in% colnames(dataSW)==FALSE,
               paste("CAUTION: the variable called",c(sample.id, trait, cell.names, batch, child.sex, main.covariates)[i],"is missing from pheno file"),
               paste("variable called",c(sample.id, trait, cell.names, batch, child.sex, main.covariates)[i],"is present in pheno file")))
}


############################################################## METHYLATION FILE + WISORIZING 2% #######################################################

####### PLEASE CHANGE load your methylation file between quotes

## This is a file with the methylation betas as you have them, with CpGs as rows and participants as columns. 
## ATTENTION: We expect this file NOT to be 3IQR-trimmed


methylation <- "..."

methylation

str(methylation)

#Load Methylation Data File
load(methylation)


n <- nrow(x)         # number of lines (genes)
p <- ncol(x)         # number of columns (conditions)

n
p

class(x)
dim(x)

methylation <- as.matrix(x)  

str(methylation) 


##################### Functions to winsorize data

winsorize <- function(methylation,pct=winsorize.pct) {
    quantiles <- matrixStats::rowQuantiles(methylation, probs=c(pct,1-pct), na.rm=T)
    low <- quantiles[,1]
    upper <- quantiles[,2]

    outliers.lower <- rowSums(methylation < low, na.rm=T)
    outliers.upper <- rowSums(methylation > upper, na.rm=T)
    
    idx <- which(methylation < low, arr.ind=T)
    methylation[idx] <- low[idx[,1]]
    
    idx <- which(methylation > upper, arr.ind=T)
    methylation[idx] <- upper[idx[,1]]

    n <- rowSums(!is.na(methylation))
    log <- data.frame(outliers.lower, outliers.upper, n)
    
    return(list(methylation=methylation, log=log))
}


###### Replace outliers using winsorizing 2% (1% each side)

replace.outliers <- winsorize(methylation, 0.01)
new.methylation <- replace.outliers$methylation
outlier.log <- replace.outliers$log



###### Save the log


write.csv(outlier.log, file=paste0(cohort, "_N_Winsorize0.01_Log", ".csv"), quote=FALSE, row.names=F)


######## 


str(new.methylation)


n <- nrow(new.methylation)         # number of lines (genes)
p <- ncol(new.methylation)         # number of columns (conditions)

n
p

class(new.methylation)
dim(new.methylation)


## Rename methylation file as betaFINAL and make sure it is a matrix
betaFINAL <- new.methylation 
betaFINAL <- as.matrix(betaFINAL)


## Check N cases with methylation data available (before excluding NA's) and order the files according to ID
index<-which(colnames(betaFINAL) %in% dataSW$sample.id)
length(index)
beta<-betaFINAL[,index]
beta<-beta[,order(colnames(beta))]
dataSW<-dataSW[order(dataSW$sample.id),] 
ncol(beta)
nrow(beta)
betaFINAL<-beta



############################################################## MODEL 1: BASIC MODEL #######################################################


## Basic model: remove cases with NA for any of mat.shiftwork and covariates
intersecting.samples <- intersect(dataSW$sample.id,colnames(betaFINAL))
dataSW_basic <- na.omit(dataSW[which(dataSW$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.shiftwork", batch, child.sex))])
k_basic = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_basic] = rowMedians(betaFINAL, na.rm=TRUE)[k_basic[,1]]

summary(dataSW_basic) ## Summarizes data for the basic model
dim(dataSW_basic) 


## Basic model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataSW_basic$sample.id)
length(index)
beta_basic<-betaFINAL[,index]
beta_basic<-beta_basic[,order(colnames(beta_basic))]
dataSW_basic<-dataSW_basic[order(dataSW_basic$sample.id),]
ncol(beta_basic)
nrow(beta_basic)
all.equal(as.character(dataSW_basic$sample.id),as.character(colnames(beta_basic))) #Return must be TRUE, not FALSE!


betaFINAL_basic <- beta_basic

## Transpose the methylation data (rows are newborns and columns are CpGs)
beta_matrix_basic<-t(betaFINAL_basic)
dim(beta_matrix_basic)
all.equal(as.character(dataSW_basic$sample.id),rownames(beta_matrix_basic))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataSW_basic <- dataSW_basic[dataSW_basic$sample.id %in% intersect(dataSW_basic$sample.id, rownames(beta_matrix_basic)),]
dataSW_basic <- dataSW_basic[match(rownames(beta_matrix_basic), dataSW_basic$sample.id),]
sum(dataSW_basic$sample.id == rownames(beta_matrix_basic))/nrow(beta_matrix_basic) #check that the IDs are in the same order, should return 1
nprobes_basic = ncol(beta_matrix_basic)

### run EWA for model 1: Basic model
results_basic = data.frame(probeID_basic=colnames(beta_matrix_basic),
                           beta=rep(0, times=nprobes_basic),
                           se=rep(0, times=nprobes_basic),
                           p_val=rep(0, times=nprobes_basic),
                           n=rep(0, times=nprobes_basic))

for(i in 1:nprobes_basic){
  tryCatch({
    CpG_basic = beta_matrix_basic[,i]
    rlm.fit = rlm(CpG_basic ~ mat.shiftwork +
                    batch +
                    child.sex,
                  data=dataSW_basic,
                  maxit=200) 
    test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
    beta_basic = test[2,"Estimate"]
    SE_basic = test[2,"Std. Error"]
    PVAL_basic = test[2,"Pr(>|z|)"]
    N_basic = length(rlm.fit$residual)
    
    set(results_basic, i, 2L, beta_basic)
    set(results_basic, i, 3L, SE_basic)
    set(results_basic, i, 4L, PVAL_basic)
    set(results_basic, i, 5L, N_basic)
  }, error = function(err) {
    message("Error in ", colnames(beta_matrix_basic)[i])
    set(results_basic, i, 2L, NA)
    set(results_basic, i, 3L, NA)
    set(results_basic, i, 4L, NA)
    set(results_basic, i, 5L, NA)
  })
  
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_basic, "%\n")}
}

cat("Sorting and saving result file for model 1 (Basic model)...\n\n")
results_basic = results_basic[order(results_basic$p_val),]
cat("Saving the EWAS results for model 1 (Basic model)...\n")

## saves results in .csv file
write.csv(results_basic, file=paste0(cohort,".",basic.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_basic)
cat("\n")
summary(results_basic) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a file:
lambda_basic = median(qchisq(results_basic$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 1 (Basic model) is: ", lambda_basic, "\n")
write.table(lambda_basic, file=paste0(cohort,".",basic.model,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summary as .csv file
basic_shiftwork.tableone <- as.data.frame(print(CreateTableOne(data=dataSW_basic[,-1],factorVars=c("mat.edu", "mat.smoking", "child.sex", "batch")),stringsAsFactors=FALSE))
write.csv(basic_shiftwork.tableone,file=paste0(cohort,".",basic.model,".summary.",analysis.date,".csv"))

print("Model 1: Basic model shift work analysis completed") ## This completes the analysis of model 1. 



############################################################## MODEL 2: MAIN MODEL #######################################################


## Main model: remove cases with NA for any of mataternal shift work and covariates
intersecting.samples <- intersect(dataSW$sample.id,colnames(betaFINAL))
dataSW_main.model <- na.omit(dataSW[which(dataSW$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.shiftwork", batch, "child.sex", main.covariates))])
k_main.model = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_main.model] = rowMedians(betaFINAL, na.rm=TRUE)[k_main.model[,1]]

summary(dataSW_main.model) ## Summarizes data for the Main model
dim(dataSW_main.model) 


## Main model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataSW_main.model$sample.id)
length(index)
beta_main.model<-betaFINAL[,index]
beta_main.model<-beta_main.model[,order(colnames(beta_main.model))]
dataSW_main.model<-dataSW_main.model[order(dataSW_main.model$sample.id),]
ncol(beta_main.model)
nrow(beta_main.model)
all.equal(as.character(dataSW_main.model$sample.id),as.character(colnames(beta_main.model))) #Return must be TRUE, not FALSE!


betaFINAL_main.model <- beta_main.model

## Transpose the methylation data (rows are newborns and columns are CpGs)
beta_matrix_main.model<-t(betaFINAL_main.model)
dim(beta_matrix_main.model)
all.equal(as.character(dataSW_main.model$sample.id),rownames(beta_matrix_main.model))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataSW_main.model <- dataSW_main.model[dataSW_main.model$sample.id %in% intersect(dataSW_main.model$sample.id, rownames(beta_matrix_main.model)),]
dataSW_main.model <- dataSW_main.model[match(rownames(beta_matrix_main.model), dataSW_main.model$sample.id),]
sum(dataSW_main.model$sample.id == rownames(beta_matrix_main.model))/nrow(beta_matrix_main.model) #check that the IDs are in the same order, should return 1
nprobes_main.model = ncol(beta_matrix_main.model)

### run EWA for model 2: main model
results_main.model = data.frame(probeID_main.model=colnames(beta_matrix_main.model),
                                beta=rep(0, times=nprobes_main.model),
                                se=rep(0, times=nprobes_main.model),
                                p_val=rep(0, times=nprobes_main.model),
                                n=rep(0, times=nprobes_main.model))

for(i in 1:nprobes_main.model){
  tryCatch({
    CpG_main.model = beta_matrix_main.model[,i]
    rlm.fit = rlm(CpG_main.model ~ mat.shiftwork +
                    batch +
                    mat.age + 
                    mat.edu +
                    mat.smoking +
                    child.sex,
                  data=dataSW_main.model,
                  maxit=200) 
    test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
    beta_main.model = test[2,"Estimate"]
    SE_main.model = test[2,"Std. Error"]
    PVAL_main.model = test[2,"Pr(>|z|)"]
    N_main.model = length(rlm.fit$residual)
    
    set(results_main.model, i, 2L, beta_main.model)
    set(results_main.model, i, 3L, SE_main.model)
    set(results_main.model, i, 4L, PVAL_main.model)
    set(results_main.model, i, 5L, N_main.model)
  }, error = function(err) {
    message("Error in ", colnames(beta_matrix_main.model)[i])
    set(results_main.model, i, 2L, NA)
    set(results_main.model, i, 3L, NA)
    set(results_main.model, i, 4L, NA)
    set(results_main.model, i, 5L, NA)
  })
  
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_main.model, "%\n")}
}

cat("Sorting and saving result file for model 2 (Main model)...\n\n")
results_main.model = results_main.model[order(results_main.model$p_val),]
cat("Saving the EWAS results for model 2 (Main model)...\n")

## saves results in .csv file
write.csv(results_main.model, file=paste0(cohort,".",main.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_main.model)
cat("\n")
summary(results_main.model) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a file:
lambda_main.model = median(qchisq(results_main.model$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 2 (main model) is: ", lambda_main.model, "\n")
write.table(lambda_main.model, file=paste0(cohort,".",main.model,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summary as .csv file
main.model_shiftwork.tableone <- as.data.frame(print(CreateTableOne(data=dataSW_main.model[,-1],factorVars=c("mat.edu", "mat.smoking", "child.sex", "batch")),stringsAsFactors=FALSE))
write.csv(main.model_shiftwork.tableone,file=paste0(cohort,".",main.model,".summary.",analysis.date,".csv"))

print("Model 2: main shift work analysis completed") ## This completes the analysis of model 2. 




############################################################## MODEL 3: MAIN MODEL + CELL COUNTS #######################################################


## Main models + cell counts : remove cases with NA for any of mat.shiftwork and covariates
intersecting.samples <- intersect(dataSW$sample.id,colnames(betaFINAL))
dataSW_main.cells <- na.omit(dataSW[which(dataSW$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.shiftwork", batch, cell.names, child.sex, main.covariates))])
k_main.cells = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_main.cells] = rowMedians(betaFINAL, na.rm=TRUE)[k_main.cells[,1]]

summary(dataSW_main.cells) ## Summarizes data for the main model + cell counts
dim(dataSW_main.cells) 


## Main model + cell counts: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataSW_main.cells$sample.id)
length(index)
beta_main.cells<-betaFINAL[,index]
beta_main.cells<-beta_main.cells[,order(colnames(beta_main.cells))]
dataSW_main.cells<-dataSW_main.cells[order(dataSW_main.cells$sample.id),]
ncol(beta_main.cells)
nrow(beta_main.cells)
all.equal(as.character(dataSW_main.cells$sample.id),as.character(colnames(beta_main.cells))) #Return must be TRUE, not FALSE!


betaFINAL_main.cells <- beta_main.cells

## Transpose the methylation data (rows are newborns and columns are CpGs)
beta_matrix_main.cells<-t(betaFINAL_main.cells)
dim(beta_matrix_main.cells)
all.equal(as.character(dataSW_main.cells$sample.id),rownames(beta_matrix_main.cells))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataSW_main.cells <- dataSW_main.cells[dataSW_main.cells$sample.id %in% intersect(dataSW_main.cells$sample.id, rownames(beta_matrix_main.cells)),]
dataSW_main.cells <- dataSW_main.cells[match(rownames(beta_matrix_main.cells), dataSW_main.cells$sample.id),]
sum(dataSW_main.cells$sample.id == rownames(beta_matrix_main.cells))/nrow(beta_matrix_main.cells) #check that the IDs are in the same order, should return 1
nprobes_main.cells = ncol(beta_matrix_main.cells)

### run EWA for model 3: main model +cell counts
results_main.cells = data.frame(probeID_main.cells=colnames(beta_matrix_main.cells),
                                            beta=rep(0, times=nprobes_main.cells),
                                            se=rep(0, times=nprobes_main.cells),
                                            p_val=rep(0, times=nprobes_main.cells),
                                            n=rep(0, times=nprobes_main.cells))

for(i in 1:nprobes_main.cells){
  tryCatch({
    CpG_main.cells = beta_matrix_main.cells[,i]
    rlm.fit = rlm(CpG_main.cells ~ mat.shiftwork +
                    batch +
                    bcell +
                    mono + 
                    cd4t +
                    cd8t + 
                    gran + 
                    nk + 
                    nRBC +
                    mat.age + 
                    mat.edu +
                    mat.smoking +
                    child.sex,
                  data=dataSW_main.cells,
                  maxit=200) 
    test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
    beta_main.cells = test[2,"Estimate"]
    SE_main.cells = test[2,"Std. Error"]
    PVAL_main.cells = test[2,"Pr(>|z|)"]
    N_main.cells = length(rlm.fit$residual)
    
    set(results_main.cells, i, 2L, beta_main.cells)
    set(results_main.cells, i, 3L, SE_main.cells)
    set(results_main.cells, i, 4L, PVAL_main.cells)
    set(results_main.cells, i, 5L, N_main.cells)
  }, error = function(err) {
    message("Error in ", colnames(beta_matrix_main.cells)[i])
    set(results_main.cells, i, 2L, NA)
    set(results_main.cells, i, 3L, NA)
    set(results_main.cells, i, 4L, NA)
    set(results_main.cells, i, 5L, NA)
  })
  
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_main.cells, "%\n")}
}

cat("Sorting and saving result file for model 3 (main.cells)...\n\n")
results_main.cells = results_main.cells[order(results_main.cells$p_val),]
cat("Saving the EWAS results for model 3 (main.cells)...\n")

## saves results in .csv file
write.csv(results_main.cells, file=paste0(cohort,".",main.cells,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_main.cells)
cat("\n")
summary(results_main.cells) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a file:
lambda_main.cells = median(qchisq(results_main.cells$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 3 (main.cells) is: ", lambda_main.cells, "\n")
write.table(lambda_main.cells, file=paste0(cohort,".",main.cells,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summary as .csv file
main.cells_shiftwork.tableone <- as.data.frame(print(CreateTableOne(data=dataSW_main.cells[,-1],factorVars=c("mat.edu", "mat.smoking", "child.sex", "batch", "cell.names")),stringsAsFactors=FALSE))
write.csv(main.cells_shiftwork.tableone,file=paste0(cohort,".",main.cells,".summary.",analysis.date,".csv"))

print("Model 3: main.cells shif twork analysis completed") ## This completes the analysis of model 3. 


####################################################################### COMPLETE ANALYSES #############################################################

print("All analyses completed")








