###########################################################################################################################################################################################3#
# Purpose: Cross validation of a meta analysis result by leave-one(study)-out 
#          for CpGs with P values < 1E-07 in the meta-analysis
# Input: 1. Meta-analysis csv file including all cohorts, post processing
#        2. cohort results csv files
# Output:  1. Pdf file with leave one out plots for the selected CpGs  
#          2. Csv file with results of leave one out analysis for the selected CpGs 
#        
# adapted from original script by J. P. Jin
#############################################################################################################################################################################################
## ---------------------------------------------------------------
### ----- LEAVE-ONE-OUT ANALYSIS METAFOR [I-squared>50] ----------
## ---------------------------------------------------------------


# set working directory 
setwd("/home/060089/LOO_EWAS_SW/")

# load ggplot2
library(ggplot2)

###################
# Read in Meta-analysis result including all cohorts
###################
MA <- read.csv("MA_M3_MatShiftWork.25.02.2022.csv", header = T)



colnames(MA)[which(names(MA) == "MarkerName")] <- "PROBEID"
colnames(MA)[which(names(MA) == "Effect")] <- "BETA"
colnames(MA)[which(names(MA) == "StdErr")] <- "SE"
colnames(MA)[which(names(MA) == "P.value")] <- "P_VAL"


dim(MA)
#[1] 429959     32


###################
# select the relevant columns (probeID, P_VAL, UCSC_RefGene_Name
###################

MA<-MA[,c("PROBEID","BETA","SE","P_VAL","Direction","HetISq","CHR","UCSC_RefGene_Name")]
colnames(MA)[8] <- "gene"
# make subset of MA dataframe, only including the CpGs with small P value
#MA_subset <-subset(MA,MA$P_VAL < (0.05/nrow(MA)))  # is this significant to /by all tests?


MA_subset <-subset(MA,MA$P_VAL < 1.e-05)  # in this case, n=1

# check
dim(MA)
#[1] 16 8



###################
## Read in the cohort results 
###################


####
# full GI - overweight
####
ALSPAC<-read.csv("ALSPAC.main.cells2.results.20210806.csv", header = T)
summary(ALSPAC)

EAGER<-read.csv("EAGeR.main.cells.results.20210716.csv", header = T, strings = F)
summary(EAGER)

GENR<-read.csv("GenR.main.cells.results.20210603.csv", header = T, strings = F)
summary(GENR)

INMA <- read.csv("INMA.main.cells.results.20210729.csv", header = T)
summary(INMA)

POSEIDON <- read.csv("POSEIDON.main.cells.results.20210730.csv", header = T)
summary(POSEIDON)

MOBA1 <- read.csv("MoBa1.main.cells.csv", header = T)
summary(MOBA1)

MOBA2 <- read.csv("MOBA2.main.cells.results.20220106.csv", header = T)
summary(MOBA2)

######### subset? change colnames??? ####################

# make as subset of each cohort result file, only including cpgs of interest in the Meta-analysis subset     
ALSPAC <- ALSPAC[ALSPAC$PROBEID %in% MA$PROBEID, ]
head(ALSPAC)
summary(ALSPAC)
# Length:427530 

EAGER <- EAGER[EAGER$PROBEID %in% MA$PROBEID, ]
head(EAGER)
summary(EAGER)
#  Length:396361  

GENR <- GENR[GENR$PROBEID %in% MA$PROBEID, ]
head(GENR)
summary(GENR)
#  Length:415786  

INMA <- INMA[INMA$PROBEID %in% MA$PROBEID, ]
head(INMA)
summary(INMA)
#  Length:422721    

POSEIDON <- POSEIDON[POSEIDON$PROBEID %in% MA$PROBEID, ]
head(POSEIDON)
summary(POSEIDON)
#  Length:429959 

MOBA1 <- MOBA1[MOBA1$PROBEID %in% MA$PROBEID, ]
head(MOBA1)
summary(MOBA1)
# Length:429941 

MOBA2 <- MOBA2[MOBA2$PROBEID %in% MA$PROBEID, ]
head(MOBA2)
summary(MOBA2)
# Length:429859  




###################
## Merge the cohort and meta-analysis data frames  
###################


###################
# full GI - overweight
MA_cohorts <- merge(MA_subset, ALSPAC[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)
MA_cohorts <- merge(MA_cohorts, EAGER[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)
MA_cohorts <- merge(MA_cohorts,GENR[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)
MA_cohorts <- merge(MA_cohorts,INMA[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)
MA_cohorts <- merge(MA_cohorts,POSEIDON[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)
MA_cohorts <- merge(MA_cohorts,MOBA1[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)
MA_cohorts <- merge(MA_cohorts,MOBA2[,c(1:5)], by.x="PROBEID", by.y="PROBEID", sort=FALSE, all.x = T, all.y = F)



# rename colnames
colnames(MA_cohorts) <- c("MarkerName","Effect","StdErr","P.value","Direction","HetISq","CHR","gene",
							"ALSPAC_beta","ALSPAC_se","ALSPAC_p","ALSPAC_n",
							"EAGER_beta","EAGER_se","EAGER_p", "EAGER_n",
							"GENR_beta","GENR_se","GENR_p","GENR_n",
							"INMA_beta","INMA_se","INMA_p","INMA_n",
							"POSEIDON_beta","POSEIDON_se","POSEIDON_p","POSEIDON_n",
							"MOBA1_beta","MOBA1_se","MOBA1_p","MOBA1_n",
							"MOBA2_beta","MOBA2_se","MOBA2_p","MOBA2_n")

# check 
summary(MA_cohorts)
head(MA_cohorts)
dim(MA_cohorts)
#[1]  16 36





###################
# Write merged files (incl cohort results)
###################
write.csv(MA_cohorts, file="MA_cohorts_results_forestplot.csv")




######################################
######################################
# Create leave-one-out for Shift work EWAS
######################################
######################################


#Studies is a vector of study names
Studies <- as.vector(c("ALSPAC","EAGER","GENR","INMA","POSEIDON", "MOBA1", "MOBA2"))


#xy is a list of CpGs of interest
xy<-as.vector(as.character(MA_cohorts$MarkerName))

# don't know if we need to do something with this
## select the test set by HetISq
#xy<-as.vector(as.character(MA_cohorts$MarkerName[MA_cohorts$HetISq<50]))
#xy<-as.vector(as.character(MA_cohorts$MarkerName[MA_cohorts$HetISq>50]))


#M_C is the ewas results data frame, it contains the betas, SEs and N for each study plus the meta-analysis results
M_C <- MA_cohorts[MA_cohorts$MarkerName %in% xy, ]
M_C$MarkerName<-factor(M_C$MarkerName)		#for all CpGs M_C == MA_cohorts

#meta-analysis function (fixed effects)
FEmeta<-function(Z){
  require(metafor)
  Studies_B<-paste(unlist(Z),"_beta",sep="")
  Studies_S<-paste(unlist(Z),"_se",sep="")
  Studies_N<-paste(unlist(Z),"_n",sep="")
  M_C<-M_C[which(M_C$MarkerName %in% xy),]
  sites<-list(M_C$MarkerName)
  M_CB<-M_C[,c("MarkerName",Studies_B)]
  colnames(M_CB)<-c("PROBEID",unlist(Z))
  M_CS<-M_C[,c("MarkerName",Studies_S)]
  colnames(M_CS)<-c("PROBEID",unlist(Z))
  M_CN<-M_C[,c("MarkerName",Studies_N)]
  colnames(M_CN)<-c("PROBEID",unlist(Z))
  require(reshape)
  Betas<-melt(M_CB)
  names(Betas)<-c("PROBEID","Study","Beta")
  SEs<-melt(M_CS)
  Ns<-melt(M_CN)
  Data<-cbind(Betas,SEs[,"value"],Ns[,"value"])
  names(Data)<-c("PROBEID","Study","Betas","SE","Weight")
  Data$Study<-as.character(Data$Study)
  Data$Study<-as.factor(Data$Study)
  Data$PROBEID<-as.factor(Data$PROBEID)
  List<-split(Data,f=Data$PROBEID)
  List.res<-list(lapply(List,function(x) rma.uni(slab=x$Study, yi=x$Betas,sei=x$SE,method="FE",weighted=TRUE)),sites,xlab="Coefficient")
  List.res<-List.res[[1]]
}

#leave one out
List.res<-FEmeta(list(Studies))         ## list the result (fixed-effects model) of a CpG with test for heterogeneity, statistics that includes ci.lb & ci.ub  
Leave.res<-lapply(List.res,leave1out)   ## repeatedly fit the specified model, leaving out one observation/study at a time.
Leave.res.df<-lapply(Leave.res,function(X) as.data.frame(print(X))) #list each of which is saved as data frame

#take the full set result (original)
result<-vector()
##options(scipen=999)	#removing scientific notation; back to default: options(scipen=0); check the status: getOption("scipen") or options()$scipen
for(i in 1:length(names(List.res))) {
    stat <- unlist(sapply(List.res[i], function(x) c(x$b,x$se,x$zval,x$pval,x$ci.lb,x$ci.ub,x$QE,x$QEp,x$I2)))
    result <- cbind(result, stat)
}
rownames(result) <- c("estimate","se","zval","pval","ci.lb","ci.ub","Q","Qp","I2")
List.res.original <- t(data.frame(result))

#Combine original (full) results with L-O-O results 
#List.res.original = lapply(List.res, function(x)unlist(c(x$b[[1]],x[c("se","zval","pval","ci.lb","ci.ub","QE","QEp")])))
#prepare a file with correct annotation order and colour indicators for plotting
loo.processing<-
function(i,X.loo,X.original){
  res = rbind(X.original[i,],X.loo[[i]])
  rownames(res)[1] <- c("NONE")
  res$leftout<-factor(c("none",row.names(X.loo[[i]])),levels=c("none",sort(row.names(X.loo[[i]])),ordered=TRUE))
  res$colour <- c("black", rep("red",nrow(X.loo[[i]])))
  res
}

Leave.res.loo <-lapply(1:length(Leave.res.df),loo.processing,X.loo=Leave.res.df, X.original=List.res.original)
names(Leave.res.loo) <- names(Leave.res)

#write out a file with the columns of interest
options(digits=4,scipen=2)
final <- list()
for(i in 1:length(rownames(List.res.original))) {
    I2 <- as.numeric(c(List.res.original[i,9], rep("NA",length(rownames(Leave.res.loo[[i]]))-1)))
    final[[i]] <- cbind(Leave.res.loo[[i]][,c(1:4,8)],I2) 
}
names(final) <- names(Leave.res.loo)
final_de <- final
#capture.output(final, file="Loo_maternal_main.txt")	##list of which each result is saved as data frame format

##Leave-one-out forest plot 20%
for(i in 1:length(Leave.res)) {
    print(dim(data.frame(Leave.res[i]))[1])
}
 
draw.loo.plot<- function(list.of.dfs,i){ 
   df <- list.of.dfs[[i]] 
   name.df <- names(list.of.dfs)[[i]] 
   ggplot(df,aes(x=leftout,y=estimate,colour=colour))+ 
geom_point()+ 
geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),width=0.5)+ 
scale_color_manual(values=c("black","red"))+ 
geom_hline(yintercept=0,linetype="dashed")+ 
geom_rect(aes(xmin=0, xmax=Inf, ymin=(estimate[1]-(estimate[1]*0.20)), ymax=(estimate[1]+(estimate[1]*0.20))), linetype="blank", col="darkblue", fill="darkblue", alpha = .1)+
theme_bw() + theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1),panel.grid.major.x = element_blank()) + 
xlab("Left out study (highlight: +/- 20% of M.A. beta)")+ylab("Effect Estimate (95% CI)")+ 
ggtitle(name.df) + theme(plot.title = element_text(hjust = 0.5)) 
} 

pdf("LOO.plots.fullGI.overweight.pdf",width=6,height=4) 
par(mfrow=c(2,2))
lapply(1:length(Leave.res.loo),draw.loo.plot,list.of.dfs=Leave.res.loo) 

dev.off() 
 
 
 
##Leave-one-out forest plot 50%
for(i in 1:length(Leave.res)) {
    print(dim(data.frame(Leave.res[i]))[1])
}
 
draw.loo.plot<- function(list.of.dfs,i){ 
   df <- list.of.dfs[[i]] 
   name.df <- names(list.of.dfs)[[i]] 
   ggplot(df,aes(x=leftout,y=estimate,colour=colour))+ 
geom_point()+ 
geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),width=0.5)+ 
scale_color_manual(values=c("black","red"))+ 
geom_hline(yintercept=0,linetype="dashed")+ 
geom_rect(aes(xmin=0, xmax=Inf, ymin=(estimate[1]-(estimate[1]*0.50)), ymax=(estimate[1]+(estimate[1]*0.50))), linetype="blank", col="darkblue", fill="darkblue", alpha = .1)+
theme_bw() + theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1),panel.grid.major.x = element_blank()) + 
xlab("Left out study (highlight: +/- 50% of M.A. beta)")+ylab("Effect Estimate (95% CI)")+ 
ggtitle(name.df) + theme(plot.title = element_text(hjust = 0.5)) 
} 

pdf("LOO.plots.fullGI.overweight.50pct.pdf",width=6,height=4) 
par(mfrow=c(2,2))
lapply(1:length(Leave.res.loo),draw.loo.plot,list.of.dfs=Leave.res.loo) 

dev.off() 
