module load apps/gcc/R/3.6.1 
qrsh -l short -V -cwd -pe smp.pe 8 R --vanilla --interactive
#load biolink libraries
####This did not work so I used the GDCRNATools for acquiring read counts
###This will give me data to look at the miRNA isoform quantification for the 437 patients

load("/mnt/iusers01/cw01/mqbpkmk3/TCGA_BLCA.RData")####2588 miRNAs inclusing both arms and 437 patients
dataRPM = apply(RNA,2,function(x){(x/sum(x,na.rm=TRUE))*1000000})
dataRPM<-log2(dataRPM+1)
seed_genes<-read.table("/mnt/iusers01/cw01/mqbpkmk3/miRNA_seed.csv",sep=",",stringsAsFactors=FALSE)
###will have to update it and make a new one as have to add 1 and 3 isoforms to it...
###need to look at the full list of miRNAs first...
###write.table(rownames(dataRPM),file="/mnt/iusers01/cw01/mqbpkmk3/miRNA_names.csv",sep=",")
###scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/miRNA_names.csv "/Users/mkhan/Documents/fifth year"
###There are a total of 417 patients, which ones do I remove...
###Look at the development of the miRNA signature...
### change 320
seed_genes[seed_genes$V1=="hsa-miR-146b",]<- "hsa-miR-146b-5p"
dataRPM<-dataRPM[as.character(unique(seed_genes$V1)),]
#####remove the normal samples
tumour<-which(grepl("-01",colnames(dataRPM)))
dataRPM<-dataRPM[,tumour]####417 patients
dataRPM<-dataRPM[,-which(duplicated(colnames(dataRPM)))]
###remove those miRNAs whose AVERAGE less than 1
dataRPM<-dataRPM[rowMeans(dataRPM)>1,]
namescounts<-sapply(colnames(dataRPM),function(x){
  NAME<-strsplit(as.character(x), fixed=T, "-")[[1]][1:3]
  NAME<-as.vector(NAME)
  paste(c(NAME[1],NAME[2],NAME[3]),collapse = "-")})
colnames(dataRPM) <-namescounts
dataRPM<-as.data.frame(dataRPM)

###load libraries
library(caret)
library(Boruta)
library("randomForest")
library("varSelRF")
library("survival")
library("survminer")


########Divide into train and test dataset using the categories.....
winter_manual<-get(load("Category_winter.RData"))
winter_scores<-get(load("scores_winter.RData"))
expression_winter<-merge(t(dataRPM),winter_manual,by="row.names")
expression_winter<-na.omit(expression_winter)
rownames(expression_winter)<-expression_winter$Row.names
colnames(expression_winter)[ncol(expression_winter)]<-"median_stratification"

###Do Boruta with train and test split!!!
set.seed(123)
trainIndex <- createDataPartition(factor(expression_winter$median_stratification),
                                  p = 0.70, list = FALSE)

train_data <- expression_winter[trainIndex, ]
test_data <-expression_winter[-trainIndex, ]
#save(train_data,file="/mnt/iusers01/cw01/mqbpkmk3/train_data.RData")
#save(test_data,file="/mnt/iusers01/cw01/mqbpkmk3/test_data.RData")
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/train_data.RData "/Users/mkhan/Documents/fifth year"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/test_data.RData "/Users/mkhan/Documents/fifth year"

train_data<-train_data[,2:61]
set.seed(123)
algo<-Boruta(median_stratification ~ .,data=train_data,doTrace=2,maxRuns=1000)###19 miRNAs, 38, 3
getConfirmedFormula(algo)###19 miRNAs
#pdf(file="/mnt/iusers01/cw01/mqbpkmk3/boruta_manual_winter_all_miRNAs.pdf")
#plot(algo,las=2,cex.axis=0.5,xlab="")
#dev.off()
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/boruta_manual_winter_all_miRNAs.pdf "/Users/mkhan/Documents/fifth year"
boruta.df <- attStats(algo)
imp_miRNA<-rownames(boruta.df[boruta.df$decision=="Confirmed",])
imp_miRNA<-gsub(pattern = "`", replacement = "", x = imp_miRNA)


###Also check same 19 using varSELRF
set.seed(123)
rf.vsl<-varSelRF(train_data[,imp_miRNA],train_data[,"median_stratification"], c.sd = 0, mtryFactor = 1, ntree = 5000,
         ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2,
         whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE,
         returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE)

####survival analysis
tr_survival<-read.table("/mnt/iusers01/cw01/mqbpkmk3/clinical_bladder_balance_edited.csv",
                        header=T,
                        sep=",",
                        stringsAsFactors = FALSE)
#overall survival values in terms of 0 and 1
OS<-sapply(tr_survival$Overall.Survival.Status,
           function(x){
             if(x=="DECEASED"){
               y=1
             }else if(x=="LIVING"){
               y=0
             }})
tr_survival$OS<-OS
#censoring dataset at 5 years for overall survival
censorship<-60
tr_survival_event_censorshipo  <- ifelse(  tr_survival$Overall.Survival..Months. <= censorship & tr_survival$OS == 1, 1 ,0)
tr_survival_time_censorshipo   <- ifelse(  tr_survival_event_censorshipo == 0 & tr_survival$Overall.Survival..Months. >= censorship, censorship , tr_survival$Overall.Survival..Months. )   
tr_survival                   <- cbind( tr_survival, tr_survival_time_censorshipo, tr_survival_event_censorshipo )
colnames( tr_survival )[ (ncol(tr_survival)-1):ncol(tr_survival) ] <- c( "censored_otime", "censored_ostatus" )
#remove duplicate patient IDs
tr_survival<-tr_survival[!duplicated(tr_survival$Patient.ID), ]
#setting rownames of clinical data
rownames(tr_survival)<-tr_survival$Patient.ID
#censor dataset at 5 years for progression free survival
tr_survival<-tr_survival[!is.na(tr_survival$Progression.Free.Status),]
tr_survival<-tr_survival[!is.na(tr_survival$Progress.Free.Survival..Months.),]
#censoring dataset at 5 years for progression free survival
#Progression free survival in terms of 0 and 1
PFS<-sapply(tr_survival$Progression.Free.Status,
           function(x){
             if(x=="PROGRESSION"){
               y=1
             }else if(x=="CENSORED"){
               y=0
             }})
tr_survival$PFS<-PFS
#censor dataset at 5 years for progression free survival
censorship<-60
tr_survival_event_censorshipo  <- ifelse(  tr_survival$Progress.Free.Survival..Months. <= censorship & tr_survival$PFS == 1, 1 ,0)
tr_survival_time_censorshipo   <- ifelse(  tr_survival_event_censorshipo == 0 & tr_survival$Progress.Free.Survival..Months. >= censorship, censorship , tr_survival$Progress.Free.Survival..Months. )   
tr_survival                   <- cbind( tr_survival, tr_survival_time_censorshipo, tr_survival_event_censorshipo )
colnames( tr_survival )[ (ncol(tr_survival)-1):ncol(tr_survival) ] <- c( "censored_ptime", "censored_pstatus" )

####Spearman correlations with winter scores
winter<-merge((train_data[,imp_miRNA]),winter_scores,by="row.names")
rownames(winter)<-winter$Row.names
colnames(winter)[ncol(winter)]<-"winter_scores" ####405 patients
winter<-na.omit(winter)
correlation<- sapply(imp_miRNA,function(x){
  cancer<-winter[, c(x,"winter_scores")]
  significance <- cor.test(cancer[,x],cancer[,"winter_scores"],use="everything", method=c("spearman") )
  pvalue<-as.numeric(unlist(significance[3]))
  cor<-as.numeric(unlist(significance[4]))
cor_winter<-cbind(pvalue,cor)})
cor<- (as.data.frame(t(correlation)))
colnames(cor)[1]<-"p-value_winter"
colnames(cor)[2]<-"estimate_winter"
FDRcorrection_winter<- p.adjust(cor[,1], method="fdr" )

pearson<- as.data.frame(cbind(cor,FDRcorrection_winter))
pearson<-pearson[order(-pearson[,2]),]
x<-rownames(pearson[pearson$FDRcorrection_winter<0.05,])
hyp<-rownames(pearson[pearson$FDRcorrection_winter<0.05&pearson$estimate_winter>0,])
nor<-rownames(pearson[pearson$FDRcorrection_winter<0.05&pearson$estimate_winter<0,])
x<-c(hyp,nor)

expression <- as.data.frame(t(train_data[,x ]))
#look at using overall survival
pts.expr<-sapply(colnames(expression),function(y){
hypoxia<-mean((expression[hyp,y]))
normoxia<-mean((expression[nor,y]))
p<-hypoxia-normoxia
})
medianexpression <- quantile(pts.expr,0.50)
pts.category <- ifelse(pts.expr > medianexpression, "Hypoxia", "Normoxia")
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,
                             row.names = colnames(expression),
                             stringsAsFactors = FALSE)
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,stringsAsFactors = FALSE)
 rownames(combinationcox)<-colnames(expression)
 combinationcox <- cbind(combinationcox,
                        tr_survival[colnames(expression), c("censored_ostatus","censored_otime","Diagnosis.Age","Sex","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","censored_ptime", "censored_pstatus")])
###look at using progression free survival

combinationcox$Category<-factor(combinationcox$Category,levels=c("Normoxia","Hypoxia"))
#Overall survival
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Category, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Score, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Sex, data = combinationcox)
summary(res.cox)

stage<-combinationcox[!is.na(combinationcox$American.Joint.Committee.on.Cancer.Tumor.Stage.Code),]
stagecategory<-sapply(stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code,function(x){
  if(x=="T0"|x=="T1"|x=="TX"|x=="T2"|x=="T2a"|x=="T2b"){y=2}else{
    y=3}})

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory, data = stage)
summary(res.cox) 
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory+Diagnosis.Age+Category, data = stage)
summary(res.cox)
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Category, data = combinationcox)
summary(res.cox)

#####multivariate analysis
#age stage sex
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Sex, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory, data = stage)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory+Category, data = stage)
summary(res.cox)
  
####test dataset
expression1 <- as.data.frame(t(test_data[,x]))
#look at using overall survival
#pts.expr <- apply(expression1,2,median)
pts.expr<-sapply(colnames(expression1),function(y){
	hypoxia<-mean((expression1[hyp,y]))
	normoxia<-mean((expression1[nor,y]))
	p<-hypoxia-normoxia
})

medianexpression <- quantile(pts.expr,0.50)
pts.category <- ifelse(pts.expr > medianexpression, "Hypoxia", "Normoxia")
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,
                             row.names = colnames(expression1),
                             stringsAsFactors = FALSE)
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,stringsAsFactors = FALSE)
 rownames(combinationcox)<-colnames(expression1)
 combinationcox <- cbind(combinationcox,
                        tr_survival[colnames(expression1), c("censored_ostatus","censored_otime","Diagnosis.Age","Sex","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","censored_ptime", "censored_pstatus")])
###look at using progression free survival

combinationcox$Category<-factor(combinationcox$Category,levels=c("Normoxia","Hypoxia"))
#Overall survival
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Category, data = combinationcox)
summary(res.cox)
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Sex, data = combinationcox)
summary(res.cox)

stage<-combinationcox[!is.na(combinationcox$American.Joint.Committee.on.Cancer.Tumor.Stage.Code),]
stagecategory<-sapply(stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code,function(x){
  if(x=="T0"|x=="T1"|x=="TX"|x=="T2"|x=="T2a"|x=="T2b"){y=2}else{
    y=3}})

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory, data = stage)
summary(res.cox) 

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~Diagnosis.Age+Category, data = stage)
summary(res.cox) 

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Category, data = combinationcox)
summary(res.cox)


#####multivariate analysis
#age stage sex
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)#0.4 1.01[0.98-1.04]

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Sex, data = combinationcox)
summary(res.cox)#0.5 1.31[0.62-2.75]

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory, data = stage)
summary(res.cox)#0.9, 0.94 [0.37-2.43]

####Km curves to be plotted on local computer

###random forest to measure Gini importance

rf.df<-train_data[,c(rownames(dataRPM),"median_stratification")]
colnames(rf.df)[1:ncol(rf.df)-1]<-gsub("-",".",colnames(rf.df[1:ncol(rf.df)-1]))
set.seed(123)
rf<- randomForest(median_stratification ~.,data=rf.df)
pdf(file="/mnt/iusers01/cw01/mqbpkmk3/varimplot_all.pdf")
varImpPlot(rf)
dev.off()
scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/varimplot_all.pdf "/Users/mkhan/Documents/fifth year"


####Test the NMIBC signature in the Ochoa cohort
NMBC<-c("hsa-miR-210-3p","hsa-miR-193b-3p","hsa-miR-125a-3p","hsa-miR-708-5p","hsa-miR-517a-3p","hsa-miR-145-5p")
dataRPM = apply(RNA,2,function(x){(x/sum(x,na.rm=TRUE))*1000000})
dataRPM<-log2(dataRPM+1)
tumour<-which(grepl("-01",colnames(dataRPM)))
dataRPM<-dataRPM[,tumour]####417 patients
dataRPM<-dataRPM[,-which(duplicated(colnames(dataRPM)))]
namescounts<-sapply(colnames(dataRPM),function(x){
  NAME<-strsplit(as.character(x), fixed=T, "-")[[1]][1:3]
  NAME<-as.vector(NAME)
  paste(c(NAME[1],NAME[2],NAME[3]),collapse = "-")})
colnames(dataRPM) <-namescounts
dataRPM<-as.data.frame(dataRPM)

signature_NMBC<-dataRPM[NMBC,]
scores<-apply(signature_NMBC,2,mean)
categories<-ifelse(scores>quantile(scores,0.75),"Hypoxia","Normoxia")
survival_scores<-merge(tr_survival,categories,by="row.names")
rownames(survival_scores)<-survival_scores$Row.names

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ y, data = survival_scores)
summary(res.cox)
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ y, data = survival_scores)
summary(res.cox)


####Look at the signature for overall survival in the whole TCGA cohort
dataRPM_imp<-dataRPM[x,c(rownames(train_data),rownames(test_data))]

expression <- as.data.frame((dataRPM_imp))
#look at using overall survival
pts.expr<-sapply(colnames(expression),function(y){
hypoxia<-mean((expression[hyp,y]))
normoxia<-mean((expression[nor,y]))
p<-hypoxia-normoxia
})
medianexpression <- quantile(pts.expr,0.50)
pts.category <- ifelse(pts.expr > medianexpression, "Hypoxia", "Normoxia")
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,
                             row.names = colnames(expression),
                             stringsAsFactors = FALSE)
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,stringsAsFactors = FALSE)
 rownames(combinationcox)<-colnames(expression)
 combinationcox <- cbind(combinationcox,
                        tr_survival[colnames(expression), c("censored_ostatus","censored_otime","Diagnosis.Age","Sex","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","censored_ptime", "censored_pstatus")])
###look at using progression free survival

combinationcox$Category<-factor(combinationcox$Category,levels=c("Normoxia","Hypoxia"))
#Overall survival
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Category, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Score, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Sex, data = combinationcox)
summary(res.cox)

stage<-combinationcox[!is.na(combinationcox$American.Joint.Committee.on.Cancer.Tumor.Stage.Code),]
stagecategory<-sapply(stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code,function(x){
  if(x=="T0"|x=="T1"|x=="TX"|x=="T2"|x=="T2a"|x=="T2b"){y=2}else{
    y=3}})

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory, data = stage)
summary(res.cox) 
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory+Diagnosis.Age+Category, data = stage)
summary(res.cox)
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Category, data = combinationcox)
summary(res.cox)

#####multivariate analysis
#age stage sex
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Sex, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory, data = stage)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory+Category, data = stage)
summary(res.cox)


########
module load apps/gcc/R/3.6.1 
qrsh -l short -V -cwd -pe smp.pe 8 R --vanilla --interactive
#load biolink libraries
####This did not work so I used the GDCRNATools for acquiring read counts
###This will give me data to look at the miRNA isoform quantification for the 437 patients
setwd("/Users/mkhan/Documents/miRNA_final_work")
load("TCGA_BLCA.RData")####2588 miRNAs inclusing both arms and 437 patients
dataRPM = apply(RNA,2,function(x){(x/sum(x,na.rm=TRUE))*1000000})
dataRPM<-log2(dataRPM+1)
seed_genes<-read.table("miRNA_seed.csv",sep=",",stringsAsFactors=FALSE)
###will have to update it and make a new one as have to add 1 and 3 isoforms to it...
###need to look at the full list of miRNAs first...
###write.table(rownames(dataRPM),file="/mnt/iusers01/cw01/mqbpkmk3/miRNA_names.csv",sep=",")
###scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/miRNA_names.csv "/Users/mkhan/Documents/fifth year"
###There are a total of 417 patients, which ones do I remove...
###Look at the development of the miRNA signature...
### change 320
seed_genes[seed_genes$V1=="hsa-miR-146b",]<- "hsa-miR-146b-5p"
dataRPM<-dataRPM[as.character(unique(seed_genes$V1)),]
#####remove the normal samples
tumour<-which(grepl("-01",colnames(dataRPM)))
dataRPM<-dataRPM[,tumour]####417 patients
dataRPM<-dataRPM[,-which(duplicated(colnames(dataRPM)))]
###remove those miRNAs whose AVERAGE less than 1
dataRPM<-dataRPM[rowMeans(dataRPM)>1,]
namescounts<-sapply(colnames(dataRPM),function(x){
  NAME<-strsplit(as.character(x), fixed=T, "-")[[1]][1:3]
  NAME<-as.vector(NAME)
  paste(c(NAME[1],NAME[2],NAME[3]),collapse = "-")})
colnames(dataRPM) <-namescounts
dataRPM<-as.data.frame(dataRPM)

###load libraries
library(caret)
library(Boruta)
library("randomForest")
library("varSelRF")
library("survival")
library("survminer")


########Divide into train and test dataset using the categories.....
winter_manual<-get(load("Category_winter.RData"))
winter_scores<-get(load("scores_winter.RData"))
expression_winter<-merge(t(dataRPM),winter_manual,by="row.names")
expression_winter<-na.omit(expression_winter)
rownames(expression_winter)<-expression_winter$Row.names
colnames(expression_winter)[ncol(expression_winter)]<-"median_stratification"

###Do Boruta with train and test split!!!
set.seed(123)
trainIndex <- createDataPartition(factor(expression_winter$median_stratification),
                                  p = 0.70, list = FALSE)

train_data <- expression_winter[trainIndex, ]
test_data <-expression_winter[-trainIndex, ]
#save(train_data,file="/mnt/iusers01/cw01/mqbpkmk3/train_data.RData")
#save(test_data,file="/mnt/iusers01/cw01/mqbpkmk3/test_data.RData")
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/train_data.RData "/Users/mkhan/Documents/fifth year"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/test_data.RData "/Users/mkhan/Documents/fifth year"

train_data<-train_data[,2:61]
set.seed(123)
algo<-Boruta(median_stratification ~ .,data=train_data,doTrace=2,maxRuns=1000)###19 miRNAs, 38, 3
getConfirmedFormula(algo)###19 miRNAs
#pdf(file="/mnt/iusers01/cw01/mqbpkmk3/boruta_manual_winter_all_miRNAs.pdf")
#plot(algo,las=2,cex.axis=0.5,xlab="")
#dev.off()
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/boruta_manual_winter_all_miRNAs.pdf "/Users/mkhan/Documents/fifth year"
boruta.df <- attStats(algo)
imp_miRNA<-rownames(boruta.df[boruta.df$decision=="Confirmed",])
imp_miRNA<-gsub(pattern = "`", replacement = "", x = imp_miRNA)


###Also check same 19 using varSELRF
set.seed(123)
rf.vsl<-varSelRF(train_data[,imp_miRNA],train_data[,"median_stratification"], c.sd = 0, mtryFactor = 1, ntree = 5000,
         ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2,
         whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE,
         returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE)

####survival analysis
tr_survival<-read.table("clinical_bladder_balance_edited.csv",
                        header=T,
                        sep=",",
                        stringsAsFactors = FALSE)
#overall survival values in terms of 0 and 1
OS<-sapply(tr_survival$Overall.Survival.Status,
           function(x){
             if(x=="DECEASED"){
               y=1
             }else if(x=="LIVING"){
               y=0
             }})
tr_survival$OS<-OS
#censoring dataset at 5 years for overall survival
censorship<-60
tr_survival_event_censorshipo  <- ifelse(  tr_survival$Overall.Survival..Months. <= censorship & tr_survival$OS == 1, 1 ,0)
tr_survival_time_censorshipo   <- ifelse(  tr_survival_event_censorshipo == 0 & tr_survival$Overall.Survival..Months. >= censorship, censorship , tr_survival$Overall.Survival..Months. )   
tr_survival                   <- cbind( tr_survival, tr_survival_time_censorshipo, tr_survival_event_censorshipo )
colnames( tr_survival )[ (ncol(tr_survival)-1):ncol(tr_survival) ] <- c( "censored_otime", "censored_ostatus" )
#remove duplicate patient IDs
tr_survival<-tr_survival[!duplicated(tr_survival$Patient.ID), ]
#setting rownames of clinical data
rownames(tr_survival)<-tr_survival$Patient.ID
#censor dataset at 5 years for progression free survival
tr_survival<-tr_survival[!is.na(tr_survival$Progression.Free.Status),]
tr_survival<-tr_survival[!is.na(tr_survival$Progress.Free.Survival..Months.),]
#censoring dataset at 5 years for progression free survival
#Progression free survival in terms of 0 and 1
PFS<-sapply(tr_survival$Progression.Free.Status,
           function(x){
             if(x=="PROGRESSION"){
               y=1
             }else if(x=="CENSORED"){
               y=0
             }})
tr_survival$PFS<-PFS
#censor dataset at 5 years for progression free survival
censorship<-60
tr_survival_event_censorshipo  <- ifelse(  tr_survival$Progress.Free.Survival..Months. <= censorship & tr_survival$PFS == 1, 1 ,0)
tr_survival_time_censorshipo   <- ifelse(  tr_survival_event_censorshipo == 0 & tr_survival$Progress.Free.Survival..Months. >= censorship, censorship , tr_survival$Progress.Free.Survival..Months. )   
tr_survival                   <- cbind( tr_survival, tr_survival_time_censorshipo, tr_survival_event_censorshipo )
colnames( tr_survival )[ (ncol(tr_survival)-1):ncol(tr_survival) ] <- c( "censored_ptime", "censored_pstatus" )

####Spearman correlations with winter scores
winter<-merge((train_data[,imp_miRNA]),winter_scores,by="row.names")
rownames(winter)<-winter$Row.names
colnames(winter)[ncol(winter)]<-"winter_scores" ####405 patients
winter<-na.omit(winter)
correlation<- sapply(imp_miRNA,function(x){
  cancer<-winter[, c(x,"winter_scores")]
  significance <- cor.test(cancer[,x],cancer[,"winter_scores"],use="everything", method=c("spearman") )
  pvalue<-as.numeric(unlist(significance[3]))
  cor<-as.numeric(unlist(significance[4]))
cor_winter<-cbind(pvalue,cor)})
cor<- (as.data.frame(t(correlation)))
colnames(cor)[1]<-"p-value_winter"
colnames(cor)[2]<-"estimate_winter"
FDRcorrection_winter<- p.adjust(cor[,1], method="fdr" )

pearson<- as.data.frame(cbind(cor,FDRcorrection_winter))
pearson<-pearson[order(-pearson[,2]),]
x<-rownames(pearson[pearson$FDRcorrection_winter<0.05,])
hyp<-rownames(pearson[pearson$FDRcorrection_winter<0.05&pearson$estimate_winter>0,])
nor<-rownames(pearson[pearson$FDRcorrection_winter<0.05&pearson$estimate_winter<0,])
x<-c(hyp,nor)

expression <- as.data.frame(t(train_data[,x ]))
#look at using overall survival
pts.expr<-sapply(colnames(expression),function(y){
hypoxia<-mean((expression[hyp,y]))
normoxia<-mean((expression[nor,y]))
p<-hypoxia-normoxia
})
medianexpression <- quantile(pts.expr,0.50)
pts.category <- ifelse(pts.expr > medianexpression, "Hypoxia", "Normoxia")
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,
                             row.names = colnames(expression),
                             stringsAsFactors = FALSE)
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,stringsAsFactors = FALSE)
 rownames(combinationcox)<-colnames(expression)
 combinationcox <- cbind(combinationcox,
                        tr_survival[colnames(expression), c("censored_ostatus","censored_otime","Diagnosis.Age","Sex","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","censored_ptime", "censored_pstatus")])
###look at using progression free survival

combinationcox$Category<-factor(combinationcox$Category,levels=c("Normoxia","Hypoxia"))
#Overall survival
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Category, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Score, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Sex, data = combinationcox)
summary(res.cox)

stage<-combinationcox[!is.na(combinationcox$American.Joint.Committee.on.Cancer.Tumor.Stage.Code),]
stagecategory<-sapply(stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code,function(x){
  if(x=="T0"|x=="T1"|x=="TX"|x=="T2"|x=="T2a"|x=="T2b"){y=2}else{
    y=3}})

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory, data = stage)
summary(res.cox) 
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory+Diagnosis.Age+Category, data = stage)
summary(res.cox)
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Category, data = combinationcox)
summary(res.cox)

#####multivariate analysis
#age stage sex
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Sex, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory, data = stage)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory+Category, data = stage)
summary(res.cox)
  
####test dataset
expression1 <- as.data.frame(t(test_data[,x]))
#look at using overall survival
#pts.expr <- apply(expression1,2,median)
pts.expr<-sapply(colnames(expression1),function(y){
	hypoxia<-mean((expression1[hyp,y]))
	normoxia<-mean((expression1[nor,y]))
	p<-hypoxia-normoxia
})

medianexpression <- quantile(pts.expr,0.50)
pts.category <- ifelse(pts.expr > medianexpression, "Hypoxia", "Normoxia")
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,
                             row.names = colnames(expression1),
                             stringsAsFactors = FALSE)
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,stringsAsFactors = FALSE)
 rownames(combinationcox)<-colnames(expression1)
 combinationcox <- cbind(combinationcox,
                        tr_survival[colnames(expression1), c("censored_ostatus","censored_otime","Diagnosis.Age","Sex","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","censored_ptime", "censored_pstatus")])
###look at using progression free survival

combinationcox$Category<-factor(combinationcox$Category,levels=c("Normoxia","Hypoxia"))
#Overall survival
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Category, data = combinationcox)
summary(res.cox)
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Sex, data = combinationcox)
summary(res.cox)

stage<-combinationcox[!is.na(combinationcox$American.Joint.Committee.on.Cancer.Tumor.Stage.Code),]
stagecategory<-sapply(stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code,function(x){
  if(x=="T0"|x=="T1"|x=="TX"|x=="T2"|x=="T2a"|x=="T2b"){y=2}else{
    y=3}})

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory, data = stage)
summary(res.cox) 

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~Diagnosis.Age+Category, data = stage)
summary(res.cox) 

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Category, data = combinationcox)
summary(res.cox)


#####multivariate analysis
#age stage sex
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)#0.4 1.01[0.98-1.04]

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Sex, data = combinationcox)
summary(res.cox)#0.5 1.31[0.62-2.75]

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory, data = stage)
summary(res.cox)#0.9, 0.94 [0.37-2.43]

####Km curves to be plotted on local computer

###random forest to measure Gini importance

rf.df<-train_data[,c(rownames(dataRPM),"median_stratification")]
colnames(rf.df)[1:ncol(rf.df)-1]<-gsub("-",".",colnames(rf.df[1:ncol(rf.df)-1]))
set.seed(123)
rf<- randomForest(median_stratification ~.,data=rf.df)
pdf(file="/mnt/iusers01/cw01/mqbpkmk3/varimplot_all.pdf")
varImpPlot(rf)
dev.off()
scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/varimplot_all.pdf "/Users/mkhan/Documents/fifth year"


####Test the NMIBC signature in the Ochoa cohort
NMBC<-c("hsa-miR-210-3p","hsa-miR-193b-3p","hsa-miR-125a-3p","hsa-miR-708-5p","hsa-miR-517a-3p","hsa-miR-145-5p")
dataRPM = apply(RNA,2,function(x){(x/sum(x,na.rm=TRUE))*1000000})
dataRPM<-log2(dataRPM+1)
tumour<-which(grepl("-01",colnames(dataRPM)))
dataRPM<-dataRPM[,tumour]####417 patients
dataRPM<-dataRPM[,-which(duplicated(colnames(dataRPM)))]
namescounts<-sapply(colnames(dataRPM),function(x){
  NAME<-strsplit(as.character(x), fixed=T, "-")[[1]][1:3]
  NAME<-as.vector(NAME)
  paste(c(NAME[1],NAME[2],NAME[3]),collapse = "-")})
colnames(dataRPM) <-namescounts
dataRPM<-as.data.frame(dataRPM)

signature_NMBC<-dataRPM[NMBC,]
scores<-apply(signature_NMBC,2,mean)
categories<-ifelse(scores>quantile(scores,0.75),"Hypoxia","Normoxia")
survival_scores<-merge(tr_survival,categories,by="row.names")
rownames(survival_scores)<-survival_scores$Row.names

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ y, data = survival_scores)
summary(res.cox)
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ y, data = survival_scores)
summary(res.cox)


####Look at the signature for overall survival in the whole TCGA cohort
dataRPM_imp<-dataRPM[x,c(rownames(train_data),rownames(test_data))]

expression <- as.data.frame((dataRPM_imp))
#look at using overall survival
pts.expr<-sapply(colnames(expression),function(y){
hypoxia<-mean((expression[hyp,y]))
normoxia<-mean((expression[nor,y]))
p<-hypoxia-normoxia
})
medianexpression <- quantile(pts.expr,0.50)
pts.category <- ifelse(pts.expr > medianexpression, "Hypoxia", "Normoxia")
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,
                             row.names = colnames(expression),
                             stringsAsFactors = FALSE)
combinationcox <- data.frame("Score" = pts.expr,
                             "Category" = pts.category,stringsAsFactors = FALSE)
 rownames(combinationcox)<-colnames(expression)
 combinationcox <- cbind(combinationcox,
                        tr_survival[colnames(expression), c("censored_ostatus","censored_otime","Diagnosis.Age","Sex","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","censored_ptime", "censored_pstatus")])
###look at using progression free survival

combinationcox$Category<-factor(combinationcox$Category,levels=c("Normoxia","Hypoxia"))
#Overall survival
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Category, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Score, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~ Sex, data = combinationcox)
summary(res.cox)

stage<-combinationcox[!is.na(combinationcox$American.Joint.Committee.on.Cancer.Tumor.Stage.Code),]
stagecategory<-sapply(stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code,function(x){
  if(x=="T0"|x=="T1"|x=="TX"|x=="T2"|x=="T2a"|x=="T2b"){y=2}else{
    y=3}})

res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory, data = stage)
summary(res.cox) 
res.cox <- coxph(Surv(censored_otime,censored_ostatus) ~stagecategory+Diagnosis.Age+Category, data = stage)
summary(res.cox)
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Category, data = combinationcox)
summary(res.cox)

#####multivariate analysis
#age stage sex
res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Diagnosis.Age, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~ Sex, data = combinationcox)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory, data = stage)
summary(res.cox)

res.cox <- coxph(Surv(censored_ptime,censored_pstatus) ~stagecategory+Category, data = stage)
summary(res.cox)





























