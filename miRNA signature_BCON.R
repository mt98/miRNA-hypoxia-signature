#analyisng nanostring data for the miRNA signature for BCON
setwd("/Users/mkhan/Documents/fifth year")
#load library
library("survival")
library("survminer")
##get the survival data for BCON
survivalanalysis<-read.table("/Users/mkhan/Documents/ third year/bladder data/bladder infofiles/bladder cohort.csv",sep=",",quote="",header=T,stringsAsFactors = FALSE)
rownames(survivalanalysis)<-survivalanalysis$X
###look at the characteristics of the BCON Data for the methods for all 192 patients
patients<-read.table("/Users/mkhan/Documents/fifth year/miRNA hypoxia signature/Nanostring_data/sampleno_patient no_BCONmiRNA.csv",sep=",",quote="",header=T,stringsAsFactors = FALSE)
selected<-survivalanalysis[as.character(patients$Patient.no),]
selected_RT<-selected[selected$treatment=="RT",]
###
selected_RT_CON<-selected[selected$treatment=="RT+CON",]

###read the processed BCON Data

BCONdata<-read.csv("miRNA hypoxia signature/Nanostring_data/log2 normalised without flags.csv",sep=",",header=T)
rownames(BCONdata)<-BCONdata$Probe.Name

hyp<-c("hsa-miR-193b-3p" , "hsa-miR-21-5p"  ,"hsa-miR-210-3p" , "hsa-miR-221-3p" , 
       "hsa-miR-224-5p","hsa-miR-27a-3p","hsa-miR-455-5p") 
nor<-c( "hsa-miR-190a-5p", "hsa-miR-191-5p", "hsa-miR-28-5p" ,  "hsa-miR-30b-5p","hsa-miR-491-5p" ,
        "hsa-miR-182-5p","hsa-miR-93-5p")
sig<-c(hyp,nor)
signature<-BCONdata[sig,-c(1:7)]
scores<-sapply(colnames(signature),function(y){
  hypoxia<-mean((signature[hyp,y]))
  normoxia<-mean((signature[nor,y]))
  p<-hypoxia-normoxia
})
names(scores)<-sub('X', '', names(scores))
categories<-ifelse(scores>quantile(scores,0.75),"Hypoxia","Normoxia")
##get the survival data for BCON
survivalanalysis<-get(load("/Users/mkhan/Documents/fourth year /miRNAs/BLCA signature/BCON_survival.RData"))
#censor the survival analysis for os and lpfs
censorship<-60
BCON_survival_event_censorship  <- ifelse(survivalanalysis$os.time <= censorship & survivalanalysis$os.status == 1 , 1 ,0  )
BCON_survival_time_censorship   <- ifelse( BCON_survival_event_censorship == 0 & survivalanalysis$os.time >= censorship , censorship ,survivalanalysis$os.time )   
finalsurvivalanalysis       <- cbind( survivalanalysis, BCON_survival_event_censorship, BCON_survival_time_censorship )
colnames(finalsurvivalanalysis)[ (ncol(finalsurvivalanalysis)-1):ncol(finalsurvivalanalysis) ] <- c("c_event","c_event_time")
BCON_survival_event_censorship  <- ifelse(finalsurvivalanalysis$lpfs.time <= censorship & finalsurvivalanalysis$lpfs.status == 1 , 1 ,0  )
BCON_survival_time_censorship   <- ifelse( BCON_survival_event_censorship == 0 & finalsurvivalanalysis$lpfs.time >= censorship , censorship ,finalsurvivalanalysis$lpfs.time )   
finalsurvivalanalysis       <- cbind( finalsurvivalanalysis, BCON_survival_event_censorship, BCON_survival_time_censorship )
colnames(finalsurvivalanalysis)[ (ncol(finalsurvivalanalysis)-1):ncol(finalsurvivalanalysis) ] <- c("c_levent","c_levent_time")
#merge score data with the finalsurvivalanalysis
survival_scores<-as.data.frame(cbind(finalsurvivalanalysis[names(categories),],categories))
survival_scores$categories<-factor(survival_scores$categories, levels=c("Normoxia","Hypoxia"))
survival_scores<-as.data.frame(cbind(survival_scores[names(scores),],scores))

#high hypoxia
survival_scoreshigh<-survival_scores[survival_scores$categories=="Hypoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh)
par(pty="s")
gplotRT<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                    font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRT$plot <- gplotRT$plot + labs(
  title    = "Hypoxia"        
)
gplotRT <- ggpar(
  gplotRT,
  font.title    = c(25, "bold"))
gplotRT$plot<-gplotRT$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRT$table <- ggpar(gplotRT$table,
                       font.title = list(size = 16))

res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh)
par(pty="s")
gplotRTOS<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                      font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOS$plot <- gplotRTOS$plot + labs(
  title    = ""        
)
gplotRTOS <- ggpar(
  gplotRTOS,
  font.title    = c(25, "bold"))
gplotRTOS$plot<-gplotRTOS$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOS$table <- ggpar(gplotRTOS$table,
                         font.title = list(size = 16))
####multivariate analysis
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ CIS, data = survival_scoreshigh)
summary(res.cox)### 
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ Necrosis, data = survival_scoreshigh)
summary(res.cox)### 0.7
###
survival_scoreshigh$Stage<-as.factor(survival_scoreshigh$Stage)
stagecategory<-ifelse(survival_scoreshigh$Stage==1|survival_scoreshigh$Stage==2,"low","high")
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ stagecategory, data = survival_scoreshigh)
summary(res.cox)##0.9
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ age, data = survival_scoreshigh)
summary(res.cox)#0.2
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ Sex, data = survival_scoreshigh)
summary(res.cox)##
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ age+treatment, data = survival_scoreshigh)
summary(res.cox)#0

res.cox <- coxph(Surv(c_event_time, c_event) ~ CIS, data = survival_scoreshigh)
summary(res.cox)### 0.1
res.cox <- coxph(Surv(c_event_time, c_event) ~ Necrosis, data = survival_scoreshigh)
summary(res.cox)### 0.8
###
res.cox <- coxph(Surv(c_event_time, c_event) ~ stagecategory, data = survival_scoreshigh)
summary(res.cox)##0.8
res.cox <- coxph(Surv(c_event_time, c_event) ~ age, data = survival_scoreshigh)
summary(res.cox)##0.3
res.cox <- coxph(Surv(c_event_time, c_event) ~ Sex, data = survival_scoreshigh)
summary(res.cox)#


res.cox <- coxph(Surv(c_event_time, c_event) ~ age+treatment, data = survival_scoreshigh)
summary(res.cox)##0.


###hypoxia low
survival_scoreslow<-survival_scores[survival_scores$categories=="Normoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreslow)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreslow)
par(pty="s")
gplotRTlow<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                       risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                       risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                       xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                       font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTlow$plot <- gplotRTlow$plot + labs(
  title    = "Normoxia"        
)
gplotRTlow <- ggpar(
  gplotRTlow,
  font.title    = c(25, "bold"))
gplotRTlow$plot<-gplotRTlow$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTlow$table <- ggpar(gplotRTlow$table,
                          font.title = list(size = 16))

res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreslow)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreslow)
par(pty="s")
gplotRTOSlow<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                         risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOSlow$plot <- gplotRTOSlow$plot + labs(
  title    = ""        
)
gplotRTOSlow <- ggpar(
  gplotRTOSlow,
  font.title    = c(25, "bold"))
gplotRTOSlow$plot<-gplotRTOSlow$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOSlow$table <- ggpar(gplotRTOSlow$table,
                            font.title = list(size = 16))
####
splots<-list(gplotRT,gplotRTOS,gplotRTlow,gplotRTOSlow)
surv<-arrange_ggsurvplots(splots, print = TRUE,ncol=2,nrow=2)
ggsave("/Users/mkhan/Documents/miRNA hypoxia BJC figures/Figure 3.pdf", surv,width=12,height=12)
####
###look at the hypoxia markers in relation to it...
###look at them!!!
protein_expression<-read.table("/Users/mkhan/Documents/link_to_3proteins/BCON Trans microarray sample log_09052016_AW.csv",sep=",",header=T,stringsAsFactors = FALSE)
####
rownames(protein_expression)<-protein_expression$X
###add in the hypoxia scores
protein_expression_scores<-merge(protein_expression,scores,by="row.names")
"remove n/a"
protein_expression_scores$CA9<-gsub("n/a","",protein_expression_scores$CA9)
protein_expression_scores$CA9<-as.numeric(protein_expression_scores$CA9)
CA9<-protein_expression_scores[!is.na(protein_expression_scores$CA9),]
###stratify CA9 by median
CA9_category<-ifelse(as.numeric(CA9$CA9)>quantile(CA9$CA9,0.50),"high","low")
CA9$CA9_category<-CA9_category#
highCA9<-CA9[CA9$CA9_category=="high","y"]
lowCA9<-CA9[CA9$CA9_category=="low","y"]
### HIF-1 alpha expression
protein_expression_scores$HIF1A<-gsub("n/a","",protein_expression_scores$HIF1A)
protein_expression_scores$HIF1A<-as.numeric(protein_expression_scores$HIF1A)
HIF1A<-protein_expression_scores[!is.na(protein_expression_scores$HIF1A),]
###stratify HIF1 by lower quartile
HIF1A_category<-ifelse(as.numeric(HIF1A$HIF1A)>quantile(HIF1A$HIF1A,0.50),"high","low")
HIF1A$HIF1A_category<-HIF1A_category##p-value =0.03269
highHIF1A<-HIF1A[HIF1A$HIF1A_category=="high","y"]
lowHIF1A<-HIF1A[HIF1A$HIF1A_category=="low","y"]
####GLUT-1 expression
protein_expression_scores$Glut1<-gsub("n/a","",protein_expression_scores$Glut1)
protein_expression_scores$Glut1<-as.numeric(protein_expression_scores$Glut1)
Glut1<-protein_expression_scores[!is.na(protein_expression_scores$Glut1),]
###stratify HIF1 by lower quartile
Glut1_category<-ifelse(as.numeric(Glut1$Glut1)>quantile(Glut1$Glut1,0.50),"high","low")
Glut1$Glut1_category<-Glut1_category##No association
highGlut1<-Glut1[Glut1$Glut1_category=="high","y"]
lowGlut1<-Glut1[Glut1$Glut1_category=="low","y"]
###make a diagram showing the HIF-1 and Glut-1 expression in Prism
###look at these markers in patients with high hypoxia versus low hypoxia
hypoxia_category<-ifelse(protein_expression_scores$y>quantile(protein_expression_scores$y,0.75),"hypoxia","normoxia")
protein_expression_scores$hypoxia_category<-hypoxia_category
###HIF-1
HIF1_expression_hyp<-protein_expression_scores[protein_expression_scores$hypoxia_category=="hypoxia","HIF1A"]
HIF1_expression_hyp<-na.omit(HIF1_expression_hyp)
HIF1_expression_nor<-protein_expression_scores[protein_expression_scores$hypoxia_category=="normoxia","HIF1A"]
HIF1_expression_nor<-na.omit(HIF1_expression_nor)
###CA9
CA9_expression_hyp<-protein_expression_scores[protein_expression_scores$hypoxia_category=="hypoxia","CA9"]
CA9_expression_hyp<-na.omit(CA9_expression_hyp)
CA9_expression_nor<-protein_expression_scores[protein_expression_scores$hypoxia_category=="normoxia","CA9"]
CA9_expression_nor<-na.omit(CA9_expression_nor)
###Glut-1
protein_expression_scores$Glut1<-gsub("n/a","",protein_expression_scores$Glut1)
protein_expression_scores$Glut1<-as.numeric(protein_expression_scores$Glut1)
glut1_expression_hyp<-protein_expression_scores[protein_expression_scores$hypoxia_category=="hypoxia","Glut1"]
glut1_expression_hyp<-na.omit(glut1_expression_hyp)
glut1_expression_nor<-protein_expression_scores[protein_expression_scores$hypoxia_category=="normoxia","Glut1"]
glut1_expression_nor<-na.omit(glut1_expression_nor)

###study in RT only
survival_scores_RT<-survival_scores[survival_scores$treatment=="RT",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~categories , data = survival_scores_RT)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ categories, data = survival_scores_RT)
par(pty="s")
gplotRT<-ggsurvplot(fit,censor=TRUE,pval=T,legend.title="",legend.labs=c("Normoxia","Hypoxia")
                    ,linetype=c(3,1),break.time.by = 20, pval.size=7,pval.coord=c(0,0.03),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(), 
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",
                    xlab="Time in months", font.x = c("bold",22), font.y =c("bold",22), font.tickslab = "bold",
                    font.legend =  c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotRT$plot <- gplotRT$plot + labs(
  title    = "RT"           
)
gplotRT <- ggpar(
  gplotRT,
  font.title    = c(25, "bold"))
gplotRT$plot<-gplotRT$plot + theme(plot.title = element_text(hjust = 0.5))
gplotRT$table <- ggpar(gplotRT$table,
                       font.title = list(size = 16))   

res.cox <- coxph(Surv(c_event_time, c_event) ~ categories, data = survival_scores_RT)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ categories, data = survival_scores_RT)
par(pty="s")
gplotRTOS<-ggsurvplot(fit,censor=TRUE,pval=T,legend.title=""
                      ,legend.labs=c("Normoxia","Hypoxia")
                      ,linetype=c(3,1),break.time.by = 20, pval.size=7,pval.coord=c(0,0.03),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(), 
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",
                      xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                      font.legend =  c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotRTOS$plot <- gplotRTOS$plot + labs(
  title    = ""           
)
gplotRTOS <- ggpar(
  gplotRTOS,
  font.title    = c(25, "bold"))
gplotRTOS$plot<-gplotRTOS$plot + theme(plot.title = element_text(hjust = 0.5))
gplotRTOS$table <- ggpar(gplotRTOS$table,
                         font.title = list(size = 16))     

res.cox <- coxph(Surv(c_event_time, c_event) ~ CIS, data = survival_scores_RT)
summary(res.cox)#0.08

res.cox <- coxph(Surv(c_event_time, c_event) ~ age, data = survival_scores_RT)
summary(res.cox)#0.00113 **

res.cox <- coxph(Surv(c_event_time, c_event) ~ age+categories, data = survival_scores_RT)
summary(res.cox)#0.00113


survival_scores_RT$Stage<-as.factor(survival_scores_RT$Stage)
stagecategory<-ifelse(survival_scores_RT$Stage==1|survival_scores_RT$Stage==2,"low","high")
survival_scores_RT$stagecategory<-stagecategory
res.cox <- coxph(Surv(c_event_time, c_event) ~ stagecategory, data = survival_scores_RT)
summary(res.cox)#0.5
####necrosis yes
RT_nec<-survival_scores_RT[survival_scores_RT$Necrosis=="y",]
res.cox <- coxph(Surv(c_event_time, c_event) ~ categories, data = RT_nec)
summary(res.cox)#0.08

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ categories, data = RT_nec)
summary(res.cox)#0.08




###RT+CON arm
survival_scores_RTCON<-survival_scores[survival_scores$treatment=="RT+CON",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~categories , data = survival_scores_RTCON)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ categories, data = survival_scores_RTCON)
par(pty="s")
gplotCON<-ggsurvplot(fit,censor=TRUE,pval=T,legend.title=""
                     ,legend.labs=c("Normoxia","Hypoxia")
                     ,linetype=c(3,1), pval.size=7,pval.coord=c(0,0.03),break.time.by = 20,
                     risk.table.x.text = FALSE,tables.theme = clean_theme(), 
                     risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",
                     xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                     font.legend =  c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotCON$plot <- gplotCON$plot + labs(
  title    = "RT+CON"           
)
gplotCON <- ggpar(
  gplotCON,
  font.title    = c(25, "bold"))
gplotCON$plot<-gplotCON$plot + theme(plot.title = element_text(hjust = 0.5))
gplotCON$table <- ggpar(gplotCON$table,
                        font.title = list(size = 16))  

res.cox <- coxph(Surv(c_event_time, c_event) ~ categories, data = survival_scores_RTCON)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ categories, data = survival_scores_RTCON)
par(pty="s")
gplotCONS<-ggsurvplot(fit,censor=TRUE,pval=T,legend.title=""
                      ,legend.labs=c("Normoxia","Hypoxia")
                      ,linetype=c(3,1), pval.size=7,pval.coord=c(0,0.03),break.time.by = 20,
                      risk.table.x.text = FALSE,tables.theme = clean_theme(), 
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",
                      xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                      font.legend =  c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotCONS$plot <- gplotCONS$plot + labs(
  title    = ""           
)
gplotCONS <- ggpar(
  gplotCONS,
  font.title    = c(25, "bold"))
gplotCONS$plot<-gplotCONS$plot + theme(plot.title = element_text(hjust = 0.5))
gplotCONS$table <- ggpar(gplotCONS$table,
                         font.title = list(size = 16))   

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ CIS, data = survival_scores_RTCON)
summary(res.cox)#0.08

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ age, data = survival_scores_RTCON)
summary(res.cox)#0.00113 **

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ Sex, data = survival_scores_RTCON)
summary(res.cox)#0.00113 *

survival_scores_RTCON$Stage<-as.factor(survival_scores_RTCON$Stage)
stagecategory<-ifelse(survival_scores_RTCON$Stage==1|survival_scores_RTCON$Stage==2,"low","high")
survival_scores_RTCON$stagecategory<-stagecategory
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ stagecategory, data = survival_scores_RTCON)
summary(res.cox)#0.5
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ Necrosis, data = survival_scores_RTCON)
summary(res.cox)#0.5

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ age+categories, data = survival_scores_RTCON)
summary(res.cox)#0.5

res.cox <- coxph(Surv(c_event_time, c_event) ~ CIS, data = survival_scores_RTCON)
summary(res.cox)#0.08

res.cox <- coxph(Surv(c_event_time, c_event) ~ Necrosis, data = survival_scores_RTCON)
summary(res.cox)#


res.cox <- coxph(Surv(c_event_time, c_event) ~ age, data = survival_scores_RTCON)
summary(res.cox)#0.00113 **

res.cox <- coxph(Surv(c_event_time, c_event) ~ Sex, data = survival_scores_RTCON)
summary(res.cox)#0.00113


res.cox <- coxph(Surv(c_event_time, c_event) ~ stagecategory, data = survival_scores_RTCON)
summary(res.cox)#0.5


res.cox <- coxph(Surv(c_event_time, c_event) ~ age+categories, data = survival_scores_RTCON)
summary(res.cox)

splots<-list(gplotRT,gplotRTOS,gplotCON,gplotCONS)
s<-arrange_ggsurvplots(splots,ncol=2,nrow=2)
ggsave("/Users/mkhan/Documents/miRNA hypoxia BJC figures/S2.pdf",plot=s,height=12,width=12)

###DOES MIr-210 WORK
miR27<-c("hsa-miR-191-5p")
signature27<-BCONdata[miR27,-c(1:7)]
categories27<-ifelse(unlist(signature27)>quantile(unlist(signature27),0.75),"hypoxia","normoxia")
names(categories27)<-sub('X', '', names(categories27))
##get the survival data for BCON
survival_scores27<-as.data.frame(cbind(finalsurvivalanalysis[names(categories27),],categories27))
#high hypoxia
survival_scoreshigh27<-survival_scores27[survival_scores27$categories=="hypoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh27)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh27)
par(pty="s")
gplotRT<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("red","blue"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                    font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRT$plot <- gplotRT$plot + labs(
  title    = "BCON (hypoxia)"        
)
gplotRT <- ggpar(
  gplotRT,
  font.title    = c(25, "bold"))
gplotRT$plot<-gplotRT$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRT$table <- ggpar(gplotRT$table,
                       font.title = list(size = 16))

res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh27)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh27)
par(pty="s")
gplotRTOS<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("red","blue"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                      font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOS$plot <- gplotRTOS$plot + labs(
  title    = ""        
)
gplotRTOS <- ggpar(
  gplotRTOS,
  font.title    = c(25, "bold"))
gplotRTOS$plot<-gplotRTOS$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOS$table <- ggpar(gplotRTOS$table,
                         font.title = list(size = 16))
#####low hypoxia
survival_scoreslow27<-survival_scores27[survival_scores27$categories=="normoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreslow27)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreslow27)
par(pty="s")
gplotRTlow<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                       risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                       risk.table = TRUE,palette=c("red","blue"),ylab="Local relapse free survival",break.time.by = 20,
                       xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                       font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTlow$plot <- gplotRTlow$plot + labs(
  title    = "BCON (normoxia)"        
)
gplotRTlow <- ggpar(
  gplotRTlow,
  font.title    = c(25, "bold"))
gplotRTlow$plot<-gplotRTlow$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTlow$table <- ggpar(gplotRTlow$table,
                          font.title = list(size = 16))

res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreslow27)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreslow27)
par(pty="s")
gplotRTOSlow<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.69",legend.title="",linetype=c(1,3),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                         font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOSlow$plot <- gplotRTOSlow$plot + labs(
  title    = ""        
)
gplotRTOSlow <- ggpar(
  gplotRTOSlow,
  font.title    = c(25, "bold"))
gplotRTOSlow$plot<-gplotRTOSlow$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOSlow$table <- ggpar(gplotRTOSlow$table,
                            font.title = list(size = 16))
splots<-list(gplotRT,gplotRTOS,gplotRTlow,gplotRTOSlow)
surv<-arrange_ggsurvplots(splots, print = TRUE,ncol=2,nrow=2)
ggsave("training-validationbladderBCONmiR-210.pdf", surv,width=12,height=12)
###look at the hypoxia markers in relation to it...
###look at them!!!
protein_expression<-read.table("/Users/mkhan/Documents/link_to_3proteins/BCON Trans microarray sample log_09052016_AW.csv",sep=",",header=T)
####
rownames(protein_expression)<-protein_expression$X
###add in the hypoxia scores


# CLINICOPATHALOGIC VARIABLE DIFFERENCE IN HYPOXIA AND NORMOXIA CA --------
survival_scoreshigh$Stage<-as.factor(survival_scoreshigh$Stage)
stagecategory<-ifelse(survival_scoreshigh$Stage==1|survival_scoreshigh$Stage==2,"low","high")
survival_scoreshigh$stagecategory<-stagecategory

survival_scoreslow$Stage<-as.factor(survival_scoreslow$Stage)
stagecategory<-ifelse(survival_scoreslow$Stage==1|survival_scoreslow$Stage==2,"low","high")
survival_scoreslow$stagecategory<-stagecategory

chisq.test(cbind(table(survival_scoreshigh$stagecategory),table(survival_scoreslow$stagecategory)))


TURBTcategory<-ifelse(survival_scoreshigh$TURBT=="Complete","low","high")
survival_scoreshigh$TURBTcategory<-TURBTcategory

TURBTcategory<-ifelse(survival_scoreslow$TURBT=="Complete","low","high")
survival_scoreslow$TURBTcategory<-TURBTcategory

chisq.test(cbind(table(survival_scoreshigh$TURBTcategory),table(survival_scoreslow$TURBTcategory)))

chisq.test(cbind(table(survival_scoreshigh$Sex),table(survival_scoreslow$Sex)))

chisq.test(cbind(table(survival_scoreshigh$Necrosis),table(survival_scoreslow$Necrosis)))

chisq.test(cbind(table(survival_scoreshigh$CIS),table(survival_scoreslow$CIS)))


chisq.test(cbind(table(survival_scores_RT$TURBTcategory),table(survival_scores_RTCON$TURBTcategory)))



chisq.test(cbind(table(survival_scores_RT$TURBTcategory),table(survival_scores_RTCON$TURBTcategory)))

chisq.test(cbind(table(survival_scores_RT$Sex),table(survival_scores_RTCON$Sex)))

chisq.test(cbind(table(survival_scores_RT$Necrosis),table(survival_scores_RTCON$Necrosis)))

chisq.test(cbind(table(survival_scores_RT$CIS),table(survival_scores_RTCON$CIS)))

survival_scores_RT$Stage[which(grepl("T4b",survival_scores_RT$Stage))]

t.test(survival_scores_RT$age,survival_scores_RTCON$age)

survival_scores_RT$Stage<-as.character(survival_scores_RT$Stage)
stages_RT<-sapply(survival_scores_RT$Stage,function(x){
  if(x=="4a"|x=="4b"){y="4"}else{y=x}})

stages_RT_CON<-sapply(survival_scores_RTCON$Stage,function(x){
  if(x=="4a"){y="4"}else{y=x}})
chisq.test(cbind(table(stages_RT),table(stages_RT_CON)))

# Yang mRNA signature -----------------------------------------------------
DataBLCAhyp<-get(load("/Users/mkhan/Documents/fifth year/BCON_aroma.RData"))
load("/Users/mkhan/Documents/fourth year /miRNAs/BLCA signature/BCON_survival.RData")
#censor BCON survival data
censorship<-60
tr_survival_event_censorship  <- ifelse(  BCON_survival$lpfs.time <= censorship & BCON_survival$lpfs.status == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & BCON_survival$lpfs.time >= censorship, censorship , BCON_survival$lpfs.time )   
BCON_survival                  <- cbind( BCON_survival, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( BCON_survival )[ (ncol(BCON_survival)-1):ncol(BCON_survival) ] <- c( "censored_time", "censored_status" )
censorship<-60
tr_survival_event_censorshipo  <- ifelse(  BCON_survival$os.time <= censorship & BCON_survival$os.status == 1, 1 ,0)
tr_survival_time_censorshipo   <- ifelse(  tr_survival_event_censorshipo == 0 & BCON_survival$os.time >= censorship, censorship , BCON_survival$os.time )   
BCON_survival                   <- cbind(BCON_survival, tr_survival_time_censorshipo, tr_survival_event_censorshipo )
colnames( BCON_survival)[ (ncol(BCON_survival)-1):ncol(BCON_survival) ] <- c( "censored_otime", "censored_ostatus" )

mRNAhypoxiasig<-c("CAV1",
                  "COL5A1",
                  "ITGA5",
                  "SLC16A1",
                  "P4HA2",
                  "TGFBI",
                  "DPYSL2",
                  "SRPX",
                  "TRAM2",
                  "SYDE1",
                  "LRP1",
                  "PDLIM2",
                  "SAV1",
                  "AHNAK2",
                  "CAD",
                  "CYP1B1",
                  "DAAM1",
                  "DSC2",
                  "SLC2A3",
                  "FUT11",
                  "GLG1",
                  "GULP1",
                  "LDLR",
                  "THBS4")
DataBLCAhyp<-DataBLCAhyp[mRNAhypoxiasig,2:153]###already log transformed
bladderid<-read.table("/Users/mkhan/Documents/ third year/bladder data/bladder infofiles/sampleiDbladder.csv",sep=",",quote="",header=T,stringsAsFactors = FALSE)
bladderpatientno<-sapply(colnames(DataBLCAhyp),function(x){
  NAME<-strsplit(x, fixed=T, ".")[[1]][[1]]
  sampleid<-bladderid[bladderid$sampleCEL==NAME,]
  q<-sampleid['sampleIDs']})
bladderpatientno<-unlist(bladderpatientno)
colnames(DataBLCAhyp)<-bladderpatientno
DataBLCAhyp<-DataBLCAhyp[,intersect(rownames(BCON_survival),colnames(DataBLCAhyp))]
#calculate the median and categories together
medianexpression<-apply(DataBLCAhyp,2,median)
names(medianexpression)<-colnames(DataBLCAhyp)
stratification<-median(medianexpression)
category<-ifelse(medianexpression>=median(medianexpression),"High score","Low score")
###names of the BCON samples to be added
combinationcox<-cbind(BCON_survival[intersect(rownames(BCON_survival),names(category)),],category[intersect(rownames(BCON_survival),names(category))])
colnames(combinationcox)[ncol(combinationcox)]<-"y"
# comparison with Yang ----------------------------------------------------
##combinationcox for Yang 
mRNA_lncRNA<-cbind(survival_scores[intersect(rownames(combinationcox),rownames(survival_scores)),],combinationcox[intersect(rownames(combinationcox),rownames(survival_scores)),"y"])
colnames(mRNA_lncRNA)[ncol(mRNA_lncRNA)]<-"Yang"
mRNA_lncRNA$Yang<-factor(mRNA_lncRNA$Yang, levels=c("Low score","High score"))
###test for interaction
###high categories
miRNA<-mRNA_lncRNA[mRNA_lncRNA$categories=="Hypoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = miRNA)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = miRNA)
par(pty="s")
gplotRTmiRNA<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                    font.legend= c( "bold", "black",12),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTmiRNA$plot <- gplotRTmiRNA$plot + labs(
  title    = "miRNA signature"        
)
gplotRTmiRNA <- ggpar(
  gplotRTmiRNA,
  font.title    = c(25, "bold"))
gplotRTmiRNA$plot<-gplotRTmiRNA$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTmiRNA$table <- ggpar(gplotRTmiRNA$table,
                       font.title = list(size = 16))
res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = miRNA)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = miRNA)
par(pty="s")
gplotRTOSmiRNA<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                      font.legend= c( "bold", "black",12),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOSmiRNA$plot <- gplotRTOSmiRNA$plot + labs(
  title    = ""        
)
gplotRTOSmiRNA <- ggpar(
  gplotRTOSmiRNA,
  font.title    = c(25, "bold"))
gplotRTOSmiRNA$plot<-gplotRTOSmiRNA$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOSmiRNA$table <- ggpar(gplotRTOSmiRNA$table,
                         font.title = list(size = 16))
#####look at the mRNA signature
mRNA<-mRNA_lncRNA[mRNA_lncRNA$Yang=="High score",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = mRNA)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = mRNA)
par(pty="s")
gplotRTmRNA<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                    font.legend= c( "bold", "black",12),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTmRNA$plot <- gplotRTmRNA$plot + labs(
  title    = "Yang mRNA signature"        
)
gplotRTmRNA <- ggpar(
  gplotRTmRNA,
  font.title    = c(25, "bold"))
gplotRTmRNA$plot<-gplotRTmRNA$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTmRNA$table <- ggpar(gplotRTmRNA$table,
                       font.title = list(size = 16))
res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = mRNA)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = mRNA)
par(pty="s")
gplotRTOSmRNA<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                      font.legend= c( "bold", "black",12),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOSmRNA$plot <- gplotRTOSmRNA$plot + labs(
  title    = ""        
)
gplotRTOSmRNA <- ggpar(
  gplotRTOSmRNA,
  font.title    = c(25, "bold"))
gplotRTOSmRNA$plot<-gplotRTOSmRNA$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOSmRNA$table <- ggpar(gplotRTOSmRNA$table,
                         font.title = list(size = 16))
####using both
miRNA_mRNA<-mRNA_lncRNA[mRNA_lncRNA$Yang=="High score"&mRNA_lncRNA$categories=="Hypoxia",]
res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = miRNA_mRNA)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = miRNA_mRNA)
par(pty="s")
gplotRT<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                    font.legend= c( "bold", "black",12),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRT$plot <- gplotRT$plot + labs(
  title    = "combined"        
)
gplotRT <- ggpar(
  gplotRT,
  font.title    = c(25, "bold"))
gplotRT$plot<-gplotRT$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRT$table <- ggpar(gplotRT$table,
                       font.title = list(size = 16))
res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = miRNA_mRNA)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = miRNA_mRNA)
par(pty="s")
gplotRTOS<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                      font.legend= c( "bold", "black",12),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOS$plot <- gplotRTOS$plot + labs(
  title    = ""        
)
gplotRTOS <- ggpar(
  gplotRTOS,
  font.title    = c(25, "bold"))
gplotRTOS$plot<-gplotRTOS$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOS$table <- ggpar(gplotRTOS$table,
                         font.title = list(size = 16))

splots<-list(gplotRTmiRNA,gplotRTOSmiRNA,gplotRTmRNA,gplotRTOSmRNA,gplotRT,gplotRTOS)
s<-arrange_ggsurvplots(splots,ncol=3,nrow=2)
ggsave("/Users/mkhan/Documents/miRNA hypoxia BJC figures/Figure 5.pdf",plot=s,height=14,width=14)
###Test for interactionL miRNA
model<- coxph(Surv(c_levent_time, c_levent) ~ categories, data = mRNA_lncRNA)
model1 <- coxph(Surv(c_levent_time, c_levent) ~ categories*treatment, data = mRNA_lncRNA)
anova(model,model1)##0.1268
model<- coxph(Surv(c_event_time, c_event) ~ categories, data = mRNA_lncRNA)
model1 <- coxph(Surv(c_event_time, c_event) ~ categories*treatment, data = mRNA_lncRNA)
anova(model,model1)###0.06677
####Test for interaction: Yang
model<- coxph(Surv(c_levent_time, c_levent) ~ Yang, data = mRNA_lncRNA)
model1 <- coxph(Surv(c_levent_time, c_levent) ~ Yang*treatment, data = mRNA_lncRNA)
anova(model,model1)##0.02453
model<- coxph(Surv(c_event_time, c_event) ~ Yang, data = mRNA_lncRNA)
model1 <- coxph(Surv(c_event_time, c_event) ~ Yang*treatment, data = mRNA_lncRNA)
anova(model,model1)###0.04379
###Test for interaction: Both
mRNA_lncRNA$prognosis<-ifelse(mRNA_lncRNA$Yang=="High score"&mRNA_lncRNA$categories=="Hypoxia","High","Low")
model<- coxph(Surv(c_event_time, c_event) ~ prognosis, data = mRNA_lncRNA)
model1 <- coxph(Surv(c_event_time, c_event) ~ prognosis*treatment, data = mRNA_lncRNA)
anova(model,model1)##0.03316
model<- coxph(Surv(c_levent_time, c_levent) ~ prognosis, data = mRNA_lncRNA)
model1 <- coxph(Surv(c_levent_time, c_levent) ~ prognosis*treatment, data = mRNA_lncRNA)
anova(model,model1)##0.0649

# Protein expression most hypoxic tumours ---------------------------------
protein_expression_scores<-merge(protein_expression,mRNA_lncRNA,by="row.names",stringsAsFactors=F)
protein_expression_scores$HIF1A<-gsub("n/a","",protein_expression_scores$HIF1A)
protein_expression_scores$HIF1A<-as.numeric(protein_expression_scores$HIF1A)
##HIF
HIF1_expression_hyp<-protein_expression_scores[protein_expression_scores$prognosis=="High","HIF1A"]
HIF1_expression_hyp<-na.omit((HIF1_expression_hyp))
HIF1_expression_nor<-protein_expression_scores[protein_expression_scores$prognosis=="Low","HIF1A"]
HIF1_expression_nor<-na.omit(HIF1_expression_nor)##0.007
###CA9
protein_expression_scores$CA9<-gsub("n/a","",protein_expression_scores$CA9)
protein_expression_scores$CA9<-as.numeric(protein_expression_scores$CA9)
CA9_expression_hyp<-protein_expression_scores[protein_expression_scores$prognosis=="High","CA9"]
CA9_expression_hyp<-na.omit(CA9_expression_hyp)
CA9_expression_nor<-protein_expression_scores[protein_expression_scores$prognosis=="Low","CA9"]
CA9_expression_nor<-na.omit(CA9_expression_nor)###0.0008374
###Glut-1
protein_expression_scores$Glut1<-gsub("n/a","",protein_expression_scores$Glut1)
protein_expression_scores$Glut1<-as.numeric(protein_expression_scores$Glut1)
glut1_expression_hyp<-protein_expression_scores[protein_expression_scores$prognosis=="High","Glut1"]
glut1_expression_hyp<-na.omit(glut1_expression_hyp)
glut1_expression_nor<-protein_expression_scores[protein_expression_scores$prognosis=="Low","Glut1"]
glut1_expression_nor<-na.omit(glut1_expression_nor)##0.1582
# Protein expression Yang hypoxic tumours ---------------------------------
HIF1_expression_hyp<-protein_expression_scores[protein_expression_scores$Yang=="High score","HIF1A"]
HIF1_expression_hyp<-na.omit((HIF1_expression_hyp))
HIF1_expression_nor<-protein_expression_scores[protein_expression_scores$Yang=="Low score","HIF1A"]
HIF1_expression_nor<-na.omit(HIF1_expression_nor)##0.67

protein_expression_scores$CA9<-gsub("n/a","",protein_expression_scores$CA9)
protein_expression_scores$CA9<-as.numeric(protein_expression_scores$CA9)
CA9_expression_hyp<-protein_expression_scores[protein_expression_scores$Yang=="High score","CA9"]
CA9_expression_hyp<-na.omit(CA9_expression_hyp)
CA9_expression_nor<-protein_expression_scores[protein_expression_scores$Yang=="Low score","CA9"]
CA9_expression_nor<-na.omit(CA9_expression_nor)###0.18

###miRNA
HIF1_expression_hyp<-protein_expression_scores[protein_expression_scores$categories=="Hypoxia","HIF1A"]
HIF1_expression_hyp<-na.omit((HIF1_expression_hyp))
HIF1_expression_nor<-protein_expression_scores[protein_expression_scores$categories=="Normoxia","HIF1A"]
HIF1_expression_nor<-na.omit(HIF1_expression_nor)##0.007

protein_expression_scores$CA9<-gsub("n/a","",protein_expression_scores$CA9)
protein_expression_scores$CA9<-as.numeric(protein_expression_scores$CA9)
CA9_expression_hyp<-protein_expression_scores[protein_expression_scores$categories=="Hypoxia","CA9"]
CA9_expression_hyp<-na.omit(CA9_expression_hyp)
CA9_expression_nor<-protein_expression_scores[protein_expression_scores$categories=="Normoxia","CA9"]
CA9_expression_nor<-na.omit(CA9_expression_nor)###0.0008










###looking at radiotherapy performance
survival_scores_RTCON_combined<-mRNA_lncRNA[mRNA_lncRNA$treatment=="RT+CON",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~categories , data = survival_scores_RTCON_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time, c_event) ~categories , data = survival_scores_RTCON_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_levent_time, c_levent) ~Yang , data = survival_scores_RTCON_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time, c_event) ~Yang , data = survival_scores_RTCON_combined)
summary(res.cox)

survival_scores_RTCON_combined$prognosis<-ifelse(survival_scores_RTCON_combined$Yang=="High score"&survival_scores_RTCON_combined$categories=="Hypoxia","High","Low")

res.cox <- coxph(Surv(c_levent_time, c_levent) ~prognosis , data = survival_scores_RTCON_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time, c_event) ~prognosis , data = survival_scores_RTCON_combined)
summary(res.cox)

###RT
survival_scores_RT_combined<-mRNA_lncRNA[mRNA_lncRNA$treatment=="RT",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~categories , data = survival_scores_RT_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time, c_event) ~categories , data = survival_scores_RT_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_levent_time, c_levent) ~Yang , data = survival_scores_RT_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time, c_event) ~Yang , data = survival_scores_RT_combined)
summary(res.cox)

survival_scores_RT_combined$prognosis<-ifelse(survival_scores_RT_combined$Yang=="High score"&survival_scores_RT_combined$categories=="Hypoxia","High","Low")

res.cox <- coxph(Surv(c_levent_time, c_levent) ~prognosis , data = survival_scores_RT_combined)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time, c_event) ~prognosis , data = survival_scores_RT_combined)
summary(res.cox)


# NMIBC signature ---------------------------------------------------------

NMBC<-c("hsa-miR-210-3p","hsa-miR-193b-3p","hsa-miR-125a-3p","hsa-miR-708-5p","hsa-miR-517a-3p","hsa-miR-145-5p")
signature_NMBC<-BCONdata[NMBC,-c(1:7)]
scores<-apply(signature_NMBC,2,mean)
names(scores)<-sub('X', '', names(scores))
categories<-ifelse(scores>quantile(scores,0.75),"Hypoxia","Normoxia")
survival_scores<-as.data.frame(cbind(finalsurvivalanalysis[names(categories),],categories))
survival_scores$categories<-factor(survival_scores$categories, levels=c("Normoxia","Hypoxia"))
survival_scores<-as.data.frame(cbind(survival_scores[names(scores),],scores))

#high hypoxia
survival_scoreshigh<-survival_scores[survival_scores$categories=="Hypoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh)
par(pty="s")
gplotRT<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                    font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRT$plot <- gplotRT$plot + labs(
  title    = "Hypoxia"        
)
gplotRT <- ggpar(
  gplotRT,
  font.title    = c(25, "bold"))
gplotRT$plot<-gplotRT$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRT$table <- ggpar(gplotRT$table,
                       font.title = list(size = 16))

res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh)
par(pty="s")
gplotRTOS<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                      font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOS$plot <- gplotRTOS$plot + labs(
  title    = ""        
)
gplotRTOS <- ggpar(
  gplotRTOS,
  font.title    = c(25, "bold"))
gplotRTOS$plot<-gplotRTOS$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOS$table <- ggpar(gplotRTOS$table,
                         font.title = list(size = 16))
#low hypoxia
survival_scoreshigh<-survival_scores[survival_scores$categories=="Normoxia",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh)
summary(res.cox)
fit <- survfit(Surv(c_levent_time, c_levent) ~ treatment, data = survival_scoreshigh)
par(pty="s")
gplotRT<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                    risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                    risk.table = TRUE,palette=c("black","gray"),ylab="Local relapse free survival",break.time.by = 20,
                    xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                    font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRT$plot <- gplotRT$plot + labs(
  title    = "Hypoxia"        
)
gplotRT <- ggpar(
  gplotRT,
  font.title    = c(25, "bold"))
gplotRT$plot<-gplotRT$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRT$table <- ggpar(gplotRT$table,
                       font.title = list(size = 16))

res.cox <- coxph(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh)
summary(res.cox)
fit <- survfit(Surv(c_event_time, c_event) ~ treatment, data = survival_scoreshigh)
par(pty="s")
gplotRTOS<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),
                      risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.03),
                      risk.table = TRUE,palette=c("black","gray"),ylab="Overall survival",break.time.by = 20,
                      xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                      font.legend= c( "bold", "black",14),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotRTOS$plot <- gplotRTOS$plot + labs(
  title    = ""        
)
gplotRTOS <- ggpar(
  gplotRTOS,
  font.title    = c(25, "bold"))
gplotRTOS$plot<-gplotRTOS$plot + theme(plot.title = element_text(hjust = 0.5))

gplotRTOS$table <- ggpar(gplotRTOS$table,
                         font.title = list(size = 16))
####RT arm

RT<-survival_scores[survival_scores$treatment=="RT",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ categories, data = RT)
summary(res.cox)
res.cox <- coxph(Surv(c_event_time, c_event) ~ categories, data = RT)
summary(res.cox)

###RT+CON
RTCON<-survival_scores[survival_scores$treatment=="RT+CON",]

res.cox <- coxph(Surv(c_levent_time, c_levent) ~ categories, data = RTCON)
summary(res.cox)
res.cox <- coxph(Surv(c_event_time, c_event) ~ categories, data = RTCON)
summary(res.cox)



