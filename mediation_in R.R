
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", "plyr", "reshape2", "RColorBrewer", "scales", "grid",
              'Rcmdr','Hmisc','lavaan')
ipak(packages)

library(ggplot2)
library(Rcmdr)
library(Hmisc)
library(reshape2)
library(lavaan)


med_master_v3 <- read.csv("MED_MASTER_v3.csv", stringsAsFactors = FALSE )
View(med_master_v3)
colnames(med_master_v3) <- tolower(colnames(med_master_v3))

med_bdi_join <- read.csv("MED_BDI_JOIN.csv", stringsAsFactors = FALSE)
colnames(med_bdi_join) <- tolower(colnames(med_bdi_join))
med_bdi <- read.csv("MED_BDI.csv", stringsAsFactors = FALSE)
colnames(med_bdi) <- tolower(colnames(med_bdi))

library(tidyr)

medbdi_wide <- spread(med_bdi,bdi_label,bdi_total)
View(medbdi_wide)
names(medbdi_wide)
medbdi_wide2 <- medbdi_wide[,c('patno','BDI_s1',	'BDI_s2',	'BDI_s3',	'BDI_s4',	'BDI_s5',	'BDI_s6',	'BDI_s7',	'BDI_s8',	'BDI_s9',	'BDI_s10',	'BDI_s11',	'BDI_s12',	'BDI_s13',	'BDI_s14',	'BDI_s15',	'BDI_s16',	'BDI_s17',	'BDI_s18',	'BDI_s19',	'BDI_s20',	'BDI_s21',	'BDI_s22',	'BDI_s23',	'BDI_s24','BDI_s25',	'BDI_s26',	'BDI_s27',	'BDI_s28')]
View(medbdi_wide2)
save(medbdi_wide2, file="medbdi_wide2.Rdata")
medbdi_wide3 <- medbdi_wide2[c(-1)]
View(medbdi_wide3)

med_master_v4 <- cbind(med_master_v3,medbdi_wide3)
View(med_master_v4)
save(med_master_v4, file="med_master_v4.Rdata")

med_preditems <- read.csv("MED_PREDitems.csv", stringsAsFactors = FALSE)
View(med_preditems)
med_master_v5 <- cbind(med_master_v4,med_preditems)
View(med_master_v5)
View(med_master_v5[,-(1:100)])
load(file = "med_master_v5.Rdata")

length(which(med_master_v5$BDI_s6 %in% 0:100))

## FITTING 2-2-1 MSEM MODEL 

# *VAR NAMES UNIQUE TO  MED

# CCavg_1
# CCavg_2 (average across all 5 sessions on item 2)
# CCavg_tot (Average across all 5 session on TOTAL score)
# CCtot_s1  (TOTAL score at session 1 specifically)
# WAIavg_2 	(average across all 5 sessions on item 2)
# WAIavg_tot (Average across all 5 session on TOTAL score)
# WAItot_s1  (TOTAL score at session 1 specifically)
# BDI_s1  (BDI TOTAL SCORE AT SESSION 1)	
  
########################################

#adding in reverse scored vars so that HIGH IV = GREATER VULNERABILITY FOR ALL INDICATORS;
names(med_master_v5)
summary(med_master_v5[,c(3:12)])

med_master_v5$efa_protrx_r <- 7.480-med_master_v5$efa_protrx
med_master_v5$efa_sociability_r <- 6.714-med_master_v5$efa_sociability
med_master_v5$imp_sit_mean_r <- 7.718-med_master_v5$imp_sit_mean
med_master_v5$t_efa_protrx_r <- 7.400-med_master_v5$t_efa_protrx
med_master_v5$t_efa_extravert_r <- 8.000-med_master_v5$t_efa_extravert  
med_master_v5$pidintake_total <- med_master_v5$pidintake_total-24
med_master_v5$iip_total <- med_master_v5$iip_total-31
head(med_master_v5$pidintake_total)
head(med_master_v5$iip_total)

View(med_master_v5)
save(med_master_v5,file="med_master_v5.Rdata")
load(file = "med_master_v5.Rdata")

#SEM 1 : PREACHER RECOMMENDATION 
#Hypothesis being tested: clients average alliance and in-session cognitive change scores 
#account for identified relations between IV at intake and subsequent slope of symptom change

#########PLOTTING INDIVIDUAL SLOPES - where x = time var 
lattice::xyplot(med, y ~ x|id, type=c("p", "l"), layout = c(5,5))
getwd()
med_nonlin <- read.csv("med_Bdi_nonlin.csv", stringsAsFactors = FALSE )
colnames(med_nonlin) <- tolower(colnames(med_nonlin))
names(med_nonlin)
med_nonlin$patno <- as.factor(med_nonlin$patno)
str(med_nonlin)
med_nonlin$bdi_total <- as.numeric(med_nonlin$bdi_total)
table(med_nonlin$time)

lattice::xyplot(bdi_total ~ time|patno, type=c("p", "r"), layout = c(5,5), data = med_nonlin)

xyplot(bdi_total ~ time|patno, type=c("p", "r"), layout = c(5,5), data = med_nonlin)
##SWITCH TO GGPLOT IF PUTTIN GIN GRAPH (aes(x = time y = output ,   )) - WRAPPER (
## GGPLOT TO TRELLIS OR PANEL REGERSSION TO GGPLOT 

# TRANSFORMATION OPTION TO TRY IF ANY :
#polynomial latent growth cruve  (degree 5 - can it smooth over and fit everyone, time up to the 5th power

# OTHER TRANSFORMATIONS
  # curvilinear line with a linear fucntion ) 
  #quadratic function  (exponential decay)
  # linear linear spline (slope and an intercept , knot, and a slope afte rthe knot - function can take on multiple forms) 



########################################################################################################
## ALLIANCE MEDIATION 
########################################################################################################

#PID ONly 
#FINAL PID ALLIANCE 

model221_LAT <- '
#a and b path 
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17


#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
## calling modifaction indices?  problem: sample specific  (driven by data)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)


#PID ONLY FINAL PID SUBSEQUENT 


model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent 
bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17


#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)



#################################################################################
## CC MEDIATION
#################################################################################

#PID ONLY


model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17


#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 
'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)
semPaths(model221_LAT,title=FALSE, curvePivot = TRUE)
library(qgraph)
library(semPlot)

.02*(-.01)
.02*
  #CC subsequent

#PID ONLY


model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5



#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent 
bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17


#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)

semPaths(model221_LAT,title=FALSE, curvePivot = TRUE)



# TRANSFORMATIONs TRIED BELOW (POLYNOMIAL AND QUADRATIC ) :
#polynomial latent growth cruve  (degree 5 - can it smooth over and fit everyone, time up to the 5th power
med_master_v5$BDI_intake <- med_master_v5$bdi_intake
#quadratic 
names(med_master_v5)

bdi_quadnames <- grep('BDI',colnames(med_master_v5[,-189]),value = T)
bdi_quadnames2<-paste((1:17)^2,'*',bdi_quadnames,sep='')
bdi_quadnames2<-paste((0:17)^2,'*',bdi_quadnames,sep='')
bdi_quadnames3<-paste(bdi_quadnames2,collapse=' + ')
bdi_quadnames3

model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M
bdi_q ~ bq*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X
bdi_q ~ cq*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs
ab_qs := a*bq

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)
total_qs := cq + (a*bq)

#M 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17

bdi_q =~ 0*bdi_intake+ 1*BDI_s1 + 4*BDI_s2 + 9*BDI_s3 + 16*BDI_s4 + 25*BDI_s5 + 36*BDI_s6 +
49*BDI_s7 + 64*BDI_s8 + 81*BDI_s9 + 100*BDI_s10 + 121*BDI_s11 + 144*BDI_s12 + 169*BDI_s13 + 
196*BDI_s14 + 225*BDI_s15 + 256*BDI_s16 + 289*BDI_s17 


#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_q ~~ bdi_q 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 
'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)


###TRYING A FIFTH DEGREE POLYNOMIAL


#polynomial latent growth cruve  (degree 5 - can it smooth over and fit everyone, time up to the 5th power
med_master_v5$BDI_intake <- med_master_v5$bdi_intake


bdi_quadnames <- grep('BDI',colnames(med_master_v5[,-189]),value = T)
bdi_quadnames2<-paste((1:17)^3,'*',bdi_quadnames,sep='')
bdi_quadnames2<-paste((1:17)^3,'*',bdi_quadnames,sep='')
bdi_quadnames3<-paste(bdi_quadnames2,collapse=' + ')
bdi_quadnames3

model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M
bdi_q ~ bq*M
bdi_p ~ bp*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X
bdi_q ~ cq*X
bdi_p ~ cp*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs
ab_qs := a*bq
ab_ps := a*bp


#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)
total_qs := cq + (a*bq)
total_ps := cp + (a*bp)

#M 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17

bdi_q =~ 0*bdi_intake+ 1*BDI_s1 + 4*BDI_s2 + 9*BDI_s3 + 16*BDI_s4 + 25*BDI_s5 + 36*BDI_s6 +
49*BDI_s7 + 64*BDI_s8 + 81*BDI_s9 + 100*BDI_s10 + 121*BDI_s11 + 144*BDI_s12 + 169*BDI_s13 + 
196*BDI_s14 + 225*BDI_s15 + 256*BDI_s16 + 289*BDI_s17 

bdi_p =~ 0*bdi_intake+ 1*BDI_s1 + 8*BDI_s2 + 27*BDI_s3 + 64*BDI_s4 + 125*BDI_s5 + 
216*BDI_s6 + 343*BDI_s7 + 512*BDI_s8 + 729*BDI_s9 + 1000*BDI_s10 + 1331*BDI_s11 + 
1728*BDI_s12 + 2197*BDI_s13 + 2744*BDI_s14 + 3375*BDI_s15 + 4096*BDI_s16 + 4913*BDI_s17

#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_q ~~ bdi_q 
bdi_p ~~ bdi_p 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 
'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)





#################################################
#ADDITIONAL MODELS EXAMINED
#################################################


#ALL SCALES DEFINE X 

model221_LAT <- '
  #a and b path - mediator effects
    M ~ a*X
    bdi_i ~ bi*M
    bdi_s ~ bs*M
  
  #c path - direct effect(s)
    bdi_i ~ ci*X
    bdi_s ~ cs*X
  
  #indirect effect (a*b)
    ab_i := a*bi
    ab_s := a*bs

  #total effect 
    total_i := ci + (a*bi)
    total_s := cs + (a*bs)

  #M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
    M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
    1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 
    
  
  #X latent  
  #fixing 1 and not the rest - - scales latent in specified manner - allows all scales to contribute differentially to x
    X =~ 1*efa_protrx_r +	efa_sociability_r + efa_neurotic +	t_efa_protrx_r +	
      t_efa_extravert_r +	t_efa_neurotic +	imp_sit_mean_r +	pidintake_total+	iip_total
  
  
  #Ylatent --- intercept average of all session scores from 0-18 simple average 
    bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
    1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
    1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 

  
    bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
              8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
              15*BDI_s15+	16*BDI_s16+	17*BDI_s17+	18*BDI_s18

  
  #residual variances and covariances
   bdi_i ~~ bdi_i 
      bdi_s ~~ bdi_s 

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)


#rmsea = .05  good, .08 = fair , 1 = bad; cfi = .95

#####DOES AVERAGE ALLIANCE ACCOUNT FOR ANY SIGNIFICANT iv-OUTCOME RELATIONS IDENFIED ABOVE

#TS-therapist neuroticism only
#START HERE QUESTIONS

model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
    X =~ 1*T_emotional_7+	1*T_worry_criticism_15+	1*T_believes_inferior_16

#Ylatent --- intercept average of all session scores from 0-18 simple average 
    bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
    1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
    1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 
    

bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
        8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
        15*BDI_s15+	16*BDI_s16+	17*BDI_s17+	18*BDI_s18


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 
BDI_s18~~errvar*BDI_s18

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)
semPaths(model221_LAT,title=FALSE, curvePivot = TRUE)
coef(fit_model221_LAT)
?paste
resid(fit_model221_LAT,type = "standardized")

noquote(paste0('BDI_s',1:18,'~~errvar*BDI_s',1:18))



#errvar = residual variance

###THE DEPENDENT  VARIABLES RESIDUAL STRUCTURE IS THE SAME AS HLM HOMOSCEDASTIC UNCORRELATED (GIVEN THE SUBJECT SEPCIFIC SLOPE AND RANDOME EFFECT CONSTANT VARIANCE UNCORRELATED) 
##QUESTION CHECK INTERPRETATION OF SIGNIFICANT Total effect BUT NOT INDIRECT EFFECT
  #signifigance of total effect coupled with non-significant indirect effect suggest
    #the relationship is not a mediational one but rather all three variables are related to one another?
    #interpretation of effect: direct effect X latent variable predicing the slope and the intercept 
     

getwd()



#SIT 
model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
X =~ 1*initiate_convo_1+	1*apologize_2+	1*thank_someone_3+	1*assert_claim_4+	1*confrontsomeone_5+	1*present_positively_6+	1*be_selfcritical_7+	1*showempathy_8+	1*reprimand_9+	1*encourage_10+	1*give_instructs_11+	1*compliment_12+	1*askfor_EMOsupport_13+	1*convince_someone_14+	1*bargain_15+	1*expressaffection_16+	1*askfor_INSTsupport_17


#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 


'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)

#################################################################################
#ALLIANCE SUBSEQUENT
#################################################################################


model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
X =~ 1*T_emotional_7+	1*T_worry_criticism_15+	1*T_believes_inferior_16

    
 #Ylatent 
bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 


bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17+	12*BDI_s18


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)
semPaths(model221_LAT,title=FALSE, curvePivot = TRUE)
coef(fit_model221_LAT)
?paste

noquote(paste0('BDI_s',1:18,'~~errvar*BDI_s',1:18))



#SIT 
model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
X =~ 1*initiate_convo_1+	1*apologize_2+	1*thank_someone_3+	1*assert_claim_4+	1*confrontsomeone_5+	1*present_positively_6+	1*be_selfcritical_7+	1*showempathy_8+	1*reprimand_9+	1*encourage_10+	1*give_instructs_11+	1*compliment_12+	1*askfor_EMOsupport_13+	1*convince_someone_14+	1*bargain_15+	1*expressaffection_16+	1*askfor_INSTsupport_17


#Ylatent 
bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 


bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17+	12*BDI_s18


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 


'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)




#rmsea = .05  good, .08 = fair , 1 = bad; cfi = .95



#TS-therapist neuroticism only

model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*T_emotional_7+	1*T_worry_criticism_15+	1*T_believes_inferior_16

#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17+	18*BDI_s18


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 
BDI_s18~~errvar*BDI_s18

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)
semPaths(model221_LAT,title=FALSE, curvePivot = TRUE)
coef(fit_model221_LAT)
?paste

noquote(paste0('BDI_s',1:18,'~~errvar*BDI_s',1:18))

getwd()



#SIT 
model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*initiate_convo_1+	1*apologize_2+	1*thank_someone_3+	1*assert_claim_4+	1*confrontsomeone_5+	1*present_positively_6+	1*be_selfcritical_7+	1*showempathy_8+	1*reprimand_9+	1*encourage_10+	1*give_instructs_11+	1*compliment_12+	1*askfor_EMOsupport_13+	1*convince_someone_14+	1*bargain_15+	1*expressaffection_16+	1*askfor_INSTsupport_17


#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 



'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)

#################################################################################
#cc SUBSEQUENT
#################################################################################


model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*T_emotional_7+	1*T_worry_criticism_15+	1*T_believes_inferior_16


#Ylatent 
bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 


bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17+	12*BDI_s18


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)
semPaths(model221_LAT,title=FALSE, curvePivot = TRUE)
coef(fit_model221_LAT)
?paste

noquote(paste0('BDI_s',1:18,'~~errvar*BDI_s',1:18))




#SIT 
model221_LAT <- '
#a and b path - mediator effects
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*ccavg_1+	1*ccavg_2+	1*ccavg_3+	1*ccavg_4+	1*ccavg_5


#X latent  
X =~ 1*initiate_convo_1+	1*apologize_2+	1*thank_someone_3+	1*assert_claim_4+	1*confrontsomeone_5+	1*present_positively_6+	1*be_selfcritical_7+	1*showempathy_8+	1*reprimand_9+	1*encourage_10+	1*give_instructs_11+	1*compliment_12+	1*askfor_EMOsupport_13+	1*convince_someone_14+	1*bargain_15+	1*expressaffection_16+	1*askfor_INSTsupport_17


#Ylatent 
bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 


bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17+	12*BDI_s18


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)



###QUESTION == WHAT IF I HAD FOUND A SIGNIFICANT AB_S EFFECT- HOW TO INTERPERT?   
#the effect of interpersonal vulnerability on slope is conditional on cognitive change such that
#????HOW TO UNDERSTAND DIRECTION -- LOOK AT REGRESSIONS? - higher interpresonal 
#vulnerabiltiy --> greater cog change which predicts steeper slope of symptoms.DOES 
#persons latent growth curve is a function of x m and and a little bit of error no indirect effects -- seems to be operating independentyl

?qgraph.lavaan
library(qgraph)
library(semPlot)
?"semPaths
#http://sachaepskamp.com/documentation/semPlot/semPaths.html
semPaths(fit_model221cog_LAT,title=FALSE, curvePivot = TRUE)
##QUESTION == WHAT ELSE MIGHT i WANT? 


## if i dont change it -- itnercept is forecasting back  to zero 

###############MAYBE IF I WANTED TO DO REPEATED MEASURES subsequent CHECK 
# for any mediation model involving at least 1 level 2 var, the indirect effect can exist only at the between level
#at least one level-2 variable in the causal chain, the indirect effect of interest is likely a between indirect effect
names(med_master_v5)

model211cog_LATsub <- '
      #a and b path - mediator effects
        M_i ~ a_i*X       #effect of IV on CC avg
        M_s ~ a_s*X       #effect of IV on CC slope
        bdi_i ~ bi_mi*M_i  #effect of CC avg on bdi avg
        bdi_i ~ bi_ms*M_s  #effect of CC slope on bdi avg === TAKE THIS OUT? 
        bdi_s ~ bs_mi*M_i  #effect of CC avg on bdi slope ==  TAKE THIS OUT ?
        bdi_s ~ bs_ms*M_s   #effect of CC slope on bdi slope
        
        #c path - direct effect(s)
        bdi_i ~ ci*X
        bdi_s ~ cs*X

        #indirect effect (a*b)
          #via average cog change
          #Ai_mi_bi := a_i*bi_mi
          #Ai_ms_bi := a_i*bi_ms
          #Ai_mi_bs := a_i*bs_mi
          #Ai_ms_bs := a_i*bs_ms

          #via slope cog change
          As_mi_bi := a_s*bi_mi
          As_ms_bi := a_s*bi_ms
          As_mi_bs := a_s*bs_mi
          As_ms_bs := a_s*bs_ms       
          
        #total effect 
        #tot_ci_a_i_bi_mi := ci + (a_i*bi_mi)
        #tot_ci_a_i_bi_ms := ci + (a_i*bi_ms)
        #tot_ci_a_i_bs_mi := ci + (a_i*bs_mi)
        #tot_ci_a_i_bs_ms := ci + (a_i*bs_ms)

        tot_cs_a_i_bi_mi := cs + (a_i*bi_mi)
        tot_cs_a_i_bi_ms := cs + (a_i*bi_ms)
        tot_cs_a_i_bs_mi := cs + (a_i*bs_mi)
        tot_cs_a_i_bs_ms := cs + (a_i*bs_ms)

        #M latent- M_i is equivalent to M in the 2-1-1 models  (Average across all sessions)
M_i =~ 1*waitot_s1+	1*waitot_s2+	1*waitot_s3+	1*waitot_s4+	1*waitot_s5

M_s =~ 0*waitot_s1+	1*waitot_s2+	2*waitot_s3+	2*waitot_s4+	4*waitot_s5
        

 #X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25


        #Ylatent 
        bdi_i =~  1*BDI_s6+	1*BDI_s7+	1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	
        1*BDI_s13+	1*BDI_s14+1*BDI_s15+	1*BDI_s16+	1*BDI_s17+	1*BDI_s18 
        
        
        bdi_s =~  0*BDI_s6+	1*BDI_s7+	2*BDI_s8+ 3*BDI_s9+	4*BDI_s10+	5*BDI_s11+	6*BDI_s12+	
        7*BDI_s13+	8*BDI_s14+	9*BDI_s15+ 10*BDI_s16+	11*BDI_s17+	12*BDI_s18
        
        #residual variances and covariances
         bdi_i ~~ bdi_i 
        bdi_s ~~ bdi_s 
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 
BDI_s18~~errvar*BDI_s18

'
fit_model211cog_LATsub <- growth(model211cog_LATsub, data = med_master_v5)
summary(model211cog_LATsub, standardized = TRUE)

#WAI 2-1-1

model211cog_LATsub <- '
#a and b path - mediator effects
M_i ~ a_i*X       #effect of IV on CC avg
M_s ~ a_s*X       #effect of IV on CC slope
bdi_i ~ bi_mi*M_i  #effect of CC avg on bdi avg
#bdi_i ~ bi_ms*M_s  effect of CC slope on bdi avg === TAKE THIS OUT? 
#bdi_s ~ bs_mi*M_i  #effect of CC avg on bdi slope ==  TAKE THIS OUT ?
bdi_s ~ bs_ms*M_s   #effect of CC slope on bdi slope

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
#via average cog change
Ai_mi_bi := a_i*bi_mi
#Ai_ms_bi := a_i*bi_ms
#Ai_mi_bs := a_i*bs_mi
Ai_ms_bs := a_i*bs_ms

#via slope cog change
As_mi_bi := a_s*bi_mi
#As_ms_bi := a_s*bi_ms
#As_mi_bs := a_s*bs_mi
As_ms_bs := a_s*bs_ms       

#total effect 
tot_ci_a_i_bi_mi := ci + (a_i*bi_mi)
#tot_ci_a_i_bi_ms := ci + (a_i*bi_ms)
#tot_ci_a_i_bs_mi := ci + (a_i*bs_mi)
tot_ci_a_i_bs_ms := ci + (a_i*bs_ms)

tot_cs_a_i_bi_mi := cs + (a_i*bi_mi)
#tot_cs_a_i_bi_ms := cs + (a_i*bi_ms)
#tot_cs_a_i_bs_mi := cs + (a_i*bs_mi)
tot_cs_a_i_bs_ms := cs + (a_i*bs_ms)

#M latent- M_i is equivalent to M in the 2-1-1 models  (Average across all sessions)
M_i =~ 1*cctot_s1+	1*cctot_s2+	1*cctot_s3+	1*cctot_s4+	1*cctot_s5
M_s =~ 0*cctot_s1+	1*cctot_s2+	2*cctot_s3+	3*cctot_s4+	4*cctot_s5


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25


#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17


bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 

'
fit_model211cog_LATsub <- growth(model211cog_LATsub, data = med_master_v5)
summary(model211cog_LATsub, standardized = TRUE)

(-.05)*(-.002)
model211cog_LATsub <- '
#a and b path - mediator effects
M_i ~ a_i*X       #effect of IV on CC avg
M_s ~ a_s*X       #effect of IV on CC slope
bdi_i ~ bi_mi*M_i  #effect of CC avg on bdi avg
bdi_i ~ bi_ms*M_s  #effect of CC slope on bdi avg === TAKE THIS OUT? 
bdi_s ~ bs_mi*M_i  #effect of CC avg on bdi slope ==  TAKE THIS OUT ?
bdi_s ~ bs_ms*M_s   #effect of CC slope on bdi slope

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
#via average cog change
Ai_mi_bi := a_i*bi_mi
Ai_ms_bi := a_i*bi_ms
Ai_mi_bs := a_i*bs_mi
Ai_ms_bs := a_i*bs_ms

#via slope cog change
As_mi_bi := a_s*bi_mi
As_ms_bi := a_s*bi_ms
As_mi_bs := a_s*bs_mi
As_ms_bs := a_s*bs_ms       

#total effect 
tot_ci_a_i_bi_mi := ci + (a_i*bi_mi)
tot_ci_a_i_bi_ms := ci + (a_i*bi_ms)
tot_ci_a_i_bs_mi := ci + (a_i*bs_mi)
tot_ci_a_i_bs_ms := ci + (a_i*bs_ms)

tot_cs_a_i_bi_mi := cs + (a_i*bi_mi)
tot_cs_a_i_bi_ms := cs + (a_i*bi_ms)
tot_cs_a_i_bs_mi := cs + (a_i*bs_mi)
tot_cs_a_i_bs_ms := cs + (a_i*bs_ms)

#M latent- M_i is equivalent to M in the 2-1-1 models  (Average across all sessions)
M_i =~ 1*waitot_s1+	1*waitot_s2+	1*waitot_s3+	1*waitot_s4+	1*waitot_s5

M_s =~ 0*waitot_s1+	1*waitot_s2+	2*waitot_s3+	2*waitot_s4+	4*waitot_s5


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25


#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s1 + 1*BDI_s1+	1*BDI_s2+	1*BDI_s3+	1*BDI_s4+	1*BDI_s5+	1*BDI_s6+	1*BDI_s7+	
1*BDI_s8+	1*BDI_s9+	1*BDI_s10+	1*BDI_s11+	1*BDI_s12+	1*BDI_s13+	1*BDI_s14+	
1*BDI_s15+	1*BDI_s16+	1*BDI_s17


bdi_s =~ 0*bdi_intake+	1*BDI_s1+	2*BDI_s2+	3*BDI_s3+	4*BDI_s4+	5*BDI_s5+	6*BDI_s6+	7*BDI_s7+	
8*BDI_s8+	9*BDI_s9+	10*BDI_s10+	11*BDI_s11+	12*BDI_s12+	13*BDI_s13+	14*BDI_s14+	
15*BDI_s15+	16*BDI_s16+	17*BDI_s17


#residual variances and covariances
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s1~~errvar*BDI_s1   
BDI_s2~~errvar*BDI_s2   
BDI_s3~~errvar*BDI_s3   
BDI_s4~~errvar*BDI_s4   
BDI_s5~~errvar*BDI_s5  
BDI_s6~~errvar*BDI_s6   
BDI_s7~~errvar*BDI_s7   
BDI_s8~~errvar*BDI_s8   
BDI_s9~~errvar*BDI_s9   
BDI_s10~~errvar*BDI_s10
BDI_s11~~errvar*BDI_s11
BDI_s12~~errvar*BDI_s12
BDI_s13~~errvar*BDI_s13
BDI_s14~~errvar*BDI_s14
BDI_s15~~errvar*BDI_s15
BDI_s16~~errvar*BDI_s16 
BDI_s17~~errvar*BDI_s17 

'
fit_model211cog_LATsub <- growth(model211cog_LATsub, data = med_master_v5)
summary(model211cog_LATsub, standardized = TRUE)

#


T_EFA_protrx 	:	
T_trx_interest_1 	
T_no_reservations_2	
T_wont_dropout_3	
T_agreeable_5	
T_conscientious_6	
T_likeable_9
T_use_nonverbalwell_17
T_speaks_well_18
T_Good_convoskills_19
T_R_hardtoworkwith_20



#PID ONly 
#FINAL PID ALLIANCE 

model221_LAT <- '
#a and b path 
M ~ a*X
bdi_i ~ bi*M
bdi_s ~ bs*M

#c path - direct effect(s)
bdi_i ~ ci*X
bdi_s ~ cs*X

#indirect effect (a*b)
ab_i := a*bi
ab_s := a*bs

#total effect 
total_i := ci + (a*bi)
total_s := cs + (a*bs)

#M latent since all loadings fixed to 1 M reflects AVERAGE ALLIANCE 
M =~ 1*waiavg_1+	1*waiavg_2+	1*waiavg_3+	1*waiavg_4+	1*waiavg_5+	1*waiavg_6 +	
1*waiavg_7+	1*waiavg_8+	1*waiavg_9+	1*waiavg_10+	1*waiavg_11+	1*waiavg_12 


#X latent  
X =~ 1*PID_1+	1*PID_2+	1*PID_3+	1*PID_4+	1*PID_5+	1*PID_6+	1*PID_7+	1*PID_8+	1*PID_9+	1*PID_10+	1*PID_11+	1*PID_12+	1*PID_13+	1*PID_14+	1*PID_15+	1*PID_16+	1*PID_17+	1*PID_18+	1*PID_19+	1*PID_20+	1*PID_21+	1*PID_22+	1*PID_23+	1*PID_24+	1*PID_25



#Ylatent --- intercept average of all session scores from 0-18 simple average 
bdi_i =~  1*bdi_intake + 1*BDI_s4 + 1*BDI_s8+	1*BDI_s12+	1*BDI_s16

bdi_s =~  0*bdi_intake + 1*BDI_s4 + 2*BDI_s8+	3*BDI_s12+	4*BDI_s16


#residual variances and covariances
#QUESTION preacher suggestion specified that I should force all residual varianceos of Y to equality- that being done here? 
bdi_i ~~ bdi_i 
bdi_s ~~ bdi_s 
bdi_intake ~~errvar*bdi_intake
BDI_s4~~errvar*BDI_s4   
BDI_s8~~errvar*BDI_s8   
BDI_s12~~errvar*BDI_s12
BDI_s16~~errvar*BDI_s16

'
fit_model221_LAT <- growth(model221_LAT, data = med_master_v5)
summary(fit_model221_LAT, standardized = TRUE)
fitmeasures(fit_model221_LAT)


?rep
###MAKING BDI FILE 

getwd()
library(plyr)


patnos <- read.csv("patnos.csv", stringsAsFactors = FALSE )
View(patnos)
patnos_rep <- data.frame(patnos[rep(seq_len(nrow(patnos)), each=23),])
names(patnos_rep)[1] <- "patnos"
View(patnos_rep)


