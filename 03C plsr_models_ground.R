setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(pls)
library(ggplot2)
library(caret)
library(reshape2)
library(RColorBrewer)
library(patchwork)

########################################
## to do:
## try Model II regressions?

#########################################
## define functions

## root mean squared deviation
RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
}

## percent RMSD (based on data quantiles)
## set min and max to 0 and 1 for range as denominator
## or to 0.25 and 0.75 for IQR as denominator
percentRMSD<-function(measured,predicted,min,max,na.rm=T){
  RMSD_data<-RMSD(measured,predicted)
  range<-unname(quantile(measured,probs=max,na.rm=na.rm)-quantile(measured,probs=min,na.rm=na.rm))
  return(RMSD_data/range)
}

## applying coefficients to validation spectra
apply.coefs<-function(coef.list,val.spec,intercept=T){
  if(sum(lapply(coef.list,length)==ncol(val.spec)+intercept) < length(coef.list)){
    stop("some coefficients have the wrong length")
  }
  
  coef.matrix<-matrix(unlist(coef.list),
                      nrow=length(coef.list),
                      byrow=T)
  
  if(intercept==T){
    pred.matrix<-t(t(as.matrix(val.spec) %*% t(coef.matrix[,-1]))+coef.matrix[,1])
  } else {
    pred.matrix<-as.matrix(val.spec) %*% t(coef.matrix)
  }
}

######################################################
## read data

ground_spec_EL_agg_train<-readRDS("ProcessedSpectralData/ground_spec_EL_agg_train.rds")
ground_spec_EL_agg_test<-readRDS("ProcessedSpectralData/ground_spec_EL_agg_test.rds")

####################################################
## building calibration models

perC_ground<-plsr(meta(ground_spec_EL_agg_train)$C~as.matrix(ground_spec_EL_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perC_ground <- selectNcomp(perC_ground, method = "onesigma", plot = FALSE)
perC_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$C))
perC_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[perC_ground_valid],
                             Species=meta(ground_spec_EL_agg_train)$Species[perC_ground_valid],
                             Project=meta(ground_spec_EL_agg_train)$Project[perC_ground_valid],
                             Stage=meta(ground_spec_EL_agg_train)$Stage[perC_ground_valid],
                             GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[perC_ground_valid],
                             measured=meta(ground_spec_EL_agg_train)$C[perC_ground_valid],
                             val_pred=perC_ground$validation$pred[,,ncomp_perC_ground])
ggplot(perC_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %C from ground-leaf spectra")

perN_ground<-plsr(meta(ground_spec_EL_agg_train)$N~as.matrix(ground_spec_EL_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perN_ground <- selectNcomp(perN_ground, method = "onesigma", plot = FALSE)
perN_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$N))
perN_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[perN_ground_valid],
                             Species=meta(ground_spec_EL_agg_train)$Species[perN_ground_valid],
                             Project=meta(ground_spec_EL_agg_train)$Project[perN_ground_valid],
                             Stage=meta(ground_spec_EL_agg_train)$Stage[perN_ground_valid],
                             GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[perN_ground_valid],
                             measured=meta(ground_spec_EL_agg_train)$N[perN_ground_valid],
                             val_pred=perN_ground$validation$pred[,,ncomp_perN_ground])
ggplot(perN_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %N from ground-leaf spectra")+guides(color=F)

LMA_ground<-plsr(meta(ground_spec_EL_agg_train)$LMA~as.matrix(ground_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_ground <- selectNcomp(LMA_ground, method = "onesigma", plot = FALSE)
LMA_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$LMA))
LMA_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[LMA_ground_valid],
                            Species=meta(ground_spec_EL_agg_train)$Species[LMA_ground_valid],
                            Project=meta(ground_spec_EL_agg_train)$Project[LMA_ground_valid],
                            Stage=meta(ground_spec_EL_agg_train)$Stage[LMA_ground_valid],
                            GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[LMA_ground_valid],
                            measured=meta(ground_spec_EL_agg_train)$LMA[LMA_ground_valid],
                            val_pred=LMA_ground$validation$pred[,,ncomp_LMA_ground])
ggplot(LMA_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.25))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting LMA from ground-leaf spectra")

LDMC_ground<-plsr(meta(ground_spec_EL_agg_train)$LDMC~as.matrix(ground_spec_EL_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_ground <- selectNcomp(LDMC_ground, method = "onesigma", plot = FALSE)
LDMC_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$LDMC))
LDMC_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[LDMC_ground_valid],
                             Species=meta(ground_spec_EL_agg_train)$Species[LDMC_ground_valid],
                             Project=meta(ground_spec_EL_agg_train)$Project[LDMC_ground_valid],
                             Stage=meta(ground_spec_EL_agg_train)$Stage[LDMC_ground_valid],
                             GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[LDMC_ground_valid],
                             measured=meta(ground_spec_EL_agg_train)$LDMC[LDMC_ground_valid],
                             val_pred=LDMC_ground$validation$pred[,,ncomp_LDMC_ground])
ggplot(LDMC_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(150,600),ylim=c(150,600))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting LDMC from ground-leaf spectra")

EWT_ground<-plsr(meta(ground_spec_EL_agg_train)$EWT~as.matrix(ground_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_ground <- selectNcomp(EWT_ground, method = "onesigma", plot = FALSE)
EWT_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$EWT))
EWT_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[EWT_ground_valid],
                            Species=meta(ground_spec_EL_agg_train)$Species[EWT_ground_valid],
                            Project=meta(ground_spec_EL_agg_train)$Project[EWT_ground_valid],
                            Stage=meta(ground_spec_EL_agg_train)$Stage[EWT_ground_valid],
                            GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[EWT_ground_valid],
                            measured=meta(ground_spec_EL_agg_train)$EWT[EWT_ground_valid],
                            val_pred=EWT_ground$validation$pred[,,ncomp_EWT_ground])
ggplot(EWT_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.03),ylim=c(0,0.03))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting EWT from ground-leaf spectra")

chlA_ground<-plsr(meta(ground_spec_EL_agg_train)$chlA~as.matrix(ground_spec_EL_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_ground <- selectNcomp(chlA_ground, method = "onesigma", plot = FALSE)
chlA_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$chlA))
chlA_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[chlA_ground_valid],
                             Species=meta(ground_spec_EL_agg_train)$Species[chlA_ground_valid],
                             Project=meta(ground_spec_EL_agg_train)$Project[chlA_ground_valid],
                             Stage=meta(ground_spec_EL_agg_train)$Stage[chlA_ground_valid],
                             GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[chlA_ground_valid],
                             measured=meta(ground_spec_EL_agg_train)$chlA[chlA_ground_valid],
                             val_pred=chlA_ground$validation$pred[,,ncomp_chlA_ground])
ggplot(chlA_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,16),ylim=c(0,16))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl a from ground-leaf spectra")

chlB_ground<-plsr(meta(ground_spec_EL_agg_train)$chlB~as.matrix(ground_spec_EL_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_ground <- selectNcomp(chlB_ground, method = "onesigma", plot = FALSE)
chlB_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$chlB))
chlB_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[chlB_ground_valid],
                             Species=meta(ground_spec_EL_agg_train)$Species[chlB_ground_valid],
                             Project=meta(ground_spec_EL_agg_train)$Project[chlB_ground_valid],
                             Stage=meta(ground_spec_EL_agg_train)$Stage[chlB_ground_valid],
                             GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[chlB_ground_valid],
                             measured=meta(ground_spec_EL_agg_train)$chlB[chlB_ground_valid],
                             val_pred=chlB_ground$validation$pred[,,ncomp_chlB_ground])
ggplot(chlB_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl b from ground-leaf spectra")

car_ground<-plsr(meta(ground_spec_EL_agg_train)$car~as.matrix(ground_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_ground <- selectNcomp(car_ground, method = "onesigma", plot = FALSE)
car_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$car))
car_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[car_ground_valid],
                            Species=meta(ground_spec_EL_agg_train)$Species[car_ground_valid],
                            Project=meta(ground_spec_EL_agg_train)$Project[car_ground_valid],
                            Stage=meta(ground_spec_EL_agg_train)$Stage[car_ground_valid],
                            GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[car_ground_valid],
                            measured=meta(ground_spec_EL_agg_train)$car[car_ground_valid],
                            val_pred=car_ground$validation$pred[,,ncomp_car_ground])
ggplot(car_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,3.5),ylim=c(0,3.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting total car from ground-leaf spectra")

solubles_ground<-plsr(meta(ground_spec_EL_agg_train)$solubles~as.matrix(ground_spec_EL_agg_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_ground <- selectNcomp(solubles_ground, method = "onesigma", plot = FALSE)
solubles_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$solubles))
solubles_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[solubles_ground_valid],
                                 Species=meta(ground_spec_EL_agg_train)$Species[solubles_ground_valid],
                                 Project=meta(ground_spec_EL_agg_train)$Project[solubles_ground_valid],
                                 Stage=meta(ground_spec_EL_agg_train)$Stage[solubles_ground_valid],
                                 GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[solubles_ground_valid],
                                 measured=meta(ground_spec_EL_agg_train)$solubles[solubles_ground_valid],
                                 val_pred=solubles_ground$validation$pred[,,ncomp_solubles_ground])
ggplot(solubles_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(35,90),ylim=c(35,90))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting solubles from ground-leaf spectra")

hemicellulose_ground<-plsr(meta(ground_spec_EL_agg_train)$hemicellulose~as.matrix(ground_spec_EL_agg_train),
                           ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_ground <- selectNcomp(hemicellulose_ground, method = "onesigma", plot = FALSE)
hemicellulose_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$hemicellulose))
hemicellulose_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[hemicellulose_ground_valid],
                                      Species=meta(ground_spec_EL_agg_train)$Species[hemicellulose_ground_valid],
                                      Project=meta(ground_spec_EL_agg_train)$Project[hemicellulose_ground_valid],
                                      Stage=meta(ground_spec_EL_agg_train)$Stage[hemicellulose_ground_valid],
                                      GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[hemicellulose_ground_valid],
                                      measured=meta(ground_spec_EL_agg_train)$hemicellulose[hemicellulose_ground_valid],
                                      val_pred=hemicellulose_ground$validation$pred[,,ncomp_hemicellulose_ground])
ggplot(hemicellulose_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,36),ylim=c(0,36))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting hemicellulose from ground-leaf spectra")

cellulose_ground<-plsr(meta(ground_spec_EL_agg_train)$cellulose~as.matrix(ground_spec_EL_agg_train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_ground <- selectNcomp(cellulose_ground, method = "onesigma", plot = FALSE)
cellulose_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$cellulose))
cellulose_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[cellulose_ground_valid],
                                  Species=meta(ground_spec_EL_agg_train)$Species[cellulose_ground_valid],
                                  Project=meta(ground_spec_EL_agg_train)$Project[cellulose_ground_valid],
                                  Stage=meta(ground_spec_EL_agg_train)$Stage[cellulose_ground_valid],
                                  GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[cellulose_ground_valid],
                                  measured=meta(ground_spec_EL_agg_train)$cellulose[cellulose_ground_valid],
                                  val_pred=cellulose_ground$validation$pred[,,ncomp_cellulose_ground])
ggplot(cellulose_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(2,28),ylim=c(2,28))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting cellulose from ground-leaf spectra")

lignin_ground<-plsr(meta(ground_spec_EL_agg_train)$lignin~as.matrix(ground_spec_EL_agg_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_ground <- selectNcomp(lignin_ground, method = "onesigma", plot = FALSE)
lignin_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$lignin))
lignin_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[lignin_ground_valid],
                               Species=meta(ground_spec_EL_agg_train)$Species[lignin_ground_valid],
                               Project=meta(ground_spec_EL_agg_train)$Project[lignin_ground_valid],
                               Stage=meta(ground_spec_EL_agg_train)$Stage[lignin_ground_valid],
                               GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[lignin_ground_valid],
                               measured=meta(ground_spec_EL_agg_train)$lignin[lignin_ground_valid],
                               val_pred=lignin_ground$validation$pred[,,ncomp_lignin_ground])
ggplot(lignin_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,22),ylim=c(0,22))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting lignin from ground-leaf spectra")


Al_ground<-plsr(meta(ground_spec_EL_agg_train)$Al~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_ground <- selectNcomp(Al_ground, method = "onesigma", plot = FALSE)
Al_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Al))
Al_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Al_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Al_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Al_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Al_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Al_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Al[Al_ground_valid],
                           val_pred=Al_ground$validation$pred[,,ncomp_Al_ground])
ggplot(Al_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.05,0.35),ylim=c(-0.05,0.35))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Al from ground-leaf spectra")

Ca_ground<-plsr(meta(ground_spec_EL_agg_train)$Ca~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_ground <- selectNcomp(Ca_ground, method = "onesigma", plot = FALSE)
Ca_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Ca))
Ca_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Ca_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Ca_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Ca_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Ca_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Ca_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Ca[Ca_ground_valid],
                           val_pred=Ca_ground$validation$pred[,,ncomp_Ca_ground])
ggplot(Ca_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-10,40),ylim=c(-10,40))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Ca from ground-leaf spectra")

Cu_ground<-plsr(meta(ground_spec_EL_agg_train)$Cu~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_ground <- selectNcomp(Cu_ground, method = "onesigma", plot = FALSE)
Cu_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Cu))
Cu_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Cu_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Cu_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Cu_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Cu_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Cu_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Cu[Cu_ground_valid],
                           val_pred=Cu_ground$validation$pred[,,ncomp_Cu_ground])
ggplot(Cu_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.005,0.055),ylim=c(-0.005,0.055))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Cu from ground-leaf spectra")

Fe_ground<-plsr(meta(ground_spec_EL_agg_train)$Fe~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_ground <- selectNcomp(Fe_ground, method = "onesigma", plot = FALSE)
Fe_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Fe))
Fe_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Fe_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Fe_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Fe_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Fe_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Fe_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Fe[Fe_ground_valid],
                           val_pred=Fe_ground$validation$pred[,,ncomp_Fe_ground])
ggplot(Fe_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Fe from ground-leaf spectra")

K_ground<-plsr(meta(ground_spec_EL_agg_train)$K~as.matrix(ground_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_ground <- selectNcomp(K_ground, method = "onesigma", plot = FALSE)
K_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$K))
K_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[K_ground_valid],
                          Species=meta(ground_spec_EL_agg_train)$Species[K_ground_valid],
                          Project=meta(ground_spec_EL_agg_train)$Project[K_ground_valid],
                          Stage=meta(ground_spec_EL_agg_train)$Stage[K_ground_valid],
                          GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[K_ground_valid],
                          measured=meta(ground_spec_EL_agg_train)$K[K_ground_valid],
                          val_pred=K_ground$validation$pred[,,ncomp_K_ground])
ggplot(K_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,35),ylim=c(0,35))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting K from ground-leaf spectra")

Mg_ground<-plsr(meta(ground_spec_EL_agg_train)$Mg~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_ground <- selectNcomp(Mg_ground, method = "onesigma", plot = FALSE)
Mg_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Mg))
Mg_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Mg_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Mg_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Mg_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Mg_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Mg_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Mg[Mg_ground_valid],
                           val_pred=Mg_ground$validation$pred[,,ncomp_Mg_ground])
ggplot(Mg_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,8),ylim=c(0,8))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mg from ground-leaf spectra")

Mn_ground<-plsr(meta(ground_spec_EL_agg_train)$Mn~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_ground <- selectNcomp(Mn_ground, method = "onesigma", plot = FALSE)
Mn_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Mn))
Mn_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Mn_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Mn_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Mn_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Mn_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Mn_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Mn[Mn_ground_valid],
                           val_pred=Mn_ground$validation$pred[,,ncomp_Mn_ground])
ggplot(Mn_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,1.1),ylim=c(-0.1,1.1))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mn from ground-leaf spectra")

Na_ground<-plsr(meta(ground_spec_EL_agg_train)$Na~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_ground <- selectNcomp(Na_ground, method = "onesigma", plot = FALSE)
Na_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Na))
Na_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Na_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Na_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Na_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Na_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Na_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Na[Na_ground_valid],
                           val_pred=Na_ground$validation$pred[,,ncomp_Na_ground])
ggplot(Na_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,5),ylim=c(-0.5,5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Na from ground-leaf spectra")

P_ground<-plsr(meta(ground_spec_EL_agg_train)$P~as.matrix(ground_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_ground <- selectNcomp(P_ground, method = "onesigma", plot = FALSE)
P_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$P))
P_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[P_ground_valid],
                          Species=meta(ground_spec_EL_agg_train)$Species[P_ground_valid],
                          Project=meta(ground_spec_EL_agg_train)$Project[P_ground_valid],
                          Stage=meta(ground_spec_EL_agg_train)$Stage[P_ground_valid],
                          GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[P_ground_valid],
                          measured=meta(ground_spec_EL_agg_train)$P[P_ground_valid],
                          val_pred=P_ground$validation$pred[,,ncomp_P_ground])
ggplot(P_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,8),ylim=c(0,8))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting P from ground-leaf spectra")

Zn_ground<-plsr(meta(ground_spec_EL_agg_train)$Zn~as.matrix(ground_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_ground <- selectNcomp(Zn_ground, method = "onesigma", plot = FALSE)
Zn_ground_valid <- which(!is.na(meta(ground_spec_EL_agg_train)$Zn))
Zn_ground_pred<-data.frame(ID=meta(ground_spec_EL_agg_train)$ID[Zn_ground_valid],
                           Species=meta(ground_spec_EL_agg_train)$Species[Zn_ground_valid],
                           Project=meta(ground_spec_EL_agg_train)$Project[Zn_ground_valid],
                           Stage=meta(ground_spec_EL_agg_train)$Stage[Zn_ground_valid],
                           GrowthForm=meta(ground_spec_EL_agg_train)$GrowthForm[Zn_ground_valid],
                           measured=meta(ground_spec_EL_agg_train)$Zn[Zn_ground_valid],
                           val_pred=Zn_ground$validation$pred[,,ncomp_Zn_ground])
ggplot(Zn_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,0.7),ylim=c(-0.1,0.7))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Zn from ground-leaf spectra")

###############################################
## VIP plots

source("VIP.R")

VIP_ground<-data.frame(LMA=VIP(LMA_ground)[ncomp_LMA_ground,],
                      LDMC=VIP(LDMC_ground)[ncomp_LDMC_ground,],
                      EWT=VIP(EWT_ground)[ncomp_EWT_ground,],
                      sol=VIP(solubles_ground)[ncomp_solubles_ground,],
                      hemi=VIP(hemicellulose_ground)[ncomp_hemicellulose_ground,],
                      cell=VIP(cellulose_ground)[ncomp_cellulose_ground,],
                      lign=VIP(lignin_ground)[ncomp_lignin_ground,],
                      chlA=VIP(chlA_ground)[ncomp_chlA_ground,],
                      chlB=VIP(chlB_ground)[ncomp_chlB_ground,],
                      car=VIP(car_ground)[ncomp_car_ground,],
                      C=VIP(perC_ground)[ncomp_perC_ground,],
                      N=VIP(perN_ground)[ncomp_perN_ground,],
                      Al=VIP(Al_ground)[ncomp_Al_ground,],
                      Ca=VIP(Ca_ground)[ncomp_Ca_ground,],
                      Cu=VIP(Cu_ground)[ncomp_Cu_ground,],
                      Fe=VIP(Fe_ground)[ncomp_Fe_ground,],
                      K=VIP(K_ground)[ncomp_K_ground,],
                      Mg=VIP(Mg_ground)[ncomp_Mg_ground,],
                      Mn=VIP(Mn_ground)[ncomp_Mn_ground,],
                      Na=VIP(Na_ground)[ncomp_Na_ground,],
                      P=VIP(P_ground)[ncomp_P_ground,],
                      Zn=VIP(Zn_ground)[ncomp_Zn_ground,],
                      wavelength=400:2400)

saveRDS(VIP_ground,"SavedResults/VIP_ground.rds")

###############################################
## validation: ground leaves

solubles_jack_coefs_ground<-list()
hemicellulose_jack_coefs_ground<-list()
cellulose_jack_coefs_ground<-list()
lignin_jack_coefs_ground<-list()
perC_jack_coefs_ground<-list()
perN_jack_coefs_ground<-list()
LMA_jack_coefs_ground<-list()
LDMC_jack_coefs_ground<-list()
EWT_jack_coefs_ground<-list()
chlA_jack_coefs_ground<-list()
chlB_jack_coefs_ground<-list()
car_jack_coefs_ground<-list()
Al_jack_coefs_ground<-list()
Ca_jack_coefs_ground<-list()
Cu_jack_coefs_ground<-list()
Fe_jack_coefs_ground<-list()
K_jack_coefs_ground<-list()
Mg_jack_coefs_ground<-list()
Mn_jack_coefs_ground<-list()
Na_jack_coefs_ground<-list()
P_jack_coefs_ground<-list()
Zn_jack_coefs_ground<-list()

solubles_jack_stats_ground<-list()
hemicellulose_jack_stats_ground<-list()
cellulose_jack_stats_ground<-list()
lignin_jack_stats_ground<-list()
perC_jack_stats_ground<-list()
perN_jack_stats_ground<-list()
LMA_jack_stats_ground<-list()
LDMC_jack_stats_ground<-list()
EWT_jack_stats_ground<-list()
chlA_jack_stats_ground<-list()
chlB_jack_stats_ground<-list()
car_jack_stats_ground<-list()
Al_jack_stats_ground<-list()
Ca_jack_stats_ground<-list()
Cu_jack_stats_ground<-list()
Fe_jack_stats_ground<-list()
K_jack_stats_ground<-list()
Mg_jack_stats_ground<-list()
Mn_jack_stats_ground<-list()
Na_jack_stats_ground<-list()
P_jack_stats_ground<-list()
Zn_jack_stats_ground<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_ground<-nrow(ground_spec_EL_agg_train)
  train_jack_ground<-sample(1:n_cal_spec_ground,floor(0.7*n_cal_spec_ground))
  test_jack_ground<-setdiff(1:n_cal_spec_ground,train_jack_ground)
  
  calib_jack_ground<-ground_spec_EL_agg_train[train_jack_ground]
  val_jack_ground<-ground_spec_EL_agg_train[test_jack_ground]
  
  solubles_ground_jack<-plsr(meta(calib_jack_ground)$solubles~as.matrix(calib_jack_ground),
                             ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_ground_jack<-plsr(meta(calib_jack_ground)$hemicellulose~as.matrix(calib_jack_ground),
                                  ncomp=30,method = "oscorespls",validation="none")
  cellulose_ground_jack<-plsr(meta(calib_jack_ground)$cellulose~as.matrix(calib_jack_ground),
                              ncomp=30,method = "oscorespls",validation="none")
  lignin_ground_jack<-plsr(meta(calib_jack_ground)$lignin~as.matrix(calib_jack_ground),
                           ncomp=30,method = "oscorespls",validation="none")
  perC_ground_jack<-plsr(meta(calib_jack_ground)$C~as.matrix(calib_jack_ground),
                         ncomp=30,method = "oscorespls",validation="none")
  perN_ground_jack<-plsr(meta(calib_jack_ground)$N~as.matrix(calib_jack_ground),
                         ncomp=30,method = "oscorespls",validation="none")
  LMA_ground_jack<-plsr(meta(calib_jack_ground)$LMA~as.matrix(calib_jack_ground),
                        ncomp=30,method = "oscorespls",validation="none")
  LDMC_ground_jack<-plsr(meta(calib_jack_ground)$LDMC~as.matrix(calib_jack_ground),
                         ncomp=30,method = "oscorespls",validation="none")
  EWT_ground_jack<-plsr(meta(calib_jack_ground)$EWT~as.matrix(calib_jack_ground),
                        ncomp=30,method = "oscorespls",validation="none")
  chlA_ground_jack<-plsr(meta(calib_jack_ground)$chlA~as.matrix(calib_jack_ground),
                         ncomp=30,method = "oscorespls",validation="none")
  chlB_ground_jack<-plsr(meta(calib_jack_ground)$chlB~as.matrix(calib_jack_ground),
                         ncomp=30,method = "oscorespls",validation="none")
  car_ground_jack<-plsr(meta(calib_jack_ground)$car~as.matrix(calib_jack_ground),
                        ncomp=30,method = "oscorespls",validation="none")
  Al_ground_jack<-plsr(meta(calib_jack_ground)$Al~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  Ca_ground_jack<-plsr(meta(calib_jack_ground)$Ca~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  Cu_ground_jack<-plsr(meta(calib_jack_ground)$Cu~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  Fe_ground_jack<-plsr(meta(calib_jack_ground)$Fe~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  K_ground_jack<-plsr(meta(calib_jack_ground)$K~as.matrix(calib_jack_ground),
                      ncomp=30,method = "oscorespls",validation="none")
  Mg_ground_jack<-plsr(meta(calib_jack_ground)$Mg~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  Mn_ground_jack<-plsr(meta(calib_jack_ground)$Mn~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  Na_ground_jack<-plsr(meta(calib_jack_ground)$Na~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  P_ground_jack<-plsr(meta(calib_jack_ground)$P~as.matrix(calib_jack_ground),
                      ncomp=30,method = "oscorespls",validation="none")
  Zn_ground_jack<-plsr(meta(calib_jack_ground)$Zn~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  
  solubles_jack_val_pred_ground<-as.vector(predict(solubles_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_solubles_ground)[,,1])
  solubles_jack_val_fit_ground<-lm(meta(val_jack_ground)$solubles~solubles_jack_val_pred_ground)
  solubles_jack_stats_ground[[i]]<-c(R2=summary(solubles_jack_val_fit_ground)$r.squared,
                                     RMSE=RMSD(meta(val_jack_ground)$solubles,solubles_jack_val_pred_ground),
                                     perRMSE=percentRMSD(meta(val_jack_ground)$solubles,solubles_jack_val_pred_ground,0.025,0.975),
                                     bias=mean(solubles_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$solubles,na.rm=T))
  
  hemicellulose_jack_val_pred_ground<-as.vector(predict(hemicellulose_ground_jack,newdata=as.matrix(val_jack_ground),
                                                        ncomp=ncomp_hemicellulose_ground)[,,1])
  hemicellulose_jack_val_fit_ground<-lm(meta(val_jack_ground)$hemicellulose~hemicellulose_jack_val_pred_ground)
  hemicellulose_jack_stats_ground[[i]]<-c(R2=summary(hemicellulose_jack_val_fit_ground)$r.squared,
                                          RMSE=RMSD(meta(val_jack_ground)$hemicellulose,hemicellulose_jack_val_pred_ground),
                                          perRMSE=percentRMSD(meta(val_jack_ground)$hemicellulose,hemicellulose_jack_val_pred_ground,0.025,0.975),
                                          bias=mean(hemicellulose_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$hemicellulose,na.rm=T))
  
  cellulose_jack_val_pred_ground<-as.vector(predict(cellulose_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_cellulose_ground)[,,1])
  cellulose_jack_val_fit_ground<-lm(meta(val_jack_ground)$cellulose~cellulose_jack_val_pred_ground)
  cellulose_jack_stats_ground[[i]]<-c(R2=summary(cellulose_jack_val_fit_ground)$r.squared,
                                      RMSE=RMSD(meta(val_jack_ground)$cellulose,cellulose_jack_val_pred_ground),
                                      perRMSE=percentRMSD(meta(val_jack_ground)$cellulose,cellulose_jack_val_pred_ground,0.025,0.975),
                                      bias=mean(cellulose_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$cellulose,na.rm=T))
  
  lignin_jack_val_pred_ground<-as.vector(predict(lignin_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_lignin_ground)[,,1])
  lignin_jack_val_fit_ground<-lm(meta(val_jack_ground)$lignin~lignin_jack_val_pred_ground)
  lignin_jack_stats_ground[[i]]<-c(R2=summary(lignin_jack_val_fit_ground)$r.squared,
                                   RMSE=RMSD(meta(val_jack_ground)$lignin,lignin_jack_val_pred_ground),
                                   perRMSE=percentRMSD(meta(val_jack_ground)$lignin,lignin_jack_val_pred_ground,0.025,0.975),
                                   bias=mean(lignin_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$lignin,na.rm=T))
  
  perC_jack_val_pred_ground<-as.vector(predict(perC_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_perC_ground)[,,1])
  perC_jack_val_fit_ground<-lm(meta(val_jack_ground)$C~perC_jack_val_pred_ground)
  perC_jack_stats_ground[[i]]<-c(R2=summary(perC_jack_val_fit_ground)$r.squared,
                                 RMSE=RMSD(meta(val_jack_ground)$C,perC_jack_val_pred_ground),
                                 perRMSE=percentRMSD(meta(val_jack_ground)$C,perC_jack_val_pred_ground,0.025,0.975),
                                 bias=mean(perC_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$C,na.rm=T))
  
  perN_jack_val_pred_ground<-as.vector(predict(perN_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_perN_ground)[,,1])
  perN_jack_val_fit_ground<-lm(meta(val_jack_ground)$N~perN_jack_val_pred_ground)
  perN_jack_stats_ground[[i]]<-c(R2=summary(perN_jack_val_fit_ground)$r.squared,
                                 RMSE=RMSD(meta(val_jack_ground)$N,perN_jack_val_pred_ground),
                                 perRMSE=percentRMSD(meta(val_jack_ground)$N,perN_jack_val_pred_ground,0.025,0.975),
                                 bias=mean(perN_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$N,na.rm=T))
  
  LMA_jack_val_pred_ground<-as.vector(predict(LMA_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_LMA_ground)[,,1])
  LMA_jack_val_fit_ground<-lm(meta(val_jack_ground)$LMA~LMA_jack_val_pred_ground)
  LMA_jack_stats_ground[[i]]<-c(R2=summary(LMA_jack_val_fit_ground)$r.squared,
                                RMSE=RMSD(meta(val_jack_ground)$LMA,LMA_jack_val_pred_ground),
                                perRMSE=percentRMSD(meta(val_jack_ground)$LMA,LMA_jack_val_pred_ground,0.025,0.975),
                                bias=mean(LMA_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$LMA,na.rm=T))
  
  LDMC_jack_val_pred_ground<-as.vector(predict(LDMC_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_LDMC_ground)[,,1])
  LDMC_jack_val_fit_ground<-lm(meta(val_jack_ground)$LDMC~LDMC_jack_val_pred_ground)
  LDMC_jack_stats_ground[[i]]<-c(R2=summary(LDMC_jack_val_fit_ground)$r.squared,
                                 RMSE=RMSD(meta(val_jack_ground)$LDMC,LDMC_jack_val_pred_ground),
                                 perRMSE=percentRMSD(meta(val_jack_ground)$LDMC,LDMC_jack_val_pred_ground,0.025,0.975),
                                 bias=mean(LDMC_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$LDMC,na.rm=T))
  
  EWT_jack_val_pred_ground<-as.vector(predict(EWT_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_EWT_ground)[,,1])
  EWT_jack_val_fit_ground<-lm(meta(val_jack_ground)$EWT~EWT_jack_val_pred_ground)
  EWT_jack_stats_ground[[i]]<-c(R2=summary(EWT_jack_val_fit_ground)$r.squared,
                                RMSE=RMSD(meta(val_jack_ground)$EWT,EWT_jack_val_pred_ground),
                                perRMSE=percentRMSD(meta(val_jack_ground)$EWT,EWT_jack_val_pred_ground,0.025,0.975),
                                bias=mean(EWT_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$EWT,na.rm=T))
  
  chlA_jack_val_pred_ground<-as.vector(predict(chlA_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_chlA_ground)[,,1])
  chlA_jack_val_fit_ground<-lm(meta(val_jack_ground)$chlA~chlA_jack_val_pred_ground)
  chlA_jack_stats_ground[[i]]<-c(R2=summary(chlA_jack_val_fit_ground)$r.squared,
                                 RMSE=RMSD(meta(val_jack_ground)$chlA,chlA_jack_val_pred_ground),
                                 perRMSE=percentRMSD(meta(val_jack_ground)$chlA,chlA_jack_val_pred_ground,0.025,0.975),
                                 bias=mean(chlA_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$chlA,na.rm=T))
  
  chlB_jack_val_pred_ground<-as.vector(predict(chlB_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_chlB_ground)[,,1])
  chlB_jack_val_fit_ground<-lm(meta(val_jack_ground)$chlB~chlB_jack_val_pred_ground)
  chlB_jack_stats_ground[[i]]<-c(R2=summary(chlB_jack_val_fit_ground)$r.squared,
                                 RMSE=RMSD(meta(val_jack_ground)$chlB,chlB_jack_val_pred_ground),
                                 perRMSE=percentRMSD(meta(val_jack_ground)$chlB,chlB_jack_val_pred_ground,0.025,0.975),
                                 bias=mean(chlB_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$chlB,na.rm=T))
  
  car_jack_val_pred_ground<-as.vector(predict(car_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_car_ground)[,,1])
  car_jack_val_fit_ground<-lm(meta(val_jack_ground)$car~car_jack_val_pred_ground)
  car_jack_stats_ground[[i]]<-c(R2=summary(car_jack_val_fit_ground)$r.squared,
                                RMSE=RMSD(meta(val_jack_ground)$car,car_jack_val_pred_ground),
                                perRMSE=percentRMSD(meta(val_jack_ground)$car,car_jack_val_pred_ground,0.025,0.975),
                                bias=mean(car_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$car,na.rm=T))
  
  Al_jack_val_pred_ground<-as.vector(predict(Al_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Al_ground)[,,1])
  Al_jack_val_fit_ground<-lm(meta(val_jack_ground)$Al~Al_jack_val_pred_ground)
  Al_jack_stats_ground[[i]]<-c(R2=summary(Al_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Al,Al_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Al,Al_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Al_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Al,na.rm=T))
  
  Ca_jack_val_pred_ground<-as.vector(predict(Ca_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Ca_ground)[,,1])
  Ca_jack_val_fit_ground<-lm(meta(val_jack_ground)$Ca~Ca_jack_val_pred_ground)
  Ca_jack_stats_ground[[i]]<-c(R2=summary(Ca_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Ca,Ca_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Ca,Ca_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Ca_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Ca,na.rm=T))
  
  Cu_jack_val_pred_ground<-as.vector(predict(Cu_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Cu_ground)[,,1])
  Cu_jack_val_fit_ground<-lm(meta(val_jack_ground)$Cu~Cu_jack_val_pred_ground)
  Cu_jack_stats_ground[[i]]<-c(R2=summary(Cu_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Cu,Cu_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Cu,Cu_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Cu_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Cu,na.rm=T))
  
  Fe_jack_val_pred_ground<-as.vector(predict(Fe_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Fe_ground)[,,1])
  Fe_jack_val_fit_ground<-lm(meta(val_jack_ground)$Fe~Fe_jack_val_pred_ground)
  Fe_jack_stats_ground[[i]]<-c(R2=summary(Fe_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Fe,Fe_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Fe,Fe_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Fe_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Fe,na.rm=T))
  
  K_jack_val_pred_ground<-as.vector(predict(K_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_K_ground)[,,1])
  K_jack_val_fit_ground<-lm(meta(val_jack_ground)$K~K_jack_val_pred_ground)
  K_jack_stats_ground[[i]]<-c(R2=summary(K_jack_val_fit_ground)$r.squared,
                              RMSE=RMSD(meta(val_jack_ground)$K,K_jack_val_pred_ground),
                              perRMSE=percentRMSD(meta(val_jack_ground)$K,K_jack_val_pred_ground,0.025,0.975),
                              bias=mean(K_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$K,na.rm=T))
  
  Mg_jack_val_pred_ground<-as.vector(predict(Mg_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Mg_ground)[,,1])
  Mg_jack_val_fit_ground<-lm(meta(val_jack_ground)$Mg~Mg_jack_val_pred_ground)
  Mg_jack_stats_ground[[i]]<-c(R2=summary(Mg_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Mg,Mg_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Mg,Mg_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Mg_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Mg,na.rm=T))
  
  Mn_jack_val_pred_ground<-as.vector(predict(Mn_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Mn_ground)[,,1])
  Mn_jack_val_fit_ground<-lm(meta(val_jack_ground)$Mn~Mn_jack_val_pred_ground)
  Mn_jack_stats_ground[[i]]<-c(R2=summary(Mn_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Mn,Mn_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Mn,Mn_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Mn_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Mn,na.rm=T))
  
  Na_jack_val_pred_ground<-as.vector(predict(Na_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Na_ground)[,,1])
  Na_jack_val_fit_ground<-lm(meta(val_jack_ground)$Na~Na_jack_val_pred_ground)
  Na_jack_stats_ground[[i]]<-c(R2=summary(Na_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Na,Na_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Na,Na_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Na_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Na,na.rm=T))
  
  P_jack_val_pred_ground<-as.vector(predict(P_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_P_ground)[,,1])
  P_jack_val_fit_ground<-lm(meta(val_jack_ground)$P~P_jack_val_pred_ground)
  P_jack_stats_ground[[i]]<-c(R2=summary(P_jack_val_fit_ground)$r.squared,
                              RMSE=RMSD(meta(val_jack_ground)$P,P_jack_val_pred_ground),
                              perRMSE=percentRMSD(meta(val_jack_ground)$P,P_jack_val_pred_ground,0.025,0.975),
                              bias=mean(P_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$P,na.rm=T))
  
  Zn_jack_val_pred_ground<-as.vector(predict(Zn_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_Zn_ground)[,,1])
  Zn_jack_val_fit_ground<-lm(meta(val_jack_ground)$Zn~Zn_jack_val_pred_ground)
  Zn_jack_stats_ground[[i]]<-c(R2=summary(Zn_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$Zn,Zn_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$Zn,Zn_jack_val_pred_ground,0.025,0.975),
                               bias=mean(Zn_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$Zn,na.rm=T))
  
  solubles_jack_coefs_ground[[i]]<-as.vector(coef(solubles_ground_jack,ncomp=ncomp_solubles_ground,intercept=TRUE))
  hemicellulose_jack_coefs_ground[[i]]<-as.vector(coef(hemicellulose_ground_jack,ncomp=ncomp_hemicellulose_ground,intercept=TRUE))
  cellulose_jack_coefs_ground[[i]]<-as.vector(coef(cellulose_ground_jack,ncomp=ncomp_cellulose_ground,intercept=TRUE))
  lignin_jack_coefs_ground[[i]]<-as.vector(coef(lignin_ground_jack,ncomp=ncomp_lignin_ground,intercept=TRUE))
  perC_jack_coefs_ground[[i]]<-as.vector(coef(perC_ground_jack,ncomp=ncomp_perC_ground,intercept=TRUE))
  perN_jack_coefs_ground[[i]]<-as.vector(coef(perN_ground_jack,ncomp=ncomp_perN_ground,intercept=TRUE))
  LMA_jack_coefs_ground[[i]]<-as.vector(coef(LMA_ground_jack,ncomp=ncomp_LMA_ground,intercept=TRUE))
  LDMC_jack_coefs_ground[[i]]<-as.vector(coef(LDMC_ground_jack,ncomp=ncomp_LDMC_ground,intercept=TRUE))
  EWT_jack_coefs_ground[[i]]<-as.vector(coef(EWT_ground_jack,ncomp=ncomp_EWT_ground,intercept=TRUE))
  chlA_jack_coefs_ground[[i]]<-as.vector(coef(chlA_ground_jack,ncomp=ncomp_chlA_ground,intercept=TRUE))
  chlB_jack_coefs_ground[[i]]<-as.vector(coef(chlB_ground_jack,ncomp=ncomp_chlB_ground,intercept=TRUE))
  car_jack_coefs_ground[[i]]<-as.vector(coef(car_ground_jack,ncomp=ncomp_car_ground,intercept=TRUE))
  Al_jack_coefs_ground[[i]]<-as.vector(coef(Al_ground_jack,ncomp=ncomp_Al_ground,intercept=TRUE))
  Ca_jack_coefs_ground[[i]]<-as.vector(coef(Ca_ground_jack,ncomp=ncomp_Ca_ground,intercept=TRUE))
  Cu_jack_coefs_ground[[i]]<-as.vector(coef(Cu_ground_jack,ncomp=ncomp_Cu_ground,intercept=TRUE))
  Fe_jack_coefs_ground[[i]]<-as.vector(coef(Fe_ground_jack,ncomp=ncomp_Fe_ground,intercept=TRUE))
  K_jack_coefs_ground[[i]]<-as.vector(coef(K_ground_jack,ncomp=ncomp_K_ground,intercept=TRUE))
  Mg_jack_coefs_ground[[i]]<-as.vector(coef(Mg_ground_jack,ncomp=ncomp_Mg_ground,intercept=TRUE))
  Mn_jack_coefs_ground[[i]]<-as.vector(coef(Mn_ground_jack,ncomp=ncomp_Mn_ground,intercept=TRUE))
  Na_jack_coefs_ground[[i]]<-as.vector(coef(Na_ground_jack,ncomp=ncomp_Na_ground,intercept=TRUE))
  P_jack_coefs_ground[[i]]<-as.vector(coef(P_ground_jack,ncomp=ncomp_P_ground,intercept=TRUE))
  Zn_jack_coefs_ground[[i]]<-as.vector(coef(Zn_ground_jack,ncomp=ncomp_Zn_ground,intercept=TRUE))
}

solubles_jack_pred_ground<-apply.coefs(solubles_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
solubles_jack_stat_ground<-t(apply(solubles_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_jack_df_ground<-data.frame(pred_mean=solubles_jack_stat_ground[,1],
                                    pred_low=solubles_jack_stat_ground[,2],
                                    pred_high=solubles_jack_stat_ground[,3],
                                    Measured=meta(ground_spec_EL_agg_test)$solubles,
                                    ncomp=ncomp_solubles_ground,
                                    Project=meta(ground_spec_EL_agg_test)$Project,
                                    GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                    Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                    ID=meta(ground_spec_EL_agg_test)$ID)

hemicellulose_jack_pred_ground<-apply.coefs(hemicellulose_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
hemicellulose_jack_stat_ground<-t(apply(hemicellulose_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_jack_df_ground<-data.frame(pred_mean=hemicellulose_jack_stat_ground[,1],
                                         pred_low=hemicellulose_jack_stat_ground[,2],
                                         pred_high=hemicellulose_jack_stat_ground[,3],
                                         Measured=meta(ground_spec_EL_agg_test)$hemicellulose,
                                         ncomp=ncomp_hemicellulose_ground,
                                         Project=meta(ground_spec_EL_agg_test)$Project,
                                         GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                         Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                         ID=meta(ground_spec_EL_agg_test)$ID)

cellulose_jack_pred_ground<-apply.coefs(cellulose_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
cellulose_jack_stat_ground<-t(apply(cellulose_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_jack_df_ground<-data.frame(pred_mean=cellulose_jack_stat_ground[,1],
                                     pred_low=cellulose_jack_stat_ground[,2],
                                     pred_high=cellulose_jack_stat_ground[,3],
                                     Measured=meta(ground_spec_EL_agg_test)$cellulose,
                                     ncomp=ncomp_cellulose_ground,
                                     Project=meta(ground_spec_EL_agg_test)$Project,
                                     GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                     Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                     ID=meta(ground_spec_EL_agg_test)$ID)

lignin_jack_pred_ground<-apply.coefs(lignin_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
lignin_jack_stat_ground<-t(apply(lignin_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_jack_df_ground<-data.frame(pred_mean=lignin_jack_stat_ground[,1],
                                  pred_low=lignin_jack_stat_ground[,2],
                                  pred_high=lignin_jack_stat_ground[,3],
                                  Measured=meta(ground_spec_EL_agg_test)$lignin,
                                  ncomp=ncomp_lignin_ground,
                                  Project=meta(ground_spec_EL_agg_test)$Project,
                                  GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                  Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                  ID=meta(ground_spec_EL_agg_test)$ID)

perC_jack_pred_ground<-apply.coefs(perC_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
perC_jack_stat_ground<-t(apply(perC_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_jack_df_ground<-data.frame(pred_mean=perC_jack_stat_ground[,1],
                                pred_low=perC_jack_stat_ground[,2],
                                pred_high=perC_jack_stat_ground[,3],
                                Measured=meta(ground_spec_EL_agg_test)$C,
                                ncomp=ncomp_perC_ground,
                                Project=meta(ground_spec_EL_agg_test)$Project,
                                GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                ID=meta(ground_spec_EL_agg_test)$ID)

perN_jack_pred_ground<-apply.coefs(perN_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
perN_jack_stat_ground<-t(apply(perN_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_jack_df_ground<-data.frame(pred_mean=perN_jack_stat_ground[,1],
                                pred_low=perN_jack_stat_ground[,2],
                                pred_high=perN_jack_stat_ground[,3],
                                Measured=meta(ground_spec_EL_agg_test)$N,
                                ncomp=ncomp_perN_ground,
                                Project=meta(ground_spec_EL_agg_test)$Project,
                                GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                ID=meta(ground_spec_EL_agg_test)$ID)

LMA_jack_pred_ground<-apply.coefs(LMA_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
LMA_jack_stat_ground<-t(apply(LMA_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df_ground<-data.frame(pred_mean=LMA_jack_stat_ground[,1],
                               pred_low=LMA_jack_stat_ground[,2],
                               pred_high=LMA_jack_stat_ground[,3],
                               Measured=meta(ground_spec_EL_agg_test)$LMA,
                               ncomp=ncomp_LMA_ground,
                               Project=meta(ground_spec_EL_agg_test)$Project,
                               GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                               ID=meta(ground_spec_EL_agg_test)$ID)

LDMC_jack_pred_ground<-apply.coefs(LDMC_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
LDMC_jack_stat_ground<-t(apply(LDMC_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC_jack_df_ground<-data.frame(pred_mean=LDMC_jack_stat_ground[,1],
                                pred_low=LDMC_jack_stat_ground[,2],
                                pred_high=LDMC_jack_stat_ground[,3],
                                Measured=meta(ground_spec_EL_agg_test)$LDMC,
                                ncomp=ncomp_LDMC_ground,
                                Project=meta(ground_spec_EL_agg_test)$Project,
                                GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                ID=meta(ground_spec_EL_agg_test)$ID)

EWT_jack_pred_ground<-apply.coefs(EWT_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
EWT_jack_stat_ground<-t(apply(EWT_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_jack_df_ground<-data.frame(pred_mean=EWT_jack_stat_ground[,1],
                               pred_low=EWT_jack_stat_ground[,2],
                               pred_high=EWT_jack_stat_ground[,3],
                               Measured=meta(ground_spec_EL_agg_test)$EWT,
                               ncomp=ncomp_EWT_ground,
                               Project=meta(ground_spec_EL_agg_test)$Project,
                               GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                               ID=meta(ground_spec_EL_agg_test)$ID)

chlA_jack_pred_ground<-apply.coefs(chlA_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
chlA_jack_stat_ground<-t(apply(chlA_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_jack_df_ground<-data.frame(pred_mean=chlA_jack_stat_ground[,1],
                                pred_low=chlA_jack_stat_ground[,2],
                                pred_high=chlA_jack_stat_ground[,3],
                                Measured=meta(ground_spec_EL_agg_test)$chlA,
                                ncomp=ncomp_chlA_ground,
                                Project=meta(ground_spec_EL_agg_test)$Project,
                                GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                ID=meta(ground_spec_EL_agg_test)$ID)

chlB_jack_pred_ground<-apply.coefs(chlB_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
chlB_jack_stat_ground<-t(apply(chlB_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_jack_df_ground<-data.frame(pred_mean=chlB_jack_stat_ground[,1],
                                pred_low=chlB_jack_stat_ground[,2],
                                pred_high=chlB_jack_stat_ground[,3],
                                Measured=meta(ground_spec_EL_agg_test)$chlB,
                                ncomp=ncomp_chlB_ground,
                                Project=meta(ground_spec_EL_agg_test)$Project,
                                GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                                Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                                ID=meta(ground_spec_EL_agg_test)$ID)

car_jack_pred_ground<-apply.coefs(car_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
car_jack_stat_ground<-t(apply(car_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_jack_df_ground<-data.frame(pred_mean=car_jack_stat_ground[,1],
                               pred_low=car_jack_stat_ground[,2],
                               pred_high=car_jack_stat_ground[,3],
                               Measured=meta(ground_spec_EL_agg_test)$car,
                               ncomp=ncomp_car_ground,
                               Project=meta(ground_spec_EL_agg_test)$Project,
                               GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                               ID=meta(ground_spec_EL_agg_test)$ID)

Al_jack_pred_ground<-apply.coefs(Al_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Al_jack_stat_ground<-t(apply(Al_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_jack_df_ground<-data.frame(pred_mean=Al_jack_stat_ground[,1],
                              pred_low=Al_jack_stat_ground[,2],
                              pred_high=Al_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Al,
                              ncomp=ncomp_Al_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

Ca_jack_pred_ground<-apply.coefs(Ca_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Ca_jack_stat_ground<-t(apply(Ca_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_jack_df_ground<-data.frame(pred_mean=Ca_jack_stat_ground[,1],
                              pred_low=Ca_jack_stat_ground[,2],
                              pred_high=Ca_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Ca,
                              ncomp=ncomp_Ca_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

Cu_jack_pred_ground<-apply.coefs(Cu_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Cu_jack_stat_ground<-t(apply(Cu_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_jack_df_ground<-data.frame(pred_mean=Cu_jack_stat_ground[,1],
                              pred_low=Cu_jack_stat_ground[,2],
                              pred_high=Cu_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Cu,
                              ncomp=ncomp_Cu_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

Fe_jack_pred_ground<-apply.coefs(Fe_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Fe_jack_stat_ground<-t(apply(Fe_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_jack_df_ground<-data.frame(pred_mean=Fe_jack_stat_ground[,1],
                              pred_low=Fe_jack_stat_ground[,2],
                              pred_high=Fe_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Fe,
                              ncomp=ncomp_Fe_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

K_jack_pred_ground<-apply.coefs(K_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
K_jack_stat_ground<-t(apply(K_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_jack_df_ground<-data.frame(pred_mean=K_jack_stat_ground[,1],
                             pred_low=K_jack_stat_ground[,2],
                             pred_high=K_jack_stat_ground[,3],
                             Measured=meta(ground_spec_EL_agg_test)$K,
                             ncomp=ncomp_K_ground,
                             Project=meta(ground_spec_EL_agg_test)$Project,
                             GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                             ID=meta(ground_spec_EL_agg_test)$ID)

Mg_jack_pred_ground<-apply.coefs(Mg_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Mg_jack_stat_ground<-t(apply(Mg_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_jack_df_ground<-data.frame(pred_mean=Mg_jack_stat_ground[,1],
                              pred_low=Mg_jack_stat_ground[,2],
                              pred_high=Mg_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Mg,
                              ncomp=ncomp_Mg_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

Mn_jack_pred_ground<-apply.coefs(Mn_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Mn_jack_stat_ground<-t(apply(Mn_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_jack_df_ground<-data.frame(pred_mean=Mn_jack_stat_ground[,1],
                              pred_low=Mn_jack_stat_ground[,2],
                              pred_high=Mn_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Mn,
                              ncomp=ncomp_Mn_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

Na_jack_pred_ground<-apply.coefs(Na_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Na_jack_stat_ground<-t(apply(Na_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_jack_df_ground<-data.frame(pred_mean=Na_jack_stat_ground[,1],
                              pred_low=Na_jack_stat_ground[,2],
                              pred_high=Na_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Na,
                              ncomp=ncomp_Na_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

P_jack_pred_ground<-apply.coefs(P_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
P_jack_stat_ground<-t(apply(P_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_jack_df_ground<-data.frame(pred_mean=P_jack_stat_ground[,1],
                             pred_low=P_jack_stat_ground[,2],
                             pred_high=P_jack_stat_ground[,3],
                             Measured=meta(ground_spec_EL_agg_test)$P,
                             ncomp=ncomp_P_ground,
                             Project=meta(ground_spec_EL_agg_test)$Project,
                             GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                             ID=meta(ground_spec_EL_agg_test)$ID)

Zn_jack_pred_ground<-apply.coefs(Zn_jack_coefs_ground,as.matrix(ground_spec_EL_agg_test))
Zn_jack_stat_ground<-t(apply(Zn_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_jack_df_ground<-data.frame(pred_mean=Zn_jack_stat_ground[,1],
                              pred_low=Zn_jack_stat_ground[,2],
                              pred_high=Zn_jack_stat_ground[,3],
                              Measured=meta(ground_spec_EL_agg_test)$Zn,
                              ncomp=ncomp_Zn_ground,
                              Project=meta(ground_spec_EL_agg_test)$Project,
                              GrowthForm=meta(ground_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(ground_spec_EL_agg_test)$Discoloration,
                              ID=meta(ground_spec_EL_agg_test)$ID)

##########################################################
## save output

ground_jack_coef_list<-list(LMA=LMA_jack_coefs_ground,
                             LDMC=LDMC_jack_coefs_ground,
                             EWT=EWT_jack_coefs_ground,
                             C=perC_jack_coefs_ground,
                             N=perN_jack_coefs_ground,
                             sol=solubles_jack_coefs_ground,
                             hemi=hemicellulose_jack_coefs_ground,
                             cell=cellulose_jack_coefs_ground,
                             lign=lignin_jack_coefs_ground,
                             chlA=chlA_jack_coefs_ground,
                             chlB=chlB_jack_coefs_ground,
                             car=car_jack_coefs_ground,
                             Al=Al_jack_coefs_ground,
                             Ca=Ca_jack_coefs_ground,
                             Cu=Cu_jack_coefs_ground,
                             Fe=Fe_jack_coefs_ground,
                             K=K_jack_coefs_ground,
                             Mg=Mg_jack_coefs_ground,
                             Mn=Mn_jack_coefs_ground,
                             Na=Na_jack_coefs_ground,
                             P=P_jack_coefs_ground,
                             Zn=Zn_jack_coefs_ground)
saveRDS(ground_jack_coef_list,"SavedResults/ground_jack_coefs_list.rds")

ground_jack_df_list<-list(LMA=LMA_jack_df_ground,
                           LDMC=LDMC_jack_df_ground,
                           EWT=EWT_jack_df_ground,
                           C=perC_jack_df_ground,
                           N=perN_jack_df_ground,
                           sol=solubles_jack_df_ground,
                           hemi=hemicellulose_jack_df_ground,
                           cell=cellulose_jack_df_ground,
                           lign=lignin_jack_df_ground,
                           chlA=chlA_jack_df_ground,
                           chlB=chlB_jack_df_ground,
                           car=car_jack_df_ground,
                           Al=Al_jack_df_ground,
                           Ca=Ca_jack_df_ground,
                           Cu=Cu_jack_df_ground,
                           Fe=Fe_jack_df_ground,
                           K=K_jack_df_ground,
                           Mg=Mg_jack_df_ground,
                           Mn=Mn_jack_df_ground,
                           Na=Na_jack_df_ground,
                           P=P_jack_df_ground,
                           Zn=Zn_jack_df_ground)
saveRDS(ground_jack_df_list,"SavedResults/ground_jack_df_list.rds")

######################################
## violin plots

ground_val_summary<-data.frame(variable=names(ground_jack_df_list),
                                perRMSE=unlist(lapply(ground_jack_df_list,
                                                      function(x) percentRMSD(x$Measured,x$pred_mean,0.025,0.975))),
                                R2=unlist(lapply(ground_jack_df_list,
                                                 function(x) summary(lm(Measured~pred_mean,data=x))$r.squared)))

R2.df_ground<-data.frame(LMA=unlist(lapply(LMA_jack_stats_ground,function(x) x[["R2"]])),
                         LDMC=unlist(lapply(LDMC_jack_stats_ground,function(x) x[["R2"]])),
                         EWT=unlist(lapply(EWT_jack_stats_ground,function(x) x[["R2"]])),
                         sol=unlist(lapply(solubles_jack_stats_ground,function(x) x[["R2"]])),
                        hemi=unlist(lapply(hemicellulose_jack_stats_ground,function(x) x[["R2"]])),
                        cell=unlist(lapply(cellulose_jack_stats_ground,function(x) x[["R2"]])),
                        lign=unlist(lapply(lignin_jack_stats_ground,function(x) x[["R2"]])),
                        C=unlist(lapply(perC_jack_stats_ground,function(x) x[["R2"]])),
                        N=unlist(lapply(perN_jack_stats_ground,function(x) x[["R2"]])),
                        LMA=unlist(lapply(LMA_jack_stats_ground,function(x) x[["R2"]])),
                        LDMC=unlist(lapply(LDMC_jack_stats_ground,function(x) x[["R2"]])),
                        EWT=unlist(lapply(EWT_jack_stats_ground,function(x) x[["R2"]])),
                        chlA=unlist(lapply(chlA_jack_stats_ground,function(x) x[["R2"]])),
                        chlB=unlist(lapply(chlB_jack_stats_ground,function(x) x[["R2"]])),
                        car=unlist(lapply(car_jack_stats_ground,function(x) x[["R2"]])),
                        Al=unlist(lapply(Al_jack_stats_ground,function(x) x[["R2"]])),
                        Ca=unlist(lapply(Ca_jack_stats_ground,function(x) x[["R2"]])),
                        Cu=unlist(lapply(Cu_jack_stats_ground,function(x) x[["R2"]])),
                        Fe=unlist(lapply(Fe_jack_stats_ground,function(x) x[["R2"]])),
                        K=unlist(lapply(K_jack_stats_ground,function(x) x[["R2"]])),
                        Mg=unlist(lapply(Mg_jack_stats_ground,function(x) x[["R2"]])),
                        Mn=unlist(lapply(Mn_jack_stats_ground,function(x) x[["R2"]])),
                        Na=unlist(lapply(Na_jack_stats_ground,function(x) x[["R2"]])),
                        P=unlist(lapply(P_jack_stats_ground,function(x) x[["R2"]])),
                        Zn=unlist(lapply(Zn_jack_stats_ground,function(x) x[["R2"]])))

R2.long_ground<-melt(R2.df_ground)
ground_val_R2<-ggplot(R2.long_ground,aes(y=value,x=variable))+
  geom_violin()+
  geom_point(data=ground_val_summary,
             aes(x=variable,y=R2),color="red",size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y=expression(italic("R"^2)))+
  ggtitle("Ground-leaf spectra")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))

perRMSE.df_ground<-data.frame(LMA=unlist(lapply(LMA_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                              LDMC=unlist(lapply(LDMC_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                              EWT=unlist(lapply(EWT_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                              sol=unlist(lapply(solubles_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             hemi=unlist(lapply(hemicellulose_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             cell=unlist(lapply(cellulose_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             lign=unlist(lapply(lignin_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             C=unlist(lapply(perC_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             N=unlist(lapply(perN_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             chlA=unlist(lapply(chlA_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             chlB=unlist(lapply(chlB_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             car=unlist(lapply(car_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Al=unlist(lapply(Al_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Ca=unlist(lapply(Ca_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Cu=unlist(lapply(Cu_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Fe=unlist(lapply(Fe_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             K=unlist(lapply(K_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Mg=unlist(lapply(Mg_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Mn=unlist(lapply(Mn_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Na=unlist(lapply(Na_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             P=unlist(lapply(P_jack_stats_ground,function(x) 100*x[["perRMSE"]])),
                             Zn=unlist(lapply(Zn_jack_stats_ground,function(x) 100*x[["perRMSE"]])))

perRMSE.long_ground<-melt(perRMSE.df_ground)
ground_val_perRMSE<-ggplot(perRMSE.long_ground,aes(y=value,x=variable))+
  geom_violin()+
  geom_point(data=ground_val_summary,
             aes(x=variable,y=perRMSE*100),color="red",size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  labs(y="%RMSE")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,max(perRMSE.long_ground$value)*1.1))

pdf("Manuscript/FigS13.pdf",width=8,height=6)
egg::ggarrange(plots = list(ground_val_R2,ground_val_perRMSE),
               nrow=2,ncol=1)
dev.off()
