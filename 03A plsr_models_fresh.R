setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(pls)
library(ggplot2)
library(caret)
library(reshape2)
library(RColorBrewer)
library(patchwork)

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

fresh_spec_EL_agg_train<-readRDS("ProcessedSpectralData/fresh_spec_EL_agg_train.rds")
fresh_spec_EL_agg_test<-readRDS("ProcessedSpectralData/fresh_spec_EL_agg_test.rds")

################################################
## building calibration models

perC_fresh<-plsr(meta(fresh_spec_EL_agg_train)$C~as.matrix(fresh_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perC_fresh <- selectNcomp(perC_fresh, method = "onesigma", plot = FALSE)
perC_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$C))
perC_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[perC_fresh_valid],
                            Species=meta(fresh_spec_EL_agg_train)$Species[perC_fresh_valid],
                            Project=meta(fresh_spec_EL_agg_train)$Project[perC_fresh_valid],
                            Stage=meta(fresh_spec_EL_agg_train)$Stage[perC_fresh_valid],
                            GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[perC_fresh_valid],
                            measured=meta(fresh_spec_EL_agg_train)$C[perC_fresh_valid],
                            val_pred=perC_fresh$validation$pred[,,ncomp_perC_fresh])
ggplot(perC_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %C from fresh-leaf spectra")

perN_fresh<-plsr(meta(fresh_spec_EL_agg_train)$N~as.matrix(fresh_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perN_fresh <- selectNcomp(perN_fresh, method = "onesigma", plot = FALSE)
perN_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$N))
perN_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[perN_fresh_valid],
                            Species=meta(fresh_spec_EL_agg_train)$Species[perN_fresh_valid],
                            Project=meta(fresh_spec_EL_agg_train)$Project[perN_fresh_valid],
                            Stage=meta(fresh_spec_EL_agg_train)$Stage[perN_fresh_valid],
                            GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[perC_fresh_valid],
                            measured=meta(fresh_spec_EL_agg_train)$N[perN_fresh_valid],
                            val_pred=perN_fresh$validation$pred[,,ncomp_perN_fresh])
ggplot(perN_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %N from fresh-leaf spectra")+guides(color=F)

LMA_fresh<-plsr(meta(fresh_spec_EL_agg_train)$LMA~as.matrix(fresh_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_fresh <- selectNcomp(LMA_fresh, method = "onesigma", plot = FALSE)
LMA_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$LMA))
LMA_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[LMA_fresh_valid],
                           Species=meta(fresh_spec_EL_agg_train)$Species[LMA_fresh_valid],
                           Project=meta(fresh_spec_EL_agg_train)$Project[LMA_fresh_valid],
                           Stage=meta(fresh_spec_EL_agg_train)$Stage[LMA_fresh_valid],
                           GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[LMA_fresh_valid],
                           measured=meta(fresh_spec_EL_agg_train)$LMA[LMA_fresh_valid],
                           val_pred=LMA_fresh$validation$pred[,,ncomp_LMA_fresh])
ggplot(LMA_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.25))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting LMA from fresh-leaf spectra")

LDMC_fresh<-plsr(meta(fresh_spec_EL_agg_train)$LDMC~as.matrix(fresh_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_fresh <- selectNcomp(LDMC_fresh, method = "onesigma", plot = FALSE)
LDMC_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$LDMC))
LDMC_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[LDMC_fresh_valid],
                            Species=meta(fresh_spec_EL_agg_train)$Species[LDMC_fresh_valid],
                            Project=meta(fresh_spec_EL_agg_train)$Project[LDMC_fresh_valid],
                            Stage=meta(fresh_spec_EL_agg_train)$Stage[LDMC_fresh_valid],
                            GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[LDMC_fresh_valid],
                            measured=meta(fresh_spec_EL_agg_train)$LDMC[LDMC_fresh_valid],
                            val_pred=LDMC_fresh$validation$pred[,,ncomp_LDMC_fresh])
ggplot(LDMC_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(150,600),ylim=c(150,600))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting LDMC from fresh-leaf spectra")

EWT_fresh<-plsr(meta(fresh_spec_EL_agg_train)$EWT~as.matrix(fresh_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_fresh <- selectNcomp(EWT_fresh, method = "onesigma", plot = FALSE)
EWT_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$EWT))
EWT_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[EWT_fresh_valid],
                           Species=meta(fresh_spec_EL_agg_train)$Species[EWT_fresh_valid],
                           Project=meta(fresh_spec_EL_agg_train)$Project[EWT_fresh_valid],
                           Stage=meta(fresh_spec_EL_agg_train)$Stage[EWT_fresh_valid],
                           GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[EWT_fresh_valid],
                           measured=meta(fresh_spec_EL_agg_train)$EWT[EWT_fresh_valid],
                           val_pred=EWT_fresh$validation$pred[,,ncomp_EWT_fresh])
ggplot(EWT_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.03),ylim=c(0,0.03))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting EWT from fresh-leaf spectra")

chlA_fresh<-plsr(meta(fresh_spec_EL_agg_train)$chlA~as.matrix(fresh_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_fresh <- selectNcomp(chlA_fresh, method = "onesigma", plot = FALSE)
chlA_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$chlA))
chlA_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[chlA_fresh_valid],
                            Species=meta(fresh_spec_EL_agg_train)$Species[chlA_fresh_valid],
                            Project=meta(fresh_spec_EL_agg_train)$Project[chlA_fresh_valid],
                            Stage=meta(fresh_spec_EL_agg_train)$Stage[chlA_fresh_valid],
                            GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[chlA_fresh_valid],
                            measured=meta(fresh_spec_EL_agg_train)$chlA[chlA_fresh_valid],
                            val_pred=chlA_fresh$validation$pred[,,ncomp_chlA_fresh])
ggplot(chlA_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,16),ylim=c(0,16))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl a from fresh-leaf spectra")

chlB_fresh<-plsr(meta(fresh_spec_EL_agg_train)$chlB~as.matrix(fresh_spec_EL_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_fresh <- selectNcomp(chlB_fresh, method = "onesigma", plot = FALSE)
chlB_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$chlB))
chlB_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[chlB_fresh_valid],
                            Species=meta(fresh_spec_EL_agg_train)$Species[chlB_fresh_valid],
                            Project=meta(fresh_spec_EL_agg_train)$Project[chlB_fresh_valid],
                            Stage=meta(fresh_spec_EL_agg_train)$Stage[chlB_fresh_valid],
                            GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[chlB_fresh_valid],
                            measured=meta(fresh_spec_EL_agg_train)$chlB[chlB_fresh_valid],
                            val_pred=chlB_fresh$validation$pred[,,ncomp_chlB_fresh])
ggplot(chlB_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl b from fresh-leaf spectra")

car_fresh<-plsr(meta(fresh_spec_EL_agg_train)$car~as.matrix(fresh_spec_EL_agg_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_fresh <- selectNcomp(car_fresh, method = "onesigma", plot = FALSE)
car_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$car))
car_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[car_fresh_valid],
                           Species=meta(fresh_spec_EL_agg_train)$Species[car_fresh_valid],
                           Project=meta(fresh_spec_EL_agg_train)$Project[car_fresh_valid],
                           Stage=meta(fresh_spec_EL_agg_train)$Stage[car_fresh_valid],
                           GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[car_fresh_valid],
                           measured=meta(fresh_spec_EL_agg_train)$car[car_fresh_valid],
                           val_pred=car_fresh$validation$pred[,,ncomp_car_fresh])
ggplot(car_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,3.5),ylim=c(0,3.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting total car from fresh-leaf spectra")

solubles_fresh<-plsr(meta(fresh_spec_EL_agg_train)$solubles~as.matrix(fresh_spec_EL_agg_train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_fresh <- selectNcomp(solubles_fresh, method = "onesigma", plot = FALSE)
solubles_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$solubles))
solubles_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[solubles_fresh_valid],
                                Species=meta(fresh_spec_EL_agg_train)$Species[solubles_fresh_valid],
                                Project=meta(fresh_spec_EL_agg_train)$Project[solubles_fresh_valid],
                                Stage=meta(fresh_spec_EL_agg_train)$Stage[solubles_fresh_valid],
                                GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[solubles_fresh_valid],
                                measured=meta(fresh_spec_EL_agg_train)$solubles[solubles_fresh_valid],
                                val_pred=solubles_fresh$validation$pred[,,ncomp_solubles_fresh])
ggplot(solubles_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(35,90),ylim=c(35,90))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting solubles from fresh-leaf spectra")

hemicellulose_fresh<-plsr(meta(fresh_spec_EL_agg_train)$hemicellulose~as.matrix(fresh_spec_EL_agg_train),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_fresh <- selectNcomp(hemicellulose_fresh, method = "onesigma", plot = FALSE)
hemicellulose_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$hemicellulose))
hemicellulose_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[hemicellulose_fresh_valid],
                                     Species=meta(fresh_spec_EL_agg_train)$Species[hemicellulose_fresh_valid],
                                     Project=meta(fresh_spec_EL_agg_train)$Project[hemicellulose_fresh_valid],
                                     Stage=meta(fresh_spec_EL_agg_train)$Stage[hemicellulose_fresh_valid],
                                     GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[hemicellulose_fresh_valid],
                                     measured=meta(fresh_spec_EL_agg_train)$hemicellulose[hemicellulose_fresh_valid],
                                     val_pred=hemicellulose_fresh$validation$pred[,,ncomp_hemicellulose_fresh])
ggplot(hemicellulose_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,36),ylim=c(0,36))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting hemicellulose from fresh-leaf spectra")

cellulose_fresh<-plsr(meta(fresh_spec_EL_agg_train)$cellulose~as.matrix(fresh_spec_EL_agg_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_fresh <- selectNcomp(cellulose_fresh, method = "onesigma", plot = FALSE)
cellulose_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$cellulose))
cellulose_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[cellulose_fresh_valid],
                                 Species=meta(fresh_spec_EL_agg_train)$Species[cellulose_fresh_valid],
                                 Project=meta(fresh_spec_EL_agg_train)$Project[cellulose_fresh_valid],
                                 Stage=meta(fresh_spec_EL_agg_train)$Stage[cellulose_fresh_valid],
                                 GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[cellulose_fresh_valid],
                                 measured=meta(fresh_spec_EL_agg_train)$cellulose[cellulose_fresh_valid],
                                 val_pred=cellulose_fresh$validation$pred[,,ncomp_cellulose_fresh])
ggplot(cellulose_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(2,28),ylim=c(2,28))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting cellulose from fresh-leaf spectra")

lignin_fresh<-plsr(meta(fresh_spec_EL_agg_train)$lignin~as.matrix(fresh_spec_EL_agg_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_fresh <- selectNcomp(lignin_fresh, method = "onesigma", plot = FALSE)
lignin_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$lignin))
lignin_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[lignin_fresh_valid],
                              Species=meta(fresh_spec_EL_agg_train)$Species[lignin_fresh_valid],
                              Project=meta(fresh_spec_EL_agg_train)$Project[lignin_fresh_valid],
                              Stage=meta(fresh_spec_EL_agg_train)$Stage[lignin_fresh_valid],
                              GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[lignin_fresh_valid],
                              measured=meta(fresh_spec_EL_agg_train)$lignin[lignin_fresh_valid],
                              val_pred=lignin_fresh$validation$pred[,,ncomp_lignin_fresh])
ggplot(lignin_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,22),ylim=c(0,22))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting lignin from fresh-leaf spectra")

Al_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Al~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_fresh <- selectNcomp(Al_fresh, method = "onesigma", plot = FALSE)
Al_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Al))
Al_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Al_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Al_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Al_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Al_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Al_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Al[Al_fresh_valid],
                          val_pred=Al_fresh$validation$pred[,,ncomp_Al_fresh])
ggplot(Al_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.05,0.35),ylim=c(-0.05,0.35))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Al from fresh-leaf spectra")

Ca_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Ca~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_fresh <- selectNcomp(Ca_fresh, method = "onesigma", plot = FALSE)
Ca_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Ca))
Ca_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Ca_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Ca_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Ca_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Ca_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Ca_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Ca[Ca_fresh_valid],
                          val_pred=Ca_fresh$validation$pred[,,ncomp_Ca_fresh])
ggplot(Ca_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-10,40),ylim=c(-10,40))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Ca from fresh-leaf spectra")

Cu_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Cu~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_fresh <- selectNcomp(Cu_fresh, method = "onesigma", plot = FALSE)
Cu_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Cu))
Cu_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Cu_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Cu_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Cu_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Cu_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Cu_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Cu[Cu_fresh_valid],
                          val_pred=Cu_fresh$validation$pred[,,ncomp_Cu_fresh])
ggplot(Cu_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.005,0.055),ylim=c(-0.005,0.055))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Cu from fresh-leaf spectra")

Fe_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Fe~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_fresh <- selectNcomp(Fe_fresh, method = "onesigma", plot = FALSE)
Fe_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Fe))
Fe_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Fe_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Fe_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Fe_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Fe_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Fe_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Fe[Fe_fresh_valid],
                          val_pred=Fe_fresh$validation$pred[,,ncomp_Fe_fresh])
ggplot(Fe_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Fe from fresh-leaf spectra")

K_fresh<-plsr(meta(fresh_spec_EL_agg_train)$K~as.matrix(fresh_spec_EL_agg_train),
              ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_fresh <- selectNcomp(K_fresh, method = "onesigma", plot = FALSE)
K_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$K))
K_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[K_fresh_valid],
                         Species=meta(fresh_spec_EL_agg_train)$Species[K_fresh_valid],
                         Project=meta(fresh_spec_EL_agg_train)$Project[K_fresh_valid],
                         Stage=meta(fresh_spec_EL_agg_train)$Stage[K_fresh_valid],
                         GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[K_fresh_valid],
                         measured=meta(fresh_spec_EL_agg_train)$K[K_fresh_valid],
                         val_pred=K_fresh$validation$pred[,,ncomp_K_fresh])
ggplot(K_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,35),ylim=c(0,35))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting K from fresh-leaf spectra")

Mg_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Mg~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_fresh <- selectNcomp(Mg_fresh, method = "onesigma", plot = FALSE)
Mg_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Mg))
Mg_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Mg_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Mg_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Mg_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Mg_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Mg_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Mg[Mg_fresh_valid],
                          val_pred=Mg_fresh$validation$pred[,,ncomp_Mg_fresh])
ggplot(Mg_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,8),ylim=c(0,8))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mg from fresh-leaf spectra")

Mn_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Mn~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_fresh <- selectNcomp(Mn_fresh, method = "onesigma", plot = FALSE)
Mn_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Mn))
Mn_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Mn_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Mn_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Mn_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Mn_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Mn_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Mn[Mn_fresh_valid],
                          val_pred=Mn_fresh$validation$pred[,,ncomp_Mn_fresh])
ggplot(Mn_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,1.1),ylim=c(-0.1,1.1))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mn from fresh-leaf spectra")

Na_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Na~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_fresh <- selectNcomp(Na_fresh, method = "onesigma", plot = FALSE)
Na_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Na))
Na_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Na_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Na_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Na_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Na_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Na_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Na[Na_fresh_valid],
                          val_pred=Na_fresh$validation$pred[,,ncomp_Na_fresh])
ggplot(Na_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,5),ylim=c(-0.5,5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Na from fresh-leaf spectra")

P_fresh<-plsr(meta(fresh_spec_EL_agg_train)$P~as.matrix(fresh_spec_EL_agg_train),
              ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_fresh <- selectNcomp(P_fresh, method = "onesigma", plot = FALSE)
P_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$P))
P_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[P_fresh_valid],
                         Species=meta(fresh_spec_EL_agg_train)$Species[P_fresh_valid],
                         Project=meta(fresh_spec_EL_agg_train)$Project[P_fresh_valid],
                         Stage=meta(fresh_spec_EL_agg_train)$Stage[P_fresh_valid],
                         GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[P_fresh_valid],
                         measured=meta(fresh_spec_EL_agg_train)$P[P_fresh_valid],
                         val_pred=P_fresh$validation$pred[,,ncomp_P_fresh])
ggplot(P_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,8),ylim=c(0,8))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting P from fresh-leaf spectra")

Zn_fresh<-plsr(meta(fresh_spec_EL_agg_train)$Zn~as.matrix(fresh_spec_EL_agg_train),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_fresh <- selectNcomp(Zn_fresh, method = "onesigma", plot = FALSE)
Zn_fresh_valid <- which(!is.na(meta(fresh_spec_EL_agg_train)$Zn))
Zn_fresh_pred<-data.frame(ID=meta(fresh_spec_EL_agg_train)$ID[Zn_fresh_valid],
                          Species=meta(fresh_spec_EL_agg_train)$Species[Zn_fresh_valid],
                          Project=meta(fresh_spec_EL_agg_train)$Project[Zn_fresh_valid],
                          Stage=meta(fresh_spec_EL_agg_train)$Stage[Zn_fresh_valid],
                          GrowthForm=meta(fresh_spec_EL_agg_train)$GrowthForm[Zn_fresh_valid],
                          measured=meta(fresh_spec_EL_agg_train)$Zn[Zn_fresh_valid],
                          val_pred=Zn_fresh$validation$pred[,,ncomp_Zn_fresh])
ggplot(Zn_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,0.7),ylim=c(-0.1,0.7))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Zn from fresh-leaf spectra")

###############################################
## VIP plots

source("VIP.R")

VIP_fresh<-data.frame(LMA=VIP(LMA_fresh)[ncomp_LMA_fresh,],
                      LDMC=VIP(LDMC_fresh)[ncomp_LDMC_fresh,],
                      EWT=VIP(EWT_fresh)[ncomp_EWT_fresh,],
                      sol=VIP(solubles_fresh)[ncomp_solubles_fresh,],
                      hemi=VIP(hemicellulose_fresh)[ncomp_hemicellulose_fresh,],
                      cell=VIP(cellulose_fresh)[ncomp_cellulose_fresh,],
                      lign=VIP(lignin_fresh)[ncomp_lignin_fresh,],
                      chlA=VIP(chlA_fresh)[ncomp_chlA_fresh,],
                      chlB=VIP(chlB_fresh)[ncomp_chlB_fresh,],
                      car=VIP(car_fresh)[ncomp_car_fresh,],
                      C=VIP(perC_fresh)[ncomp_perC_fresh,],
                      N=VIP(perN_fresh)[ncomp_perN_fresh,],
                      Al=VIP(Al_fresh)[ncomp_Al_fresh,],
                      Ca=VIP(Ca_fresh)[ncomp_Ca_fresh,],
                      Cu=VIP(Cu_fresh)[ncomp_Cu_fresh,],
                      Fe=VIP(Fe_fresh)[ncomp_Fe_fresh,],
                      K=VIP(K_fresh)[ncomp_K_fresh,],
                      Mg=VIP(Mg_fresh)[ncomp_Mg_fresh,],
                      Mn=VIP(Mn_fresh)[ncomp_Mn_fresh,],
                      Na=VIP(Na_fresh)[ncomp_Na_fresh,],
                      P=VIP(P_fresh)[ncomp_P_fresh,],
                      Zn=VIP(Zn_fresh)[ncomp_Zn_fresh,],
                      wavelength=400:2400)

saveRDS(VIP_fresh,"SavedResults/VIP_fresh.rds")

###############################################
## jackknife tests + prediction of validation data
## fresh leaves

solubles_jack_coefs_fresh<-list()
hemicellulose_jack_coefs_fresh<-list()
cellulose_jack_coefs_fresh<-list()
lignin_jack_coefs_fresh<-list()
perC_jack_coefs_fresh<-list()
perN_jack_coefs_fresh<-list()
LMA_jack_coefs_fresh<-list()
LDMC_jack_coefs_fresh<-list()
EWT_jack_coefs_fresh<-list()
chlA_jack_coefs_fresh<-list()
chlB_jack_coefs_fresh<-list()
car_jack_coefs_fresh<-list()
Al_jack_coefs_fresh<-list()
Ca_jack_coefs_fresh<-list()
Cu_jack_coefs_fresh<-list()
Fe_jack_coefs_fresh<-list()
K_jack_coefs_fresh<-list()
Mg_jack_coefs_fresh<-list()
Mn_jack_coefs_fresh<-list()
Na_jack_coefs_fresh<-list()
P_jack_coefs_fresh<-list()
Zn_jack_coefs_fresh<-list()

solubles_jack_stats_fresh<-list()
hemicellulose_jack_stats_fresh<-list()
cellulose_jack_stats_fresh<-list()
lignin_jack_stats_fresh<-list()
perC_jack_stats_fresh<-list()
perN_jack_stats_fresh<-list()
LMA_jack_stats_fresh<-list()
LDMC_jack_stats_fresh<-list()
EWT_jack_stats_fresh<-list()
chlA_jack_stats_fresh<-list()
chlB_jack_stats_fresh<-list()
car_jack_stats_fresh<-list()
Al_jack_stats_fresh<-list()
Ca_jack_stats_fresh<-list()
Cu_jack_stats_fresh<-list()
Fe_jack_stats_fresh<-list()
K_jack_stats_fresh<-list()
Mg_jack_stats_fresh<-list()
Mn_jack_stats_fresh<-list()
Na_jack_stats_fresh<-list()
P_jack_stats_fresh<-list()
Zn_jack_stats_fresh<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_fresh<-nrow(fresh_spec_EL_agg_train)
  train_jack_fresh<-sample(1:n_cal_spec_fresh,floor(0.7*n_cal_spec_fresh))
  test_jack_fresh<-setdiff(1:n_cal_spec_fresh,train_jack_fresh)
  
  calib_jack_fresh<-fresh_spec_EL_agg_train[train_jack_fresh]
  val_jack_fresh<-fresh_spec_EL_agg_train[test_jack_fresh]
  
  solubles_fresh_jack<-plsr(meta(calib_jack_fresh)$solubles~as.matrix(calib_jack_fresh),
                            ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_fresh_jack<-plsr(meta(calib_jack_fresh)$hemicellulose~as.matrix(calib_jack_fresh),
                                 ncomp=30,method = "oscorespls",validation="none")
  cellulose_fresh_jack<-plsr(meta(calib_jack_fresh)$cellulose~as.matrix(calib_jack_fresh),
                             ncomp=30,method = "oscorespls",validation="none")
  lignin_fresh_jack<-plsr(meta(calib_jack_fresh)$lignin~as.matrix(calib_jack_fresh),
                          ncomp=30,method = "oscorespls",validation="none")
  perC_fresh_jack<-plsr(meta(calib_jack_fresh)$C~as.matrix(calib_jack_fresh),
                        ncomp=30,method = "oscorespls",validation="none")
  perN_fresh_jack<-plsr(meta(calib_jack_fresh)$N~as.matrix(calib_jack_fresh),
                        ncomp=30,method = "oscorespls",validation="none")
  LMA_fresh_jack<-plsr(meta(calib_jack_fresh)$LMA~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  LDMC_fresh_jack<-plsr(meta(calib_jack_fresh)$LDMC~as.matrix(calib_jack_fresh),
                        ncomp=30,method = "oscorespls",validation="none")
  EWT_fresh_jack<-plsr(meta(calib_jack_fresh)$EWT~as.matrix(calib_jack_fresh),
                        ncomp=30,method = "oscorespls",validation="none")
  chlA_fresh_jack<-plsr(meta(calib_jack_fresh)$chlA~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  chlB_fresh_jack<-plsr(meta(calib_jack_fresh)$chlB~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  car_fresh_jack<-plsr(meta(calib_jack_fresh)$car~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Al_fresh_jack<-plsr(meta(calib_jack_fresh)$Al~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Ca_fresh_jack<-plsr(meta(calib_jack_fresh)$Ca~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Cu_fresh_jack<-plsr(meta(calib_jack_fresh)$Cu~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Fe_fresh_jack<-plsr(meta(calib_jack_fresh)$Fe~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  K_fresh_jack<-plsr(meta(calib_jack_fresh)$K~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Mg_fresh_jack<-plsr(meta(calib_jack_fresh)$Mg~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Mn_fresh_jack<-plsr(meta(calib_jack_fresh)$Mn~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Na_fresh_jack<-plsr(meta(calib_jack_fresh)$Na~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  P_fresh_jack<-plsr(meta(calib_jack_fresh)$P~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  Zn_fresh_jack<-plsr(meta(calib_jack_fresh)$Zn~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")
  
  solubles_jack_val_pred_fresh<-as.vector(predict(solubles_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_solubles_fresh)[,,1])
  solubles_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$solubles~solubles_jack_val_pred_fresh)
  solubles_jack_stats_fresh[[i]]<-c(R2=summary(solubles_jack_val_fit_fresh)$r.squared,
                                    RMSE=RMSD(meta(val_jack_fresh)$solubles,solubles_jack_val_pred_fresh),
                                    perRMSE=percentRMSD(meta(val_jack_fresh)$solubles,solubles_jack_val_pred_fresh,0.025,0.975),
                                    bias=mean(solubles_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$solubles,na.rm=T))
  
  hemicellulose_jack_val_pred_fresh<-as.vector(predict(hemicellulose_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_hemicellulose_fresh)[,,1])
  hemicellulose_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$hemicellulose~hemicellulose_jack_val_pred_fresh)
  hemicellulose_jack_stats_fresh[[i]]<-c(R2=summary(hemicellulose_jack_val_fit_fresh)$r.squared,
                                         RMSE=RMSD(meta(val_jack_fresh)$hemicellulose,hemicellulose_jack_val_pred_fresh),
                                         perRMSE=percentRMSD(meta(val_jack_fresh)$hemicellulose,hemicellulose_jack_val_pred_fresh,0.025,0.975),
                                         bias=mean(hemicellulose_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$hemicellulose,na.rm=T))
  
  cellulose_jack_val_pred_fresh<-as.vector(predict(cellulose_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_cellulose_fresh)[,,1])
  cellulose_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$cellulose~cellulose_jack_val_pred_fresh)
  cellulose_jack_stats_fresh[[i]]<-c(R2=summary(cellulose_jack_val_fit_fresh)$r.squared,
                                     RMSE=RMSD(meta(val_jack_fresh)$cellulose,cellulose_jack_val_pred_fresh),
                                     perRMSE=percentRMSD(meta(val_jack_fresh)$cellulose,cellulose_jack_val_pred_fresh,0.025,0.975),
                                     bias=mean(cellulose_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$cellulose,na.rm=T))
  
  lignin_jack_val_pred_fresh<-as.vector(predict(lignin_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_lignin_fresh)[,,1])
  lignin_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$lignin~lignin_jack_val_pred_fresh)
  lignin_jack_stats_fresh[[i]]<-c(R2=summary(lignin_jack_val_fit_fresh)$r.squared,
                                  RMSE=RMSD(meta(val_jack_fresh)$lignin,lignin_jack_val_pred_fresh),
                                  perRMSE=percentRMSD(meta(val_jack_fresh)$lignin,lignin_jack_val_pred_fresh,0.025,0.975),
                                  bias=mean(lignin_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$lignin,na.rm=T))
  
  perC_jack_val_pred_fresh<-as.vector(predict(perC_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_perC_fresh)[,,1])
  perC_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$C~perC_jack_val_pred_fresh)
  perC_jack_stats_fresh[[i]]<-c(R2=summary(perC_jack_val_fit_fresh)$r.squared,
                                RMSE=RMSD(meta(val_jack_fresh)$C,perC_jack_val_pred_fresh),
                                perRMSE=percentRMSD(meta(val_jack_fresh)$C,perC_jack_val_pred_fresh,0.025,0.975),
                                bias=mean(perC_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$C,na.rm=T))
  
  perN_jack_val_pred_fresh<-as.vector(predict(perN_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_perN_fresh)[,,1])
  perN_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$N~perN_jack_val_pred_fresh)
  perN_jack_stats_fresh[[i]]<-c(R2=summary(perN_jack_val_fit_fresh)$r.squared,
                                RMSE=RMSD(meta(val_jack_fresh)$N,perN_jack_val_pred_fresh),
                                perRMSE=percentRMSD(meta(val_jack_fresh)$N,perN_jack_val_pred_fresh,0.025,0.975),
                                bias=mean(perN_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$N,na.rm=T))
  
  LMA_jack_val_pred_fresh<-as.vector(predict(LMA_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_LMA_fresh)[,,1])
  LMA_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$LMA~LMA_jack_val_pred_fresh)
  LMA_jack_stats_fresh[[i]]<-c(R2=summary(LMA_jack_val_fit_fresh)$r.squared,
                               RMSE=RMSD(meta(val_jack_fresh)$LMA,LMA_jack_val_pred_fresh),
                               perRMSE=percentRMSD(meta(val_jack_fresh)$LMA,LMA_jack_val_pred_fresh,0.025,0.975),
                               bias=mean(LMA_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$LMA,na.rm=T))
  
  LDMC_jack_val_pred_fresh<-as.vector(predict(LDMC_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_LDMC_fresh)[,,1])
  LDMC_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$LDMC~LDMC_jack_val_pred_fresh)
  LDMC_jack_stats_fresh[[i]]<-c(R2=summary(LDMC_jack_val_fit_fresh)$r.squared,
                                RMSE=RMSD(meta(val_jack_fresh)$LDMC,LDMC_jack_val_pred_fresh),
                                perRMSE=percentRMSD(meta(val_jack_fresh)$LDMC,LDMC_jack_val_pred_fresh,0.025,0.975),
                                bias=mean(LDMC_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$LDMC,na.rm=T))

  EWT_jack_val_pred_fresh<-as.vector(predict(EWT_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_EWT_fresh)[,,1])
  EWT_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$EWT~EWT_jack_val_pred_fresh)
  EWT_jack_stats_fresh[[i]]<-c(R2=summary(EWT_jack_val_fit_fresh)$r.squared,
                               RMSE=RMSD(meta(val_jack_fresh)$EWT,EWT_jack_val_pred_fresh),
                               perRMSE=percentRMSD(meta(val_jack_fresh)$EWT,EWT_jack_val_pred_fresh,0.025,0.975),
                               bias=mean(EWT_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$EWT,na.rm=T))
  
  chlA_jack_val_pred_fresh<-as.vector(predict(chlA_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_chlA_fresh)[,,1])
  chlA_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$chlA~chlA_jack_val_pred_fresh)
  chlA_jack_stats_fresh[[i]]<-c(R2=summary(chlA_jack_val_fit_fresh)$r.squared,
                               RMSE=RMSD(meta(val_jack_fresh)$chlA,chlA_jack_val_pred_fresh),
                               perRMSE=percentRMSD(meta(val_jack_fresh)$chlA,chlA_jack_val_pred_fresh,0.025,0.975),
                               bias=mean(chlA_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$chlA,na.rm=T))
  
  chlB_jack_val_pred_fresh<-as.vector(predict(chlB_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_chlB_fresh)[,,1])
  chlB_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$chlB~chlB_jack_val_pred_fresh)
  chlB_jack_stats_fresh[[i]]<-c(R2=summary(chlB_jack_val_fit_fresh)$r.squared,
                               RMSE=RMSD(meta(val_jack_fresh)$chlB,chlB_jack_val_pred_fresh),
                               perRMSE=percentRMSD(meta(val_jack_fresh)$chlB,chlB_jack_val_pred_fresh,0.025,0.975),
                               bias=mean(chlB_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$chlB,na.rm=T))

  car_jack_val_pred_fresh<-as.vector(predict(car_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_car_fresh)[,,1])
  car_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$car~car_jack_val_pred_fresh)
  car_jack_stats_fresh[[i]]<-c(R2=summary(car_jack_val_fit_fresh)$r.squared,
                               RMSE=RMSD(meta(val_jack_fresh)$car,car_jack_val_pred_fresh),
                               perRMSE=percentRMSD(meta(val_jack_fresh)$car,car_jack_val_pred_fresh,0.025,0.975),
                               bias=mean(car_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$car,na.rm=T))
  
  Al_jack_val_pred_fresh<-as.vector(predict(Al_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Al_fresh)[,,1])
  Al_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Al~Al_jack_val_pred_fresh)
  Al_jack_stats_fresh[[i]]<-c(R2=summary(Al_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Al,Al_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Al,Al_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Al_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Al,na.rm=T))
  
  Ca_jack_val_pred_fresh<-as.vector(predict(Ca_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Ca_fresh)[,,1])
  Ca_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Ca~Ca_jack_val_pred_fresh)
  Ca_jack_stats_fresh[[i]]<-c(R2=summary(Ca_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Ca,Ca_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Ca,Ca_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Ca_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Ca,na.rm=T))
  
  Cu_jack_val_pred_fresh<-as.vector(predict(Cu_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Cu_fresh)[,,1])
  Cu_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Cu~Cu_jack_val_pred_fresh)
  Cu_jack_stats_fresh[[i]]<-c(R2=summary(Cu_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Cu,Cu_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Cu,Cu_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Cu_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Cu,na.rm=T))
  
  Fe_jack_val_pred_fresh<-as.vector(predict(Fe_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Fe_fresh)[,,1])
  Fe_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Fe~Fe_jack_val_pred_fresh)
  Fe_jack_stats_fresh[[i]]<-c(R2=summary(Fe_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Fe,Fe_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Fe,Fe_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Fe_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Fe,na.rm=T))
  
  K_jack_val_pred_fresh<-as.vector(predict(K_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_K_fresh)[,,1])
  K_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$K~K_jack_val_pred_fresh)
  K_jack_stats_fresh[[i]]<-c(R2=summary(K_jack_val_fit_fresh)$r.squared,
                             RMSE=RMSD(meta(val_jack_fresh)$K,K_jack_val_pred_fresh),
                             perRMSE=percentRMSD(meta(val_jack_fresh)$K,K_jack_val_pred_fresh,0.025,0.975),
                             bias=mean(K_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$K,na.rm=T))
  
  Mg_jack_val_pred_fresh<-as.vector(predict(Mg_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Mg_fresh)[,,1])
  Mg_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Mg~Mg_jack_val_pred_fresh)
  Mg_jack_stats_fresh[[i]]<-c(R2=summary(Mg_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Mg,Mg_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Mg,Mg_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Mg_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Mg,na.rm=T))
  
  Mn_jack_val_pred_fresh<-as.vector(predict(Mn_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Mn_fresh)[,,1])
  Mn_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Mn~Mn_jack_val_pred_fresh)
  Mn_jack_stats_fresh[[i]]<-c(R2=summary(Mn_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Mn,Mn_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Mn,Mn_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Mn_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Mn,na.rm=T))
  
  Na_jack_val_pred_fresh<-as.vector(predict(Na_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Na_fresh)[,,1])
  Na_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Na~Na_jack_val_pred_fresh)
  Na_jack_stats_fresh[[i]]<-c(R2=summary(Na_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Na,Na_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Na,Na_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Na_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Na,na.rm=T))
  
  P_jack_val_pred_fresh<-as.vector(predict(P_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_P_fresh)[,,1])
  P_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$P~P_jack_val_pred_fresh)
  P_jack_stats_fresh[[i]]<-c(R2=summary(P_jack_val_fit_fresh)$r.squared,
                             RMSE=RMSD(meta(val_jack_fresh)$P,P_jack_val_pred_fresh),
                             perRMSE=percentRMSD(meta(val_jack_fresh)$P,P_jack_val_pred_fresh,0.025,0.975),
                             bias=mean(P_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$P,na.rm=T))
  
  Zn_jack_val_pred_fresh<-as.vector(predict(Zn_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Zn_fresh)[,,1])
  Zn_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Zn~Zn_jack_val_pred_fresh)
  Zn_jack_stats_fresh[[i]]<-c(R2=summary(Zn_jack_val_fit_fresh)$r.squared,
                              RMSE=RMSD(meta(val_jack_fresh)$Zn,Zn_jack_val_pred_fresh),
                              perRMSE=percentRMSD(meta(val_jack_fresh)$Zn,Zn_jack_val_pred_fresh,0.025,0.975),
                              bias=mean(Zn_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Zn,na.rm=T))
  
  solubles_jack_coefs_fresh[[i]]<-as.vector(coef(solubles_fresh_jack,ncomp=ncomp_solubles_fresh,intercept=TRUE))
  hemicellulose_jack_coefs_fresh[[i]]<-as.vector(coef(hemicellulose_fresh_jack,ncomp=ncomp_hemicellulose_fresh,intercept=TRUE))
  cellulose_jack_coefs_fresh[[i]]<-as.vector(coef(cellulose_fresh_jack,ncomp=ncomp_cellulose_fresh,intercept=TRUE))
  lignin_jack_coefs_fresh[[i]]<-as.vector(coef(lignin_fresh_jack,ncomp=ncomp_lignin_fresh,intercept=TRUE))
  perC_jack_coefs_fresh[[i]]<-as.vector(coef(perC_fresh_jack,ncomp=ncomp_perC_fresh,intercept=TRUE))
  perN_jack_coefs_fresh[[i]]<-as.vector(coef(perN_fresh_jack,ncomp=ncomp_perN_fresh,intercept=TRUE))
  LMA_jack_coefs_fresh[[i]]<-as.vector(coef(LMA_fresh_jack,ncomp=ncomp_LMA_fresh,intercept=TRUE))
  LDMC_jack_coefs_fresh[[i]]<-as.vector(coef(LDMC_fresh_jack,ncomp=ncomp_LDMC_fresh,intercept=TRUE))
  EWT_jack_coefs_fresh[[i]]<-as.vector(coef(EWT_fresh_jack,ncomp=ncomp_EWT_fresh,intercept=TRUE))
  chlA_jack_coefs_fresh[[i]]<-as.vector(coef(chlA_fresh_jack,ncomp=ncomp_chlA_fresh,intercept=TRUE))
  chlB_jack_coefs_fresh[[i]]<-as.vector(coef(chlB_fresh_jack,ncomp=ncomp_chlB_fresh,intercept=TRUE))
  car_jack_coefs_fresh[[i]]<-as.vector(coef(car_fresh_jack,ncomp=ncomp_car_fresh,intercept=TRUE))
  Al_jack_coefs_fresh[[i]]<-as.vector(coef(Al_fresh_jack,ncomp=ncomp_Al_fresh,intercept=TRUE))
  Ca_jack_coefs_fresh[[i]]<-as.vector(coef(Ca_fresh_jack,ncomp=ncomp_Ca_fresh,intercept=TRUE))
  Cu_jack_coefs_fresh[[i]]<-as.vector(coef(Cu_fresh_jack,ncomp=ncomp_Cu_fresh,intercept=TRUE))
  Fe_jack_coefs_fresh[[i]]<-as.vector(coef(Fe_fresh_jack,ncomp=ncomp_Fe_fresh,intercept=TRUE))
  K_jack_coefs_fresh[[i]]<-as.vector(coef(K_fresh_jack,ncomp=ncomp_K_fresh,intercept=TRUE))
  Mg_jack_coefs_fresh[[i]]<-as.vector(coef(Mg_fresh_jack,ncomp=ncomp_Mg_fresh,intercept=TRUE))
  Mn_jack_coefs_fresh[[i]]<-as.vector(coef(Mn_fresh_jack,ncomp=ncomp_Mn_fresh,intercept=TRUE))
  Na_jack_coefs_fresh[[i]]<-as.vector(coef(Na_fresh_jack,ncomp=ncomp_Na_fresh,intercept=TRUE))
  P_jack_coefs_fresh[[i]]<-as.vector(coef(P_fresh_jack,ncomp=ncomp_P_fresh,intercept=TRUE))
  Zn_jack_coefs_fresh[[i]]<-as.vector(coef(Zn_fresh_jack,ncomp=ncomp_Zn_fresh,intercept=TRUE))
}

solubles_jack_pred_fresh<-apply.coefs(solubles_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
solubles_jack_stat_fresh<-t(apply(solubles_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_jack_df_fresh<-data.frame(pred_mean=solubles_jack_stat_fresh[,1],
                                   pred_low=solubles_jack_stat_fresh[,2],
                                   pred_high=solubles_jack_stat_fresh[,3],
                                   Measured=meta(fresh_spec_EL_agg_test)$solubles,
                                   ncomp=ncomp_solubles_fresh,
                                   Project=meta(fresh_spec_EL_agg_test)$Project,
                                   GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                                   Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                                   ID=meta(fresh_spec_EL_agg_test)$ID)

hemicellulose_jack_pred_fresh<-apply.coefs(hemicellulose_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
hemicellulose_jack_stat_fresh<-t(apply(hemicellulose_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_jack_df_fresh<-data.frame(pred_mean=hemicellulose_jack_stat_fresh[,1],
                                        pred_low=hemicellulose_jack_stat_fresh[,2],
                                        pred_high=hemicellulose_jack_stat_fresh[,3],
                                        Measured=meta(fresh_spec_EL_agg_test)$hemicellulose,
                                        ncomp=ncomp_hemicellulose_fresh,
                                        Project=meta(fresh_spec_EL_agg_test)$Project,
                                        GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                                        Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                                        ID=meta(fresh_spec_EL_agg_test)$ID)

cellulose_jack_pred_fresh<-apply.coefs(cellulose_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
cellulose_jack_stat_fresh<-t(apply(cellulose_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_jack_df_fresh<-data.frame(pred_mean=cellulose_jack_stat_fresh[,1],
                                    pred_low=cellulose_jack_stat_fresh[,2],
                                    pred_high=cellulose_jack_stat_fresh[,3],
                                    Measured=meta(fresh_spec_EL_agg_test)$cellulose,
                                    ncomp=ncomp_cellulose_fresh,
                                    Project=meta(fresh_spec_EL_agg_test)$Project,
                                    GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                                    Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                                    ID=meta(fresh_spec_EL_agg_test)$ID)

lignin_jack_pred_fresh<-apply.coefs(lignin_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
lignin_jack_stat_fresh<-t(apply(lignin_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_jack_df_fresh<-data.frame(pred_mean=lignin_jack_stat_fresh[,1],
                                 pred_low=lignin_jack_stat_fresh[,2],
                                 pred_high=lignin_jack_stat_fresh[,3],
                                 Measured=meta(fresh_spec_EL_agg_test)$lignin,
                                 ncomp=ncomp_lignin_fresh,
                                 Project=meta(fresh_spec_EL_agg_test)$Project,
                                 GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                                 Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                                 ID=meta(fresh_spec_EL_agg_test)$ID)

perC_jack_pred_fresh<-apply.coefs(perC_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
perC_jack_stat_fresh<-t(apply(perC_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_jack_df_fresh<-data.frame(pred_mean=perC_jack_stat_fresh[,1],
                               pred_low=perC_jack_stat_fresh[,2],
                               pred_high=perC_jack_stat_fresh[,3],
                               Measured=meta(fresh_spec_EL_agg_test)$C,
                               ncomp=ncomp_perC_fresh,
                               Project=meta(fresh_spec_EL_agg_test)$Project,
                               GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                               ID=meta(fresh_spec_EL_agg_test)$ID)

perN_jack_pred_fresh<-apply.coefs(perN_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
perN_jack_stat_fresh<-t(apply(perN_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_jack_df_fresh<-data.frame(pred_mean=perN_jack_stat_fresh[,1],
                               pred_low=perN_jack_stat_fresh[,2],
                               pred_high=perN_jack_stat_fresh[,3],
                               Measured=meta(fresh_spec_EL_agg_test)$N,
                               ncomp=ncomp_perN_fresh,
                               Project=meta(fresh_spec_EL_agg_test)$Project,
                               GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                               ID=meta(fresh_spec_EL_agg_test)$ID)

LMA_jack_pred_fresh<-apply.coefs(LMA_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
LMA_jack_stat_fresh<-t(apply(LMA_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df_fresh<-data.frame(pred_mean=LMA_jack_stat_fresh[,1],
                              pred_low=LMA_jack_stat_fresh[,2],
                              pred_high=LMA_jack_stat_fresh[,3],
                              Measured=meta(fresh_spec_EL_agg_test)$LMA,
                              ncomp=ncomp_LMA_fresh,
                              Project=meta(fresh_spec_EL_agg_test)$Project,
                              GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                              ID=meta(fresh_spec_EL_agg_test)$ID)

LDMC_jack_pred_fresh<-apply.coefs(LDMC_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
LDMC_jack_stat_fresh<-t(apply(LDMC_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC_jack_df_fresh<-data.frame(pred_mean=LDMC_jack_stat_fresh[,1],
                               pred_low=LDMC_jack_stat_fresh[,2],
                               pred_high=LDMC_jack_stat_fresh[,3],
                               Measured=meta(fresh_spec_EL_agg_test)$LDMC,
                               ncomp=ncomp_LDMC_fresh,
                               Project=meta(fresh_spec_EL_agg_test)$Project,
                               GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                               ID=meta(fresh_spec_EL_agg_test)$ID)

EWT_jack_pred_fresh<-apply.coefs(EWT_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
EWT_jack_stat_fresh<-t(apply(EWT_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_jack_df_fresh<-data.frame(pred_mean=EWT_jack_stat_fresh[,1],
                              pred_low=EWT_jack_stat_fresh[,2],
                              pred_high=EWT_jack_stat_fresh[,3],
                              Measured=meta(fresh_spec_EL_agg_test)$EWT,
                              ncomp=ncomp_EWT_fresh,
                              Project=meta(fresh_spec_EL_agg_test)$Project,
                              GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                              ID=meta(fresh_spec_EL_agg_test)$ID)

chlA_jack_pred_fresh<-apply.coefs(chlA_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
chlA_jack_stat_fresh<-t(apply(chlA_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_jack_df_fresh<-data.frame(pred_mean=chlA_jack_stat_fresh[,1],
                              pred_low=chlA_jack_stat_fresh[,2],
                              pred_high=chlA_jack_stat_fresh[,3],
                              Measured=meta(fresh_spec_EL_agg_test)$chlA,
                              ncomp=ncomp_chlA_fresh,
                              Project=meta(fresh_spec_EL_agg_test)$Project,
                              GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                              ID=meta(fresh_spec_EL_agg_test)$ID)

chlB_jack_pred_fresh<-apply.coefs(chlB_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
chlB_jack_stat_fresh<-t(apply(chlB_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_jack_df_fresh<-data.frame(pred_mean=chlB_jack_stat_fresh[,1],
                               pred_low=chlB_jack_stat_fresh[,2],
                               pred_high=chlB_jack_stat_fresh[,3],
                               Measured=meta(fresh_spec_EL_agg_test)$chlB,
                               ncomp=ncomp_chlB_fresh,
                               Project=meta(fresh_spec_EL_agg_test)$Project,
                               GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                               Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                               ID=meta(fresh_spec_EL_agg_test)$ID)

car_jack_pred_fresh<-apply.coefs(car_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
car_jack_stat_fresh<-t(apply(car_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_jack_df_fresh<-data.frame(pred_mean=car_jack_stat_fresh[,1],
                              pred_low=car_jack_stat_fresh[,2],
                              pred_high=car_jack_stat_fresh[,3],
                              Measured=meta(fresh_spec_EL_agg_test)$car,
                              ncomp=ncomp_car_fresh,
                              Project=meta(fresh_spec_EL_agg_test)$Project,
                              GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                              Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                              ID=meta(fresh_spec_EL_agg_test)$ID)

Al_jack_pred_fresh<-apply.coefs(Al_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Al_jack_stat_fresh<-t(apply(Al_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_jack_df_fresh<-data.frame(pred_mean=Al_jack_stat_fresh[,1],
                             pred_low=Al_jack_stat_fresh[,2],
                             pred_high=Al_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Al,
                             ncomp=ncomp_Al_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

Ca_jack_pred_fresh<-apply.coefs(Ca_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Ca_jack_stat_fresh<-t(apply(Ca_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_jack_df_fresh<-data.frame(pred_mean=Ca_jack_stat_fresh[,1],
                             pred_low=Ca_jack_stat_fresh[,2],
                             pred_high=Ca_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Ca,
                             ncomp=ncomp_Ca_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

Cu_jack_pred_fresh<-apply.coefs(Cu_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Cu_jack_stat_fresh<-t(apply(Cu_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_jack_df_fresh<-data.frame(pred_mean=Cu_jack_stat_fresh[,1],
                             pred_low=Cu_jack_stat_fresh[,2],
                             pred_high=Cu_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Cu,
                             ncomp=ncomp_Cu_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

Fe_jack_pred_fresh<-apply.coefs(Fe_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Fe_jack_stat_fresh<-t(apply(Fe_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_jack_df_fresh<-data.frame(pred_mean=Fe_jack_stat_fresh[,1],
                             pred_low=Fe_jack_stat_fresh[,2],
                             pred_high=Fe_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Fe,
                             ncomp=ncomp_Fe_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

K_jack_pred_fresh<-apply.coefs(K_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
K_jack_stat_fresh<-t(apply(K_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_jack_df_fresh<-data.frame(pred_mean=K_jack_stat_fresh[,1],
                            pred_low=K_jack_stat_fresh[,2],
                            pred_high=K_jack_stat_fresh[,3],
                            Measured=meta(fresh_spec_EL_agg_test)$K,
                            ncomp=ncomp_K_fresh,
                            Project=meta(fresh_spec_EL_agg_test)$Project,
                            GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                            Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                            ID=meta(fresh_spec_EL_agg_test)$ID)

Mg_jack_pred_fresh<-apply.coefs(Mg_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Mg_jack_stat_fresh<-t(apply(Mg_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_jack_df_fresh<-data.frame(pred_mean=Mg_jack_stat_fresh[,1],
                             pred_low=Mg_jack_stat_fresh[,2],
                             pred_high=Mg_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Mg,
                             ncomp=ncomp_Mg_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

Mn_jack_pred_fresh<-apply.coefs(Mn_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Mn_jack_stat_fresh<-t(apply(Mn_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_jack_df_fresh<-data.frame(pred_mean=Mn_jack_stat_fresh[,1],
                             pred_low=Mn_jack_stat_fresh[,2],
                             pred_high=Mn_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Mn,
                             ncomp=ncomp_Mn_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

Na_jack_pred_fresh<-apply.coefs(Na_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Na_jack_stat_fresh<-t(apply(Na_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_jack_df_fresh<-data.frame(pred_mean=Na_jack_stat_fresh[,1],
                             pred_low=Na_jack_stat_fresh[,2],
                             pred_high=Na_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Na,
                             ncomp=ncomp_Na_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

P_jack_pred_fresh<-apply.coefs(P_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
P_jack_stat_fresh<-t(apply(P_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_jack_df_fresh<-data.frame(pred_mean=P_jack_stat_fresh[,1],
                            pred_low=P_jack_stat_fresh[,2],
                            pred_high=P_jack_stat_fresh[,3],
                            Measured=meta(fresh_spec_EL_agg_test)$P,
                            ncomp=ncomp_P_fresh,
                            Project=meta(fresh_spec_EL_agg_test)$Project,
                            GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                            Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                            ID=meta(fresh_spec_EL_agg_test)$ID)

Zn_jack_pred_fresh<-apply.coefs(Zn_jack_coefs_fresh,as.matrix(fresh_spec_EL_agg_test))
Zn_jack_stat_fresh<-t(apply(Zn_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_jack_df_fresh<-data.frame(pred_mean=Zn_jack_stat_fresh[,1],
                             pred_low=Zn_jack_stat_fresh[,2],
                             pred_high=Zn_jack_stat_fresh[,3],
                             Measured=meta(fresh_spec_EL_agg_test)$Zn,
                             ncomp=ncomp_Zn_fresh,
                             Project=meta(fresh_spec_EL_agg_test)$Project,
                             GrowthForm=meta(fresh_spec_EL_agg_test)$GrowthForm,
                             Discoloration=meta(fresh_spec_EL_agg_test)$Discoloration,
                             ID=meta(fresh_spec_EL_agg_test)$ID)

##########################################################
## save output

fresh_jack_coef_list<-list(LMA=LMA_jack_coefs_fresh,
                             LDMC=LDMC_jack_coefs_fresh,
                             EWT=EWT_jack_coefs_fresh,
                             C=perC_jack_coefs_fresh,
                             N=perN_jack_coefs_fresh,
                             sol=solubles_jack_coefs_fresh,
                             hemi=hemicellulose_jack_coefs_fresh,
                             cell=cellulose_jack_coefs_fresh,
                             lign=lignin_jack_coefs_fresh,
                             chlA=chlA_jack_coefs_fresh,
                             chlB=chlB_jack_coefs_fresh,
                             car=car_jack_coefs_fresh,
                             Al=Al_jack_coefs_fresh,
                             Ca=Ca_jack_coefs_fresh,
                             Cu=Cu_jack_coefs_fresh,
                             Fe=Fe_jack_coefs_fresh,
                             K=K_jack_coefs_fresh,
                             Mg=Mg_jack_coefs_fresh,
                             Mn=Mn_jack_coefs_fresh,
                             Na=Na_jack_coefs_fresh,
                             P=P_jack_coefs_fresh,
                             Zn=Zn_jack_coefs_fresh)
saveRDS(fresh_jack_coef_list,"SavedResults/fresh_jack_coefs_list.rds")

fresh_jack_df_list<-list(LMA=LMA_jack_df_fresh,
                           LDMC=LDMC_jack_df_fresh,
                           EWT=EWT_jack_df_fresh,
                           C=perC_jack_df_fresh,
                           N=perN_jack_df_fresh,
                           sol=solubles_jack_df_fresh,
                           hemi=hemicellulose_jack_df_fresh,
                           cell=cellulose_jack_df_fresh,
                           lign=lignin_jack_df_fresh,
                           chlA=chlA_jack_df_fresh,
                           chlB=chlB_jack_df_fresh,
                           car=car_jack_df_fresh,
                           Al=Al_jack_df_fresh,
                           Ca=Ca_jack_df_fresh,
                           Cu=Cu_jack_df_fresh,
                           Fe=Fe_jack_df_fresh,
                           K=K_jack_df_fresh,
                           Mg=Mg_jack_df_fresh,
                           Mn=Mn_jack_df_fresh,
                           Na=Na_jack_df_fresh,
                           P=P_jack_df_fresh,
                           Zn=Zn_jack_df_fresh)
saveRDS(fresh_jack_df_list,"SavedResults/fresh_jack_df_list.rds")

####################################################
## violin plots

fresh_val_summary<-data.frame(variable=names(fresh_jack_df_list),
                                perRMSE=unlist(lapply(fresh_jack_df_list,
                                                      function(x) percentRMSD(x$Measured,x$pred_mean,0.025,0.975))),
                                R2=unlist(lapply(fresh_jack_df_list,
                                                 function(x) summary(lm(Measured~pred_mean,data=x))$r.squared)))

R2.df_fresh<-data.frame(LMA=unlist(lapply(LMA_jack_stats_fresh,function(x) x[["R2"]])),
                        LDMC=unlist(lapply(LDMC_jack_stats_fresh,function(x) x[["R2"]])),
                        EWT=unlist(lapply(EWT_jack_stats_fresh,function(x) x[["R2"]])),
                        sol=unlist(lapply(solubles_jack_stats_fresh,function(x) x[["R2"]])),
                        hemi=unlist(lapply(hemicellulose_jack_stats_fresh,function(x) x[["R2"]])),
                        cell=unlist(lapply(cellulose_jack_stats_fresh,function(x) x[["R2"]])),
                        lign=unlist(lapply(lignin_jack_stats_fresh,function(x) x[["R2"]])),
                        C=unlist(lapply(perC_jack_stats_fresh,function(x) x[["R2"]])),
                        N=unlist(lapply(perN_jack_stats_fresh,function(x) x[["R2"]])),
                        chlA=unlist(lapply(chlA_jack_stats_fresh,function(x) x[["R2"]])),
                        chlB=unlist(lapply(chlB_jack_stats_fresh,function(x) x[["R2"]])),
                        car=unlist(lapply(car_jack_stats_fresh,function(x) x[["R2"]])),
                        Al=unlist(lapply(Al_jack_stats_fresh,function(x) x[["R2"]])),
                        Ca=unlist(lapply(Ca_jack_stats_fresh,function(x) x[["R2"]])),
                        Cu=unlist(lapply(Cu_jack_stats_fresh,function(x) x[["R2"]])),
                        Fe=unlist(lapply(Fe_jack_stats_fresh,function(x) x[["R2"]])),
                        K=unlist(lapply(K_jack_stats_fresh,function(x) x[["R2"]])),
                        Mg=unlist(lapply(Mg_jack_stats_fresh,function(x) x[["R2"]])),
                        Mn=unlist(lapply(Mn_jack_stats_fresh,function(x) x[["R2"]])),
                        Na=unlist(lapply(Na_jack_stats_fresh,function(x) x[["R2"]])),
                        P=unlist(lapply(P_jack_stats_fresh,function(x) x[["R2"]])),
                        Zn=unlist(lapply(Zn_jack_stats_fresh,function(x) x[["R2"]])))

R2.long_fresh<-melt(R2.df_fresh)
fresh_val_R2<-ggplot(R2.long_fresh,aes(y=value,x=variable))+
  geom_violin()+
  geom_point(data=fresh_val_summary,
             aes(x=variable,y=R2),color="red",size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y=expression(italic("R"^2)))+
  ggtitle("Fresh-leaf spectra")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))

perRMSE.df_fresh<-data.frame(LMA=unlist(lapply(LMA_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             LDMC=unlist(lapply(LDMC_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             EWT=unlist(lapply(EWT_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             sol=unlist(lapply(solubles_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             hemi=unlist(lapply(hemicellulose_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             cell=unlist(lapply(cellulose_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             lign=unlist(lapply(lignin_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             C=unlist(lapply(perC_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             N=unlist(lapply(perN_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             chlA=unlist(lapply(chlA_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             chlB=unlist(lapply(chlB_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             car=unlist(lapply(car_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Al=unlist(lapply(Al_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Ca=unlist(lapply(Ca_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Cu=unlist(lapply(Cu_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Fe=unlist(lapply(Fe_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             K=unlist(lapply(K_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Mg=unlist(lapply(Mg_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Mn=unlist(lapply(Mn_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Na=unlist(lapply(Na_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             P=unlist(lapply(P_jack_stats_fresh,function(x) 100*x[["perRMSE"]])),
                             Zn=unlist(lapply(Zn_jack_stats_fresh,function(x) 100*x[["perRMSE"]])))

perRMSE.long_fresh<-melt(perRMSE.df_fresh)
fresh_val_perRMSE<-ggplot(perRMSE.long_fresh,aes(y=value,x=variable))+
  geom_violin()+
  geom_point(data=fresh_val_summary,
             aes(x=variable,y=perRMSE*100),color="red",size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  labs(y="%RMSE")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,max(perRMSE.long_fresh$value)*1.1))

pdf("Manuscript/FigS11.pdf",width=8,height=6)
egg::ggarrange(plots = list(fresh_val_R2,fresh_val_perRMSE),
               nrow=2,ncol=1)
dev.off()
