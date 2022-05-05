setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(pls)
library(ggplot2)
library(caret)
library(reshape2)
library(RColorBrewer)
library(patchwork)
source("Scripts/pressed-leaf-models/00 useful_functions.R")

######################################################
## read data

pressed_spec_all_train<-readRDS("ProcessedSpectralData/pressed_spec_all_train.rds")
pressed_spec_all_test<-readRDS("ProcessedSpectralData/pressed_spec_all_test.rds")
pressed_1300_spec_all_train<-pressed_spec_all_train[,1300:2400]
pressed_1300_spec_all_test<-pressed_spec_all_test[,1300:2400]

###################################################
## building calibration models

perC_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$C~as.matrix(pressed_1300_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perC_pressed_1300 <- selectNcomp(perC_pressed_1300, method = "onesigma", plot = FALSE)
perC_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$N))
perC_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[perC_pressed_1300_valid],
                              Species=meta(pressed_1300_spec_all_train)$Species[perC_pressed_1300_valid],
                              Project=meta(pressed_1300_spec_all_train)$Project[perC_pressed_1300_valid],
                              Stage=meta(pressed_1300_spec_all_train)$Stage[perC_pressed_1300_valid],
                              GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[perC_pressed_1300_valid],
                              measured=meta(pressed_1300_spec_all_train)$C[perC_pressed_1300_valid],
                              val_pred=perC_pressed_1300$validation$pred[,,ncomp_perC_pressed_1300])
ggplot(perC_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %C from pressed-leaf spectra")

perN_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$N~as.matrix(pressed_1300_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perN_pressed_1300 <- selectNcomp(perN_pressed_1300, method = "onesigma", plot = FALSE)
perN_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$N))
perN_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[perN_pressed_1300_valid],
                              Species=meta(pressed_1300_spec_all_train)$Species[perN_pressed_1300_valid],
                              Project=meta(pressed_1300_spec_all_train)$Project[perN_pressed_1300_valid],
                              Stage=meta(pressed_1300_spec_all_train)$Stage[perN_pressed_1300_valid],
                              GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[perN_pressed_1300_valid],
                              measured=meta(pressed_1300_spec_all_train)$N[perN_pressed_1300_valid],
                              val_pred=perN_pressed_1300$validation$pred[,,ncomp_perN_pressed_1300])
ggplot(perN_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %N from pressed-leaf spectra")


LMA_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$LMA~as.matrix(pressed_1300_spec_all_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_pressed_1300 <- selectNcomp(LMA_pressed_1300, method = "onesigma", plot = FALSE)
LMA_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$LMA))
LMA_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[LMA_pressed_1300_valid],
                             Species=meta(pressed_1300_spec_all_train)$Species[LMA_pressed_1300_valid],
                             Project=meta(pressed_1300_spec_all_train)$Project[LMA_pressed_1300_valid],
                             Stage=meta(pressed_1300_spec_all_train)$Stage[LMA_pressed_1300_valid],
                             GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[LMA_pressed_1300_valid],
                             measured=meta(pressed_1300_spec_all_train)$LMA[LMA_pressed_1300_valid],
                             val_pred=LMA_pressed_1300$validation$pred[,,ncomp_LMA_pressed_1300])
ggplot(LMA_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.25))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting LMA from pressed-leaf spectra")

LDMC_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$LDMC~as.matrix(pressed_1300_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_pressed_1300 <- selectNcomp(LDMC_pressed_1300, method = "onesigma", plot = FALSE)
LDMC_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$LDMC))
LDMC_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[LDMC_pressed_1300_valid],
                              Species=meta(pressed_1300_spec_all_train)$Species[LDMC_pressed_1300_valid],
                              Project=meta(pressed_1300_spec_all_train)$Project[LDMC_pressed_1300_valid],
                              Stage=meta(pressed_1300_spec_all_train)$Stage[LDMC_pressed_1300_valid],
                              GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[LDMC_pressed_1300_valid],
                              measured=meta(pressed_1300_spec_all_train)$LDMC[LDMC_pressed_1300_valid],
                              val_pred=LDMC_pressed_1300$validation$pred[,,ncomp_LDMC_pressed_1300])
ggplot(LDMC_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(150,600),ylim=c(150,600))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting LDMC from pressed-leaf spectra")

EWT_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$EWT~as.matrix(pressed_1300_spec_all_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_pressed_1300 <- selectNcomp(EWT_pressed_1300, method = "onesigma", plot = FALSE)
EWT_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$EWT))
EWT_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[EWT_pressed_1300_valid],
                             Species=meta(pressed_1300_spec_all_train)$Species[EWT_pressed_1300_valid],
                             Project=meta(pressed_1300_spec_all_train)$Project[EWT_pressed_1300_valid],
                             Stage=meta(pressed_1300_spec_all_train)$Stage[EWT_pressed_1300_valid],
                             GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[EWT_pressed_1300_valid],
                             measured=meta(pressed_1300_spec_all_train)$EWT[EWT_pressed_1300_valid],
                             val_pred=EWT_pressed_1300$validation$pred[,,ncomp_EWT_pressed_1300])
ggplot(EWT_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting EWT from pressed-leaf spectra")

chlA_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$chlA~as.matrix(pressed_1300_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_pressed_1300 <- selectNcomp(chlA_pressed_1300, method = "onesigma", plot = FALSE)
chlA_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$chlA))
chlA_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[chlA_pressed_1300_valid],
                              Species=meta(pressed_1300_spec_all_train)$Species[chlA_pressed_1300_valid],
                              Project=meta(pressed_1300_spec_all_train)$Project[chlA_pressed_1300_valid],
                              Stage=meta(pressed_1300_spec_all_train)$Stage[chlA_pressed_1300_valid],
                              GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[chlA_pressed_1300_valid],
                              measured=meta(pressed_1300_spec_all_train)$chlA[chlA_pressed_1300_valid],
                              val_pred=chlA_pressed_1300$validation$pred[,,ncomp_chlA_pressed_1300])
ggplot(chlA_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,16),ylim=c(0,16))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl a from pressed-leaf spectra")

chlB_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$chlB~as.matrix(pressed_1300_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_pressed_1300 <- selectNcomp(chlB_pressed_1300, method = "onesigma", plot = FALSE)
chlB_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$chlB))
chlB_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[chlB_pressed_1300_valid],
                              Species=meta(pressed_1300_spec_all_train)$Species[chlB_pressed_1300_valid],
                              Project=meta(pressed_1300_spec_all_train)$Project[chlB_pressed_1300_valid],
                              Stage=meta(pressed_1300_spec_all_train)$Stage[chlB_pressed_1300_valid],
                              GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[chlB_pressed_1300_valid],
                              measured=meta(pressed_1300_spec_all_train)$chlB[chlB_pressed_1300_valid],
                              val_pred=chlB_pressed_1300$validation$pred[,,ncomp_chlB_pressed_1300])
ggplot(chlB_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl b from pressed-leaf spectra")

car_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$car~as.matrix(pressed_1300_spec_all_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_pressed_1300 <- selectNcomp(car_pressed_1300, method = "onesigma", plot = FALSE)
car_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$car))
car_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[car_pressed_1300_valid],
                             Species=meta(pressed_1300_spec_all_train)$Species[car_pressed_1300_valid],
                             Project=meta(pressed_1300_spec_all_train)$Project[car_pressed_1300_valid],
                             Stage=meta(pressed_1300_spec_all_train)$Stage[car_pressed_1300_valid],
                             GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[car_pressed_1300_valid],
                             measured=meta(pressed_1300_spec_all_train)$car[car_pressed_1300_valid],
                             val_pred=car_pressed_1300$validation$pred[,,ncomp_car_pressed_1300])
ggplot(car_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,3.5),ylim=c(0,3.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting total car from pressed-leaf spectra")

solubles_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$solubles~as.matrix(pressed_1300_spec_all_train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_pressed_1300 <- selectNcomp(solubles_pressed_1300, method = "onesigma", plot = FALSE)
solubles_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$solubles))
solubles_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[solubles_pressed_1300_valid],
                                  Species=meta(pressed_1300_spec_all_train)$Species[solubles_pressed_1300_valid],
                                  Project=meta(pressed_1300_spec_all_train)$Project[solubles_pressed_1300_valid],
                                  Stage=meta(pressed_1300_spec_all_train)$Stage[solubles_pressed_1300_valid],
                                  GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[solubles_pressed_1300_valid],
                                  measured=meta(pressed_1300_spec_all_train)$solubles[solubles_pressed_1300_valid],
                                  val_pred=solubles_pressed_1300$validation$pred[,,ncomp_solubles_pressed_1300])
ggplot(solubles_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(35,90),ylim=c(35,90))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting solubles from pressed-leaf spectra")

hemicellulose_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$hemicellulose~as.matrix(pressed_1300_spec_all_train),
                            ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_pressed_1300 <- selectNcomp(hemicellulose_pressed_1300, method = "onesigma", plot = FALSE)
hemicellulose_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$hemicellulose))
hemicellulose_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[hemicellulose_pressed_1300_valid],
                                       Species=meta(pressed_1300_spec_all_train)$Species[hemicellulose_pressed_1300_valid],
                                       Project=meta(pressed_1300_spec_all_train)$Project[hemicellulose_pressed_1300_valid],
                                       Stage=meta(pressed_1300_spec_all_train)$Stage[hemicellulose_pressed_1300_valid],
                                       GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[hemicellulose_pressed_1300_valid],
                                       measured=meta(pressed_1300_spec_all_train)$hemicellulose[hemicellulose_pressed_1300_valid],
                                       val_pred=hemicellulose_pressed_1300$validation$pred[,,ncomp_hemicellulose_pressed_1300])
ggplot(hemicellulose_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,36),ylim=c(0,36))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting hemicellulose from pressed-leaf spectra")

cellulose_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$cellulose~as.matrix(pressed_1300_spec_all_train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_pressed_1300 <- selectNcomp(cellulose_pressed_1300, method = "onesigma", plot = FALSE)
cellulose_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$cellulose))
cellulose_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[cellulose_pressed_1300_valid],
                                   Species=meta(pressed_1300_spec_all_train)$Species[cellulose_pressed_1300_valid],
                                   Project=meta(pressed_1300_spec_all_train)$Project[cellulose_pressed_1300_valid],
                                   Stage=meta(pressed_1300_spec_all_train)$Stage[cellulose_pressed_1300_valid],
                                   GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[cellulose_pressed_1300_valid],
                                   measured=meta(pressed_1300_spec_all_train)$cellulose[cellulose_pressed_1300_valid],
                                   val_pred=cellulose_pressed_1300$validation$pred[,,ncomp_cellulose_pressed_1300])
ggplot(cellulose_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(2,28),ylim=c(2,28))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting cellulose from pressed-leaf spectra")

lignin_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$lignin~as.matrix(pressed_1300_spec_all_train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_pressed_1300 <- selectNcomp(lignin_pressed_1300, method = "onesigma", plot = FALSE)
lignin_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$lignin))
lignin_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[lignin_pressed_1300_valid],
                                Species=meta(pressed_1300_spec_all_train)$Species[lignin_pressed_1300_valid],
                                Project=meta(pressed_1300_spec_all_train)$Project[lignin_pressed_1300_valid],
                                Stage=meta(pressed_1300_spec_all_train)$Stage[lignin_pressed_1300_valid],
                                GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[lignin_pressed_1300_valid],
                                measured=meta(pressed_1300_spec_all_train)$lignin[lignin_pressed_1300_valid],
                                val_pred=lignin_pressed_1300$validation$pred[,,ncomp_lignin_pressed_1300])
ggplot(lignin_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,22),ylim=c(0,22))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting lignin from pressed-leaf spectra")


Al_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Al~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_pressed_1300 <- selectNcomp(Al_pressed_1300, method = "onesigma", plot = FALSE)
Al_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Al))
Al_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Al_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Al_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Al_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Al_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Al_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Al[Al_pressed_1300_valid],
                            val_pred=Al_pressed_1300$validation$pred[,,ncomp_Al_pressed_1300])
ggplot(Al_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.05,0.35),ylim=c(-0.05,0.35))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Al from pressed-leaf spectra")

Ca_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Ca~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_pressed_1300 <- selectNcomp(Ca_pressed_1300, method = "onesigma", plot = FALSE)
Ca_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Ca))
Ca_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Ca_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Ca_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Ca_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Ca_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Ca_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Ca[Ca_pressed_1300_valid],
                            val_pred=Ca_pressed_1300$validation$pred[,,ncomp_Ca_pressed_1300])
ggplot(Ca_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-10,40),ylim=c(-10,40))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Ca from pressed-leaf spectra")

Cu_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Cu~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_pressed_1300 <- selectNcomp(Cu_pressed_1300, method = "onesigma", plot = FALSE)
Cu_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Cu))
Cu_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Cu_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Cu_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Cu_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Cu_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Cu_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Cu[Cu_pressed_1300_valid],
                            val_pred=Cu_pressed_1300$validation$pred[,,ncomp_Cu_pressed_1300])
ggplot(Cu_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.005,0.055),ylim=c(-0.005,0.055))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Cu from pressed-leaf spectra")

Fe_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Fe~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_pressed_1300 <- selectNcomp(Fe_pressed_1300, method = "onesigma", plot = FALSE)
Fe_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Fe))
Fe_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Fe_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Fe_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Fe_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Fe_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Fe_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Fe[Fe_pressed_1300_valid],
                            val_pred=Fe_pressed_1300$validation$pred[,,ncomp_Fe_pressed_1300])
ggplot(Fe_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Fe from pressed-leaf spectra")

K_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$K~as.matrix(pressed_1300_spec_all_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_pressed_1300 <- selectNcomp(K_pressed_1300, method = "onesigma", plot = FALSE)
K_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$K))
K_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[K_pressed_1300_valid],
                           Species=meta(pressed_1300_spec_all_train)$Species[K_pressed_1300_valid],
                           Project=meta(pressed_1300_spec_all_train)$Project[K_pressed_1300_valid],
                           Stage=meta(pressed_1300_spec_all_train)$Stage[K_pressed_1300_valid],
                           GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[K_pressed_1300_valid],
                           measured=meta(pressed_1300_spec_all_train)$K[K_pressed_1300_valid],
                           val_pred=K_pressed_1300$validation$pred[,,ncomp_K_pressed_1300])
ggplot(K_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,35),ylim=c(0,35))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting K from pressed-leaf spectra")

Mg_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Mg~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_pressed_1300 <- selectNcomp(Mg_pressed_1300, method = "onesigma", plot = FALSE)
Mg_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Mg))
Mg_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Mg_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Mg_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Mg_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Mg_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Mg_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Mg[Mg_pressed_1300_valid],
                            val_pred=Mg_pressed_1300$validation$pred[,,ncomp_Mg_pressed_1300])
ggplot(Mg_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,8),ylim=c(0,8))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mg from pressed-leaf spectra")

Mn_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Mn~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_pressed_1300 <- selectNcomp(Mn_pressed_1300, method = "onesigma", plot = FALSE)
Mn_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Mn))
Mn_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Mn_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Mn_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Mn_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Mn_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Mn_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Mn[Mn_pressed_1300_valid],
                            val_pred=Mn_pressed_1300$validation$pred[,,ncomp_Mn_pressed_1300])
ggplot(Mn_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,1.1),ylim=c(-0.1,1.1))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mn from pressed-leaf spectra")

Na_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Na~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_pressed_1300 <- selectNcomp(Na_pressed_1300, method = "onesigma", plot = FALSE)
Na_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Na))
Na_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Na_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Na_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Na_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Na_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Na_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Na[Na_pressed_1300_valid],
                            val_pred=Na_pressed_1300$validation$pred[,,ncomp_Na_pressed_1300])
ggplot(Na_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,5),ylim=c(-0.5,5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Na from pressed-leaf spectra")

P_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$P~as.matrix(pressed_1300_spec_all_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_pressed_1300 <- selectNcomp(P_pressed_1300, method = "onesigma", plot = FALSE)
P_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$P))
P_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[P_pressed_1300_valid],
                           Species=meta(pressed_1300_spec_all_train)$Species[P_pressed_1300_valid],
                           Project=meta(pressed_1300_spec_all_train)$Project[P_pressed_1300_valid],
                           Stage=meta(pressed_1300_spec_all_train)$Stage[P_pressed_1300_valid],
                           GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[P_pressed_1300_valid],
                           measured=meta(pressed_1300_spec_all_train)$P[P_pressed_1300_valid],
                           val_pred=P_pressed_1300$validation$pred[,,ncomp_P_pressed_1300])
ggplot(P_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,8),ylim=c(0,8))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting P from pressed-leaf spectra")

Zn_pressed_1300<-plsr(meta(pressed_1300_spec_all_train)$Zn~as.matrix(pressed_1300_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_pressed_1300 <- selectNcomp(Zn_pressed_1300, method = "onesigma", plot = FALSE)
Zn_pressed_1300_valid <- which(!is.na(meta(pressed_1300_spec_all_train)$Zn))
Zn_pressed_1300_pred<-data.frame(ID=meta(pressed_1300_spec_all_train)$ID[Zn_pressed_1300_valid],
                            Species=meta(pressed_1300_spec_all_train)$Species[Zn_pressed_1300_valid],
                            Project=meta(pressed_1300_spec_all_train)$Project[Zn_pressed_1300_valid],
                            Stage=meta(pressed_1300_spec_all_train)$Stage[Zn_pressed_1300_valid],
                            GrowthForm=meta(pressed_1300_spec_all_train)$GrowthForm[Zn_pressed_1300_valid],
                            measured=meta(pressed_1300_spec_all_train)$Zn[Zn_pressed_1300_valid],
                            val_pred=Zn_pressed_1300$validation$pred[,,ncomp_Zn_pressed_1300])
ggplot(Zn_pressed_1300_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,0.7),ylim=c(-0.1,0.7))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Zn from pressed-leaf spectra")

###############################################
## VIP plots

source("VIP.R")

VIP_pressed_1300<-data.frame(LMA=VIP(LMA_pressed_1300)[ncomp_LMA_pressed_1300,],
                             LDMC=VIP(LDMC_pressed_1300)[ncomp_LDMC_pressed_1300,],
                             EWT=VIP(EWT_pressed_1300)[ncomp_EWT_pressed_1300,],
                             sol=VIP(solubles_pressed_1300)[ncomp_solubles_pressed_1300,],
                             hemi=VIP(hemicellulose_pressed_1300)[ncomp_hemicellulose_pressed_1300,],
                             cell=VIP(cellulose_pressed_1300)[ncomp_cellulose_pressed_1300,],
                             lign=VIP(lignin_pressed_1300)[ncomp_lignin_pressed_1300,],
                             chlA=VIP(chlA_pressed_1300)[ncomp_chlA_pressed_1300,],
                             chlB=VIP(chlB_pressed_1300)[ncomp_chlB_pressed_1300,],
                             car=VIP(car_pressed_1300)[ncomp_car_pressed_1300,],
                             C=VIP(perC_pressed_1300)[ncomp_perC_pressed_1300,],
                             N=VIP(perN_pressed_1300)[ncomp_perN_pressed_1300,],
                             Al=VIP(Al_pressed_1300)[ncomp_Al_pressed_1300,],
                             Ca=VIP(Ca_pressed_1300)[ncomp_Ca_pressed_1300,],
                             Cu=VIP(Cu_pressed_1300)[ncomp_Cu_pressed_1300,],
                             Fe=VIP(Fe_pressed_1300)[ncomp_Fe_pressed_1300,],
                             K=VIP(K_pressed_1300)[ncomp_K_pressed_1300,],
                             Mg=VIP(Mg_pressed_1300)[ncomp_Mg_pressed_1300,],
                             Mn=VIP(Mn_pressed_1300)[ncomp_Mn_pressed_1300,],
                             Na=VIP(Na_pressed_1300)[ncomp_Na_pressed_1300,],
                             P=VIP(P_pressed_1300)[ncomp_P_pressed_1300,],
                             Zn=VIP(Zn_pressed_1300)[ncomp_Zn_pressed_1300,],
                             wavelength=400:2400)

saveRDS(VIP_pressed_1300,"SavedResults/VIP_pressed_1300.rds")

#######################################
## validation: pressed leaves

solubles_jack_coefs_pressed_1300<-list()
hemicellulose_jack_coefs_pressed_1300<-list()
cellulose_jack_coefs_pressed_1300<-list()
lignin_jack_coefs_pressed_1300<-list()
perC_jack_coefs_pressed_1300<-list()
perN_jack_coefs_pressed_1300<-list()
LMA_jack_coefs_pressed_1300<-list()
LDMC_jack_coefs_pressed_1300<-list()
EWT_jack_coefs_pressed_1300<-list()
chlA_jack_coefs_pressed_1300<-list()
chlB_jack_coefs_pressed_1300<-list()
car_jack_coefs_pressed_1300<-list()
Al_jack_coefs_pressed_1300<-list()
Ca_jack_coefs_pressed_1300<-list()
Cu_jack_coefs_pressed_1300<-list()
Fe_jack_coefs_pressed_1300<-list()
K_jack_coefs_pressed_1300<-list()
Mg_jack_coefs_pressed_1300<-list()
Mn_jack_coefs_pressed_1300<-list()
Na_jack_coefs_pressed_1300<-list()
P_jack_coefs_pressed_1300<-list()
Zn_jack_coefs_pressed_1300<-list()

solubles_jack_stats_pressed_1300<-list()
hemicellulose_jack_stats_pressed_1300<-list()
cellulose_jack_stats_pressed_1300<-list()
lignin_jack_stats_pressed_1300<-list()
perC_jack_stats_pressed_1300<-list()
perN_jack_stats_pressed_1300<-list()
LMA_jack_stats_pressed_1300<-list()
LDMC_jack_stats_pressed_1300<-list()
EWT_jack_stats_pressed_1300<-list()
chlA_jack_stats_pressed_1300<-list()
chlB_jack_stats_pressed_1300<-list()
car_jack_stats_pressed_1300<-list()
Al_jack_stats_pressed_1300<-list()
Ca_jack_stats_pressed_1300<-list()
Cu_jack_stats_pressed_1300<-list()
Fe_jack_stats_pressed_1300<-list()
K_jack_stats_pressed_1300<-list()
Mg_jack_stats_pressed_1300<-list()
Mn_jack_stats_pressed_1300<-list()
Na_jack_stats_pressed_1300<-list()
P_jack_stats_pressed_1300<-list()
Zn_jack_stats_pressed_1300<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_pressed_1300<-nrow(pressed_1300_spec_all_train)
  train_jack_pressed_1300<-sample(1:n_cal_spec_pressed_1300,floor(0.7*n_cal_spec_pressed_1300))
  test_jack_pressed_1300<-setdiff(1:n_cal_spec_pressed_1300,train_jack_pressed_1300)
  
  calib_jack_pressed_1300<-pressed_1300_spec_all_train[train_jack_pressed_1300]
  val_jack_pressed_1300<-pressed_1300_spec_all_train[test_jack_pressed_1300]
  
  solubles_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$solubles~as.matrix(calib_jack_pressed_1300),
                              ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$hemicellulose~as.matrix(calib_jack_pressed_1300),
                                   ncomp=30,method = "oscorespls",validation="none")
  cellulose_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$cellulose~as.matrix(calib_jack_pressed_1300),
                               ncomp=30,method = "oscorespls",validation="none")
  lignin_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$lignin~as.matrix(calib_jack_pressed_1300),
                            ncomp=30,method = "oscorespls",validation="none")
  perC_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$C~as.matrix(calib_jack_pressed_1300),
                          ncomp=30,method = "oscorespls",validation="none")
  perN_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$N~as.matrix(calib_jack_pressed_1300),
                          ncomp=30,method = "oscorespls",validation="none")
  LMA_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$LMA~as.matrix(calib_jack_pressed_1300),
                         ncomp=30,method = "oscorespls",validation="none")
  LDMC_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$LDMC~as.matrix(calib_jack_pressed_1300),
                          ncomp=30,method = "oscorespls",validation="none")
  EWT_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$EWT~as.matrix(calib_jack_pressed_1300),
                         ncomp=30,method = "oscorespls",validation="none")
  chlA_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$chlA~as.matrix(calib_jack_pressed_1300),
                          ncomp=30,method = "oscorespls",validation="none")
  chlB_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$chlB~as.matrix(calib_jack_pressed_1300),
                          ncomp=30,method = "oscorespls",validation="none")
  car_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$car~as.matrix(calib_jack_pressed_1300),
                         ncomp=30,method = "oscorespls",validation="none")
  Al_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Al~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  Ca_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Ca~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  Cu_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Cu~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  Fe_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Fe~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  K_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$K~as.matrix(calib_jack_pressed_1300),
                       ncomp=30,method = "oscorespls",validation="none")
  Mg_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Mg~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  Mn_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Mn~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  Na_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Na~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  P_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$P~as.matrix(calib_jack_pressed_1300),
                       ncomp=30,method = "oscorespls",validation="none")
  Zn_pressed_1300_jack<-plsr(meta(calib_jack_pressed_1300)$Zn~as.matrix(calib_jack_pressed_1300),
                        ncomp=30,method = "oscorespls",validation="none")
  
  solubles_jack_val_pred_pressed_1300<-as.vector(predict(solubles_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_solubles_pressed_1300)[,,1])
  solubles_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$solubles~solubles_jack_val_pred_pressed_1300)
  solubles_jack_stats_pressed_1300[[i]]<-c(R2=summary(solubles_jack_val_fit_pressed_1300)$r.squared,
                                      RMSE=RMSD(meta(val_jack_pressed_1300)$solubles,solubles_jack_val_pred_pressed_1300),
                                      perRMSE=percentRMSD(meta(val_jack_pressed_1300)$solubles,solubles_jack_val_pred_pressed_1300,0.025,0.975),
                                      bias=mean(solubles_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$solubles,na.rm=T))
  
  hemicellulose_jack_val_pred_pressed_1300<-as.vector(predict(hemicellulose_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),
                                                         ncomp=ncomp_hemicellulose_pressed_1300)[,,1])
  hemicellulose_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$hemicellulose~hemicellulose_jack_val_pred_pressed_1300)
  hemicellulose_jack_stats_pressed_1300[[i]]<-c(R2=summary(hemicellulose_jack_val_fit_pressed_1300)$r.squared,
                                           RMSE=RMSD(meta(val_jack_pressed_1300)$hemicellulose,hemicellulose_jack_val_pred_pressed_1300),
                                           perRMSE=percentRMSD(meta(val_jack_pressed_1300)$hemicellulose,hemicellulose_jack_val_pred_pressed_1300,0.025,0.975),
                                           bias=mean(hemicellulose_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$hemicellulose,na.rm=T))
  
  cellulose_jack_val_pred_pressed_1300<-as.vector(predict(cellulose_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_cellulose_pressed_1300)[,,1])
  cellulose_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$cellulose~cellulose_jack_val_pred_pressed_1300)
  cellulose_jack_stats_pressed_1300[[i]]<-c(R2=summary(cellulose_jack_val_fit_pressed_1300)$r.squared,
                                       RMSE=RMSD(meta(val_jack_pressed_1300)$cellulose,cellulose_jack_val_pred_pressed_1300),
                                       perRMSE=percentRMSD(meta(val_jack_pressed_1300)$cellulose,cellulose_jack_val_pred_pressed_1300,0.025,0.975),
                                       bias=mean(cellulose_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$cellulose,na.rm=T))
  
  lignin_jack_val_pred_pressed_1300<-as.vector(predict(lignin_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_lignin_pressed_1300)[,,1])
  lignin_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$lignin~lignin_jack_val_pred_pressed_1300)
  lignin_jack_stats_pressed_1300[[i]]<-c(R2=summary(lignin_jack_val_fit_pressed_1300)$r.squared,
                                    RMSE=RMSD(meta(val_jack_pressed_1300)$lignin,lignin_jack_val_pred_pressed_1300),
                                    perRMSE=percentRMSD(meta(val_jack_pressed_1300)$lignin,lignin_jack_val_pred_pressed_1300,0.025,0.975),
                                    bias=mean(lignin_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$lignin,na.rm=T))
  
  perC_jack_val_pred_pressed_1300<-as.vector(predict(perC_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_perC_pressed_1300)[,,1])
  perC_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$C~perC_jack_val_pred_pressed_1300)
  perC_jack_stats_pressed_1300[[i]]<-c(R2=summary(perC_jack_val_fit_pressed_1300)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed_1300)$C,perC_jack_val_pred_pressed_1300),
                                  perRMSE=percentRMSD(meta(val_jack_pressed_1300)$C,perC_jack_val_pred_pressed_1300,0.025,0.975),
                                  bias=mean(perC_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$C,na.rm=T))
  
  perN_jack_val_pred_pressed_1300<-as.vector(predict(perN_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_perN_pressed_1300)[,,1])
  perN_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$N~perN_jack_val_pred_pressed_1300)
  perN_jack_stats_pressed_1300[[i]]<-c(R2=summary(perN_jack_val_fit_pressed_1300)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed_1300)$N,perN_jack_val_pred_pressed_1300),
                                  perRMSE=percentRMSD(meta(val_jack_pressed_1300)$N,perN_jack_val_pred_pressed_1300,0.025,0.975),
                                  bias=mean(perN_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$N,na.rm=T))
  
  LMA_jack_val_pred_pressed_1300<-as.vector(predict(LMA_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_LMA_pressed_1300)[,,1])
  LMA_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$LMA~LMA_jack_val_pred_pressed_1300)
  LMA_jack_stats_pressed_1300[[i]]<-c(R2=summary(LMA_jack_val_fit_pressed_1300)$r.squared,
                                 RMSE=RMSD(meta(val_jack_pressed_1300)$LMA,LMA_jack_val_pred_pressed_1300),
                                 perRMSE=percentRMSD(meta(val_jack_pressed_1300)$LMA,LMA_jack_val_pred_pressed_1300,0.025,0.975),
                                 bias=mean(LMA_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$LMA,na.rm=T))
  
  LDMC_jack_val_pred_pressed_1300<-as.vector(predict(LDMC_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_LDMC_pressed_1300)[,,1])
  LDMC_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$LDMC~LDMC_jack_val_pred_pressed_1300)
  LDMC_jack_stats_pressed_1300[[i]]<-c(R2=summary(LDMC_jack_val_fit_pressed_1300)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed_1300)$LDMC,LDMC_jack_val_pred_pressed_1300),
                                  perRMSE=percentRMSD(meta(val_jack_pressed_1300)$LDMC,LDMC_jack_val_pred_pressed_1300,0.025,0.975),
                                  bias=mean(LDMC_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$LDMC,na.rm=T))
  
  EWT_jack_val_pred_pressed_1300<-as.vector(predict(EWT_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_EWT_pressed_1300)[,,1])
  EWT_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$EWT~EWT_jack_val_pred_pressed_1300)
  EWT_jack_stats_pressed_1300[[i]]<-c(R2=summary(EWT_jack_val_fit_pressed_1300)$r.squared,
                                 RMSE=RMSD(meta(val_jack_pressed_1300)$EWT,EWT_jack_val_pred_pressed_1300),
                                 perRMSE=percentRMSD(meta(val_jack_pressed_1300)$EWT,EWT_jack_val_pred_pressed_1300,0.025,0.975),
                                 bias=mean(EWT_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$EWT,na.rm=T))
  
  chlA_jack_val_pred_pressed_1300<-as.vector(predict(chlA_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_chlA_pressed_1300)[,,1])
  chlA_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$chlA~chlA_jack_val_pred_pressed_1300)
  chlA_jack_stats_pressed_1300[[i]]<-c(R2=summary(chlA_jack_val_fit_pressed_1300)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed_1300)$chlA,chlA_jack_val_pred_pressed_1300),
                                  perRMSE=percentRMSD(meta(val_jack_pressed_1300)$chlA,chlA_jack_val_pred_pressed_1300,0.025,0.975),
                                  bias=mean(chlA_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$chlA,na.rm=T))
  
  chlB_jack_val_pred_pressed_1300<-as.vector(predict(chlB_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_chlB_pressed_1300)[,,1])
  chlB_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$chlB~chlB_jack_val_pred_pressed_1300)
  chlB_jack_stats_pressed_1300[[i]]<-c(R2=summary(chlB_jack_val_fit_pressed_1300)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed_1300)$chlB,chlB_jack_val_pred_pressed_1300),
                                  perRMSE=percentRMSD(meta(val_jack_pressed_1300)$chlB,chlB_jack_val_pred_pressed_1300,0.025,0.975),
                                  bias=mean(chlB_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$chlB,na.rm=T))
  
  car_jack_val_pred_pressed_1300<-as.vector(predict(car_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_car_pressed_1300)[,,1])
  car_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$car~car_jack_val_pred_pressed_1300)
  car_jack_stats_pressed_1300[[i]]<-c(R2=summary(car_jack_val_fit_pressed_1300)$r.squared,
                                 RMSE=RMSD(meta(val_jack_pressed_1300)$car,car_jack_val_pred_pressed_1300),
                                 perRMSE=percentRMSD(meta(val_jack_pressed_1300)$car,car_jack_val_pred_pressed_1300,0.025,0.975),
                                 bias=mean(car_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$car,na.rm=T))
  
  Al_jack_val_pred_pressed_1300<-as.vector(predict(Al_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Al_pressed_1300)[,,1])
  Al_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Al~Al_jack_val_pred_pressed_1300)
  Al_jack_stats_pressed_1300[[i]]<-c(R2=summary(Al_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Al,Al_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Al,Al_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Al_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Al,na.rm=T))
  
  Ca_jack_val_pred_pressed_1300<-as.vector(predict(Ca_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Ca_pressed_1300)[,,1])
  Ca_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Ca~Ca_jack_val_pred_pressed_1300)
  Ca_jack_stats_pressed_1300[[i]]<-c(R2=summary(Ca_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Ca,Ca_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Ca,Ca_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Ca_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Ca,na.rm=T))
  
  Cu_jack_val_pred_pressed_1300<-as.vector(predict(Cu_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Cu_pressed_1300)[,,1])
  Cu_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Cu~Cu_jack_val_pred_pressed_1300)
  Cu_jack_stats_pressed_1300[[i]]<-c(R2=summary(Cu_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Cu,Cu_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Cu,Cu_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Cu_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Cu,na.rm=T))
  
  Fe_jack_val_pred_pressed_1300<-as.vector(predict(Fe_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Fe_pressed_1300)[,,1])
  Fe_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Fe~Fe_jack_val_pred_pressed_1300)
  Fe_jack_stats_pressed_1300[[i]]<-c(R2=summary(Fe_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Fe,Fe_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Fe,Fe_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Fe_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Fe,na.rm=T))
  
  K_jack_val_pred_pressed_1300<-as.vector(predict(K_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_K_pressed_1300)[,,1])
  K_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$K~K_jack_val_pred_pressed_1300)
  K_jack_stats_pressed_1300[[i]]<-c(R2=summary(K_jack_val_fit_pressed_1300)$r.squared,
                               RMSE=RMSD(meta(val_jack_pressed_1300)$K,K_jack_val_pred_pressed_1300),
                               perRMSE=percentRMSD(meta(val_jack_pressed_1300)$K,K_jack_val_pred_pressed_1300,0.025,0.975),
                               bias=mean(K_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$K,na.rm=T))
  
  Mg_jack_val_pred_pressed_1300<-as.vector(predict(Mg_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Mg_pressed_1300)[,,1])
  Mg_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Mg~Mg_jack_val_pred_pressed_1300)
  Mg_jack_stats_pressed_1300[[i]]<-c(R2=summary(Mg_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Mg,Mg_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Mg,Mg_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Mg_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Mg,na.rm=T))
  
  Mn_jack_val_pred_pressed_1300<-as.vector(predict(Mn_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Mn_pressed_1300)[,,1])
  Mn_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Mn~Mn_jack_val_pred_pressed_1300)
  Mn_jack_stats_pressed_1300[[i]]<-c(R2=summary(Mn_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Mn,Mn_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Mn,Mn_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Mn_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Mn,na.rm=T))
  
  Na_jack_val_pred_pressed_1300<-as.vector(predict(Na_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Na_pressed_1300)[,,1])
  Na_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Na~Na_jack_val_pred_pressed_1300)
  Na_jack_stats_pressed_1300[[i]]<-c(R2=summary(Na_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Na,Na_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Na,Na_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Na_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Na,na.rm=T))
  
  P_jack_val_pred_pressed_1300<-as.vector(predict(P_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_P_pressed_1300)[,,1])
  P_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$P~P_jack_val_pred_pressed_1300)
  P_jack_stats_pressed_1300[[i]]<-c(R2=summary(P_jack_val_fit_pressed_1300)$r.squared,
                               RMSE=RMSD(meta(val_jack_pressed_1300)$P,P_jack_val_pred_pressed_1300),
                               perRMSE=percentRMSD(meta(val_jack_pressed_1300)$P,P_jack_val_pred_pressed_1300,0.025,0.975),
                               bias=mean(P_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$P,na.rm=T))
  
  Zn_jack_val_pred_pressed_1300<-as.vector(predict(Zn_pressed_1300_jack,newdata=as.matrix(val_jack_pressed_1300),ncomp=ncomp_Zn_pressed_1300)[,,1])
  Zn_jack_val_fit_pressed_1300<-lm(meta(val_jack_pressed_1300)$Zn~Zn_jack_val_pred_pressed_1300)
  Zn_jack_stats_pressed_1300[[i]]<-c(R2=summary(Zn_jack_val_fit_pressed_1300)$r.squared,
                                RMSE=RMSD(meta(val_jack_pressed_1300)$Zn,Zn_jack_val_pred_pressed_1300),
                                perRMSE=percentRMSD(meta(val_jack_pressed_1300)$Zn,Zn_jack_val_pred_pressed_1300,0.025,0.975),
                                bias=mean(Zn_jack_val_pred_pressed_1300,na.rm=T)-mean(meta(val_jack_pressed_1300)$Zn,na.rm=T))
  
  solubles_jack_coefs_pressed_1300[[i]]<-as.vector(coef(solubles_pressed_1300_jack,ncomp=ncomp_solubles_pressed_1300,intercept=TRUE))
  hemicellulose_jack_coefs_pressed_1300[[i]]<-as.vector(coef(hemicellulose_pressed_1300_jack,ncomp=ncomp_hemicellulose_pressed_1300,intercept=TRUE))
  cellulose_jack_coefs_pressed_1300[[i]]<-as.vector(coef(cellulose_pressed_1300_jack,ncomp=ncomp_cellulose_pressed_1300,intercept=TRUE))
  lignin_jack_coefs_pressed_1300[[i]]<-as.vector(coef(lignin_pressed_1300_jack,ncomp=ncomp_lignin_pressed_1300,intercept=TRUE))
  perC_jack_coefs_pressed_1300[[i]]<-as.vector(coef(perC_pressed_1300_jack,ncomp=ncomp_perC_pressed_1300,intercept=TRUE))
  perN_jack_coefs_pressed_1300[[i]]<-as.vector(coef(perN_pressed_1300_jack,ncomp=ncomp_perN_pressed_1300,intercept=TRUE))
  LMA_jack_coefs_pressed_1300[[i]]<-as.vector(coef(LMA_pressed_1300_jack,ncomp=ncomp_LMA_pressed_1300,intercept=TRUE))
  LDMC_jack_coefs_pressed_1300[[i]]<-as.vector(coef(LDMC_pressed_1300_jack,ncomp=ncomp_LDMC_pressed_1300,intercept=TRUE))
  EWT_jack_coefs_pressed_1300[[i]]<-as.vector(coef(EWT_pressed_1300_jack,ncomp=ncomp_EWT_pressed_1300,intercept=TRUE))
  chlA_jack_coefs_pressed_1300[[i]]<-as.vector(coef(chlA_pressed_1300_jack,ncomp=ncomp_chlA_pressed_1300,intercept=TRUE))
  chlB_jack_coefs_pressed_1300[[i]]<-as.vector(coef(chlB_pressed_1300_jack,ncomp=ncomp_chlB_pressed_1300,intercept=TRUE))
  car_jack_coefs_pressed_1300[[i]]<-as.vector(coef(car_pressed_1300_jack,ncomp=ncomp_car_pressed_1300,intercept=TRUE))
  Al_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Al_pressed_1300_jack,ncomp=ncomp_Al_pressed_1300,intercept=TRUE))
  Ca_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Ca_pressed_1300_jack,ncomp=ncomp_Ca_pressed_1300,intercept=TRUE))
  Cu_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Cu_pressed_1300_jack,ncomp=ncomp_Cu_pressed_1300,intercept=TRUE))
  Fe_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Fe_pressed_1300_jack,ncomp=ncomp_Fe_pressed_1300,intercept=TRUE))
  K_jack_coefs_pressed_1300[[i]]<-as.vector(coef(K_pressed_1300_jack,ncomp=ncomp_K_pressed_1300,intercept=TRUE))
  Mg_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Mg_pressed_1300_jack,ncomp=ncomp_Mg_pressed_1300,intercept=TRUE))
  Mn_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Mn_pressed_1300_jack,ncomp=ncomp_Mn_pressed_1300,intercept=TRUE))
  Na_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Na_pressed_1300_jack,ncomp=ncomp_Na_pressed_1300,intercept=TRUE))
  P_jack_coefs_pressed_1300[[i]]<-as.vector(coef(P_pressed_1300_jack,ncomp=ncomp_P_pressed_1300,intercept=TRUE))
  Zn_jack_coefs_pressed_1300[[i]]<-as.vector(coef(Zn_pressed_1300_jack,ncomp=ncomp_Zn_pressed_1300,intercept=TRUE))
}

solubles_jack_pred_pressed_1300<-apply.coefs(solubles_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
solubles_jack_stat_pressed_1300<-t(apply(solubles_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_jack_df_pressed_1300<-data.frame(pred_mean=solubles_jack_stat_pressed_1300[,1],
                                     pred_low=solubles_jack_stat_pressed_1300[,2],
                                     pred_high=solubles_jack_stat_pressed_1300[,3],
                                     Measured=meta(pressed_1300_spec_all_test)$solubles,
                                     ncomp=ncomp_solubles_pressed_1300,
                                     Project=meta(pressed_1300_spec_all_test)$Project,
                                     Species=meta(pressed_1300_spec_all_test)$Species,
                                     GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                     Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                     NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                     deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                     ID=meta(pressed_1300_spec_all_test)$ID)

hemicellulose_jack_pred_pressed_1300<-apply.coefs(hemicellulose_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
hemicellulose_jack_stat_pressed_1300<-t(apply(hemicellulose_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_jack_df_pressed_1300<-data.frame(pred_mean=hemicellulose_jack_stat_pressed_1300[,1],
                                          pred_low=hemicellulose_jack_stat_pressed_1300[,2],
                                          pred_high=hemicellulose_jack_stat_pressed_1300[,3],
                                          Measured=meta(pressed_1300_spec_all_test)$hemicellulose,
                                          ncomp=ncomp_hemicellulose_pressed_1300,
                                          Project=meta(pressed_1300_spec_all_test)$Project,
                                          Species=meta(pressed_1300_spec_all_test)$Species,
                                          GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                          Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                          NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                          deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                          ID=meta(pressed_1300_spec_all_test)$ID)

cellulose_jack_pred_pressed_1300<-apply.coefs(cellulose_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
cellulose_jack_stat_pressed_1300<-t(apply(cellulose_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_jack_df_pressed_1300<-data.frame(pred_mean=cellulose_jack_stat_pressed_1300[,1],
                                      pred_low=cellulose_jack_stat_pressed_1300[,2],
                                      pred_high=cellulose_jack_stat_pressed_1300[,3],
                                      Measured=meta(pressed_1300_spec_all_test)$cellulose,
                                      ncomp=ncomp_cellulose_pressed_1300,
                                      Project=meta(pressed_1300_spec_all_test)$Project,
                                      Species=meta(pressed_1300_spec_all_test)$Species,
                                      GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                      Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                      NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                      deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                      ID=meta(pressed_1300_spec_all_test)$ID)

lignin_jack_pred_pressed_1300<-apply.coefs(lignin_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
lignin_jack_stat_pressed_1300<-t(apply(lignin_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_jack_df_pressed_1300<-data.frame(pred_mean=lignin_jack_stat_pressed_1300[,1],
                                   pred_low=lignin_jack_stat_pressed_1300[,2],
                                   pred_high=lignin_jack_stat_pressed_1300[,3],
                                   Measured=meta(pressed_1300_spec_all_test)$lignin,
                                   ncomp=ncomp_lignin_pressed_1300,
                                   Project=meta(pressed_1300_spec_all_test)$Project,
                                   Species=meta(pressed_1300_spec_all_test)$Species,
                                   GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                   Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                   NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                   deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                   ID=meta(pressed_1300_spec_all_test)$ID)

perC_jack_pred_pressed_1300<-apply.coefs(perC_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
perC_jack_stat_pressed_1300<-t(apply(perC_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_jack_df_pressed_1300<-data.frame(pred_mean=perC_jack_stat_pressed_1300[,1],
                                 pred_low=perC_jack_stat_pressed_1300[,2],
                                 pred_high=perC_jack_stat_pressed_1300[,3],
                                 Measured=meta(pressed_1300_spec_all_test)$C,
                                 ncomp=ncomp_perC_pressed_1300,
                                 Project=meta(pressed_1300_spec_all_test)$Project,
                                 Species=meta(pressed_1300_spec_all_test)$Species,
                                 GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                 NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                 deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                 ID=meta(pressed_1300_spec_all_test)$ID)

perN_jack_pred_pressed_1300<-apply.coefs(perN_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
perN_jack_stat_pressed_1300<-t(apply(perN_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_jack_df_pressed_1300<-data.frame(pred_mean=perN_jack_stat_pressed_1300[,1],
                                 pred_low=perN_jack_stat_pressed_1300[,2],
                                 pred_high=perN_jack_stat_pressed_1300[,3],
                                 Measured=meta(pressed_1300_spec_all_test)$N,
                                 ncomp=ncomp_perN_pressed_1300,
                                 Project=meta(pressed_1300_spec_all_test)$Project,
                                 Species=meta(pressed_1300_spec_all_test)$Species,
                                 GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                 NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                 deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                 ID=meta(pressed_1300_spec_all_test)$ID)

perN_pressed_1300_val_plot<-ggplot(perN_jack_df_pressed_1300,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured N"[mass]*" (%)"),x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

LMA_jack_pred_pressed_1300<-apply.coefs(LMA_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
LMA_jack_stat_pressed_1300<-t(apply(LMA_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df_pressed_1300<-data.frame(pred_mean=LMA_jack_stat_pressed_1300[,1],
                                pred_low=LMA_jack_stat_pressed_1300[,2],
                                pred_high=LMA_jack_stat_pressed_1300[,3],
                                Measured=meta(pressed_1300_spec_all_test)$LMA,
                                ncomp=ncomp_LMA_pressed_1300,
                                Project=meta(pressed_1300_spec_all_test)$Project,
                                Species=meta(pressed_1300_spec_all_test)$Species,
                                GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                ID=meta(pressed_1300_spec_all_test)$ID)

LDMC_jack_pred_pressed_1300<-apply.coefs(LDMC_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
LDMC_jack_stat_pressed_1300<-t(apply(LDMC_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC_jack_df_pressed_1300<-data.frame(pred_mean=LDMC_jack_stat_pressed_1300[,1],
                                 pred_low=LDMC_jack_stat_pressed_1300[,2],
                                 pred_high=LDMC_jack_stat_pressed_1300[,3],
                                 Measured=meta(pressed_1300_spec_all_test)$LDMC,
                                 ncomp=ncomp_LDMC_pressed_1300,
                                 Project=meta(pressed_1300_spec_all_test)$Project,
                                 Species=meta(pressed_1300_spec_all_test)$Species,
                                 GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                 NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                 deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                 ID=meta(pressed_1300_spec_all_test)$ID)

EWT_jack_pred_pressed_1300<-apply.coefs(EWT_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
EWT_jack_stat_pressed_1300<-t(apply(EWT_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_jack_df_pressed_1300<-data.frame(pred_mean=EWT_jack_stat_pressed_1300[,1],
                                pred_low=EWT_jack_stat_pressed_1300[,2],
                                pred_high=EWT_jack_stat_pressed_1300[,3],
                                Measured=meta(pressed_1300_spec_all_test)$EWT,
                                ncomp=ncomp_EWT_pressed_1300,
                                Project=meta(pressed_1300_spec_all_test)$Project,
                                Species=meta(pressed_1300_spec_all_test)$Species,
                                GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                ID=meta(pressed_1300_spec_all_test)$ID)

chlA_jack_pred_pressed_1300<-apply.coefs(chlA_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
chlA_jack_stat_pressed_1300<-t(apply(chlA_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_jack_df_pressed_1300<-data.frame(pred_mean=chlA_jack_stat_pressed_1300[,1],
                                 pred_low=chlA_jack_stat_pressed_1300[,2],
                                 pred_high=chlA_jack_stat_pressed_1300[,3],
                                 Measured=meta(pressed_1300_spec_all_test)$chlA,
                                 ncomp=ncomp_chlA_pressed_1300,
                                 Project=meta(pressed_1300_spec_all_test)$Project,
                                 Species=meta(pressed_1300_spec_all_test)$Species,
                                 GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                 NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                 deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                 ID=meta(pressed_1300_spec_all_test)$ID)

chlB_jack_pred_pressed_1300<-apply.coefs(chlB_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
chlB_jack_stat_pressed_1300<-t(apply(chlB_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_jack_df_pressed_1300<-data.frame(pred_mean=chlB_jack_stat_pressed_1300[,1],
                                 pred_low=chlB_jack_stat_pressed_1300[,2],
                                 pred_high=chlB_jack_stat_pressed_1300[,3],
                                 Measured=meta(pressed_1300_spec_all_test)$chlB,
                                 ncomp=ncomp_chlB_pressed_1300,
                                 Project=meta(pressed_1300_spec_all_test)$Project,
                                 Species=meta(pressed_1300_spec_all_test)$Species,
                                 GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                 NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                 deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                 ID=meta(pressed_1300_spec_all_test)$ID)

car_jack_pred_pressed_1300<-apply.coefs(car_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
car_jack_stat_pressed_1300<-t(apply(car_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_jack_df_pressed_1300<-data.frame(pred_mean=car_jack_stat_pressed_1300[,1],
                                pred_low=car_jack_stat_pressed_1300[,2],
                                pred_high=car_jack_stat_pressed_1300[,3],
                                Measured=meta(pressed_1300_spec_all_test)$car,
                                ncomp=ncomp_car_pressed_1300,
                                Project=meta(pressed_1300_spec_all_test)$Project,
                                Species=meta(pressed_1300_spec_all_test)$Species,
                                GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                                Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                                NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                                deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                                ID=meta(pressed_1300_spec_all_test)$ID)

Al_jack_pred_pressed_1300<-apply.coefs(Al_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Al_jack_stat_pressed_1300<-t(apply(Al_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_jack_df_pressed_1300<-data.frame(pred_mean=Al_jack_stat_pressed_1300[,1],
                               pred_low=Al_jack_stat_pressed_1300[,2],
                               pred_high=Al_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Al,
                               ncomp=ncomp_Al_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

Ca_jack_pred_pressed_1300<-apply.coefs(Ca_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Ca_jack_stat_pressed_1300<-t(apply(Ca_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_jack_df_pressed_1300<-data.frame(pred_mean=Ca_jack_stat_pressed_1300[,1],
                               pred_low=Ca_jack_stat_pressed_1300[,2],
                               pred_high=Ca_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Ca,
                               ncomp=ncomp_Ca_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

Cu_jack_pred_pressed_1300<-apply.coefs(Cu_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Cu_jack_stat_pressed_1300<-t(apply(Cu_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_jack_df_pressed_1300<-data.frame(pred_mean=Cu_jack_stat_pressed_1300[,1],
                               pred_low=Cu_jack_stat_pressed_1300[,2],
                               pred_high=Cu_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Cu,
                               ncomp=ncomp_Cu_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

Fe_jack_pred_pressed_1300<-apply.coefs(Fe_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Fe_jack_stat_pressed_1300<-t(apply(Fe_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_jack_df_pressed_1300<-data.frame(pred_mean=Fe_jack_stat_pressed_1300[,1],
                               pred_low=Fe_jack_stat_pressed_1300[,2],
                               pred_high=Fe_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Fe,
                               ncomp=ncomp_Fe_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

K_jack_pred_pressed_1300<-apply.coefs(K_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
K_jack_stat_pressed_1300<-t(apply(K_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_jack_df_pressed_1300<-data.frame(pred_mean=K_jack_stat_pressed_1300[,1],
                              pred_low=K_jack_stat_pressed_1300[,2],
                              pred_high=K_jack_stat_pressed_1300[,3],
                              Measured=meta(pressed_1300_spec_all_test)$K,
                              ncomp=ncomp_K_pressed_1300,
                              Project=meta(pressed_1300_spec_all_test)$Project,
                              Species=meta(pressed_1300_spec_all_test)$Species,
                              GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                              Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                              NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                              deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                              ID=meta(pressed_1300_spec_all_test)$ID)

Mg_jack_pred_pressed_1300<-apply.coefs(Mg_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Mg_jack_stat_pressed_1300<-t(apply(Mg_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_jack_df_pressed_1300<-data.frame(pred_mean=Mg_jack_stat_pressed_1300[,1],
                               pred_low=Mg_jack_stat_pressed_1300[,2],
                               pred_high=Mg_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Mg,
                               ncomp=ncomp_Mg_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

Mn_jack_pred_pressed_1300<-apply.coefs(Mn_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Mn_jack_stat_pressed_1300<-t(apply(Mn_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_jack_df_pressed_1300<-data.frame(pred_mean=Mn_jack_stat_pressed_1300[,1],
                               pred_low=Mn_jack_stat_pressed_1300[,2],
                               pred_high=Mn_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Mn,
                               ncomp=ncomp_Mn_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

Na_jack_pred_pressed_1300<-apply.coefs(Na_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Na_jack_stat_pressed_1300<-t(apply(Na_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_jack_df_pressed_1300<-data.frame(pred_mean=Na_jack_stat_pressed_1300[,1],
                               pred_low=Na_jack_stat_pressed_1300[,2],
                               pred_high=Na_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Na,
                               ncomp=ncomp_Na_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

P_jack_pred_pressed_1300<-apply.coefs(P_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
P_jack_stat_pressed_1300<-t(apply(P_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_jack_df_pressed_1300<-data.frame(pred_mean=P_jack_stat_pressed_1300[,1],
                              pred_low=P_jack_stat_pressed_1300[,2],
                              pred_high=P_jack_stat_pressed_1300[,3],
                              Measured=meta(pressed_1300_spec_all_test)$P,
                              ncomp=ncomp_P_pressed_1300,
                              Project=meta(pressed_1300_spec_all_test)$Project,
                              Species=meta(pressed_1300_spec_all_test)$Species,
                              GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                              Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                              NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                              deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                              ID=meta(pressed_1300_spec_all_test)$ID)

Zn_jack_pred_pressed_1300<-apply.coefs(Zn_jack_coefs_pressed_1300,as.matrix(pressed_1300_spec_all_test))
Zn_jack_stat_pressed_1300<-t(apply(Zn_jack_pred_pressed_1300,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_jack_df_pressed_1300<-data.frame(pred_mean=Zn_jack_stat_pressed_1300[,1],
                               pred_low=Zn_jack_stat_pressed_1300[,2],
                               pred_high=Zn_jack_stat_pressed_1300[,3],
                               Measured=meta(pressed_1300_spec_all_test)$Zn,
                               ncomp=ncomp_Zn_pressed_1300,
                               Project=meta(pressed_1300_spec_all_test)$Project,
                               Species=meta(pressed_1300_spec_all_test)$Species,
                               GrowthForm=meta(pressed_1300_spec_all_test)$GrowthForm,
                               Discoloration=meta(pressed_1300_spec_all_test)$Discoloration,
                               NDVI=meta(pressed_1300_spec_all_test)$NDVI,
                               deg_ind=meta(pressed_1300_spec_all_test)$deg_ind,
                               ID=meta(pressed_1300_spec_all_test)$ID)

##########################################################
## save output

pressed_1300_jack_coef_list<-list(LMA=LMA_jack_coefs_pressed_1300,
                                  LDMC=LDMC_jack_coefs_pressed_1300,
                                  EWT=EWT_jack_coefs_pressed_1300,
                                  sol=solubles_jack_coefs_pressed_1300,
                                  hemi=hemicellulose_jack_coefs_pressed_1300,
                                  cell=cellulose_jack_coefs_pressed_1300,
                                  lign=lignin_jack_coefs_pressed_1300,
                                  C=perC_jack_coefs_pressed_1300,
                                  N=perN_jack_coefs_pressed_1300,
                                  chlA=chlA_jack_coefs_pressed_1300,
                                  chlB=chlB_jack_coefs_pressed_1300,
                                  car=car_jack_coefs_pressed_1300,
                                  Al=Al_jack_coefs_pressed_1300,
                                  Ca=Ca_jack_coefs_pressed_1300,
                                  Cu=Cu_jack_coefs_pressed_1300,
                                  Fe=Fe_jack_coefs_pressed_1300,
                                  K=K_jack_coefs_pressed_1300,
                                  Mg=Mg_jack_coefs_pressed_1300,
                                  Mn=Mn_jack_coefs_pressed_1300,
                                  Na=Na_jack_coefs_pressed_1300,
                                  P=P_jack_coefs_pressed_1300,
                                  Zn=Zn_jack_coefs_pressed_1300)
saveRDS(pressed_1300_jack_coef_list,"SavedResults/pressed_1300_jack_coefs_list.rds")

pressed_1300_jack_df_list<-list(LMA=LMA_jack_df_pressed_1300,
                                LDMC=LDMC_jack_df_pressed_1300,
                                EWT=EWT_jack_df_pressed_1300,
                                sol=solubles_jack_df_pressed_1300,
                                hemi=hemicellulose_jack_df_pressed_1300,
                                cell=cellulose_jack_df_pressed_1300,
                                lign=lignin_jack_df_pressed_1300,
                                chlA=chlA_jack_df_pressed_1300,
                                chlB=chlB_jack_df_pressed_1300,
                                car=car_jack_df_pressed_1300,
                                C=perC_jack_df_pressed_1300,
                                N=perN_jack_df_pressed_1300,
                                Al=Al_jack_df_pressed_1300,
                                Ca=Ca_jack_df_pressed_1300,
                                Cu=Cu_jack_df_pressed_1300,
                                Fe=Fe_jack_df_pressed_1300,
                                K=K_jack_df_pressed_1300,
                                Mg=Mg_jack_df_pressed_1300,
                                Mn=Mn_jack_df_pressed_1300,
                                Na=Na_jack_df_pressed_1300,
                                P=P_jack_df_pressed_1300,
                                Zn=Zn_jack_df_pressed_1300)
saveRDS(pressed_1300_jack_df_list,"SavedResults/pressed_1300_jack_df_list.rds")

############################################
## violin plots

R2.df_pressed_1300<-data.frame(LMA=unlist(lapply(LMA_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               LDMC=unlist(lapply(LDMC_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               EWT=unlist(lapply(EWT_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               sol=unlist(lapply(solubles_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               hemi=unlist(lapply(hemicellulose_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               cell=unlist(lapply(cellulose_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               lign=unlist(lapply(lignin_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               C=unlist(lapply(perC_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               N=unlist(lapply(perN_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               LMA=unlist(lapply(LMA_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               LDMC=unlist(lapply(LDMC_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               EWT=unlist(lapply(EWT_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               chlA=unlist(lapply(chlA_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               chlB=unlist(lapply(chlB_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               car=unlist(lapply(car_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Al=unlist(lapply(Al_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Ca=unlist(lapply(Ca_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Cu=unlist(lapply(Cu_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Fe=unlist(lapply(Fe_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               K=unlist(lapply(K_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Mg=unlist(lapply(Mg_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Mn=unlist(lapply(Mn_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Na=unlist(lapply(Na_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               P=unlist(lapply(P_jack_stats_pressed_1300,function(x) x[["R2"]])),
                               Zn=unlist(lapply(Zn_jack_stats_pressed_1300,function(x) x[["R2"]])))

R2.long_pressed_1300<-melt(R2.df_pressed_1300)
pressed_1300_val_R2<-ggplot(R2.long_pressed_1300,aes(y=value,x=variable))+
  geom_violin()+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y=expression(italic("R"^2)))+
  ggtitle("Pressed-leaf spectra: 1300-2500 nm")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))

perRMSE.df_pressed_1300<-data.frame(LMA=unlist(lapply(LMA_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    LDMC=unlist(lapply(LDMC_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    EWT=unlist(lapply(EWT_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    sol=unlist(lapply(solubles_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    hemi=unlist(lapply(hemicellulose_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    cell=unlist(lapply(cellulose_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    lign=unlist(lapply(lignin_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    C=unlist(lapply(perC_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    N=unlist(lapply(perN_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    chlA=unlist(lapply(chlA_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    chlB=unlist(lapply(chlB_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    car=unlist(lapply(car_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Al=unlist(lapply(Al_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Ca=unlist(lapply(Ca_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Cu=unlist(lapply(Cu_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Fe=unlist(lapply(Fe_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    K=unlist(lapply(K_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Mg=unlist(lapply(Mg_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Mn=unlist(lapply(Mn_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Na=unlist(lapply(Na_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    P=unlist(lapply(P_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])),
                                    Zn=unlist(lapply(Zn_jack_stats_pressed_1300,function(x) 100*x[["perRMSE"]])))

perRMSE.long_pressed_1300<-melt(perRMSE.df_pressed_1300)
pressed_1300_val_perRMSE<-ggplot(perRMSE.long_pressed_1300,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  labs(y="%RMSE")
scale_y_continuous(expand = c(0, 0),
                   limits = c(0,max(perRMSE.long_pressed_1300$value)*1.1))

pdf("Manuscript/FigS12_1300.pdf",width=8,height=6,onefile=F)
egg::ggarrange(plots = list(pressed_1300_val_R2,pressed_1300_val_perRMSE),
               nrow=2,ncol=1)
dev.off()
