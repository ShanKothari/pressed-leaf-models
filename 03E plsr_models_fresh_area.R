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

fresh_spec_all_train<-readRDS("ProcessedSpectralData/fresh_spec_all_train.rds")
fresh_spec_all_test<-readRDS("ProcessedSpectralData/fresh_spec_all_test.rds")

################################################
## building calibration models

perC_area_fresh<-plsr(meta(fresh_spec_all_train)$C_area~as.matrix(fresh_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perC_area_fresh <- selectNcomp(perC_area_fresh, method = "onesigma", plot = FALSE)
perC_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$C_area))
perC_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[perC_area_fresh_valid],
                                 Species=meta(fresh_spec_all_train)$Species[perC_area_fresh_valid],
                                 Project=meta(fresh_spec_all_train)$Project[perC_area_fresh_valid],
                                 Stage=meta(fresh_spec_all_train)$Stage[perC_area_fresh_valid],
                                 GrowthForm=meta(fresh_spec_all_train)$GrowthForm[perC_area_fresh_valid],
                                 measured=meta(fresh_spec_all_train)$C_area[perC_area_fresh_valid],
                                 val_pred=perC_area_fresh$validation$pred[,,ncomp_perC_area_fresh])
ggplot(perC_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %C from fresh-leaf spectra")

perN_area_fresh<-plsr(meta(fresh_spec_all_train)$N_area~as.matrix(fresh_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perN_area_fresh <- selectNcomp(perN_area_fresh, method = "onesigma", plot = FALSE)
perN_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$N_area))
perN_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[perN_area_fresh_valid],
                                 Species=meta(fresh_spec_all_train)$Species[perN_area_fresh_valid],
                                 Project=meta(fresh_spec_all_train)$Project[perN_area_fresh_valid],
                                 Stage=meta(fresh_spec_all_train)$Stage[perN_area_fresh_valid],
                                 GrowthForm=meta(fresh_spec_all_train)$GrowthForm[perN_area_fresh_valid],
                                 measured=meta(fresh_spec_all_train)$N_area[perN_area_fresh_valid],
                                 val_pred=perN_area_fresh$validation$pred[,,ncomp_perN_area_fresh])
ggplot(perN_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %N from fresh-leaf spectra")+guides(color=F)

chlA_area_fresh<-plsr(meta(fresh_spec_all_train)$chlA_area~as.matrix(fresh_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_area_fresh <- selectNcomp(chlA_area_fresh, method = "onesigma", plot = FALSE)
chlA_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$chlA_area))
chlA_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[chlA_area_fresh_valid],
                                 Species=meta(fresh_spec_all_train)$Species[chlA_area_fresh_valid],
                                 Project=meta(fresh_spec_all_train)$Project[chlA_area_fresh_valid],
                                 Stage=meta(fresh_spec_all_train)$Stage[chlA_area_fresh_valid],
                                 GrowthForm=meta(fresh_spec_all_train)$GrowthForm[chlA_area_fresh_valid],
                                 measured=meta(fresh_spec_all_train)$chlA_area[chlA_area_fresh_valid],
                                 val_pred=chlA_area_fresh$validation$pred[,,ncomp_chlA_area_fresh])
ggplot(chlA_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl a from fresh-leaf spectra")

chlB_area_fresh<-plsr(meta(fresh_spec_all_train)$chlB_area~as.matrix(fresh_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_area_fresh <- selectNcomp(chlB_area_fresh, method = "onesigma", plot = FALSE)
chlB_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$chlB_area))
chlB_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[chlB_area_fresh_valid],
                                 Species=meta(fresh_spec_all_train)$Species[chlB_area_fresh_valid],
                                 Project=meta(fresh_spec_all_train)$Project[chlB_area_fresh_valid],
                                 Stage=meta(fresh_spec_all_train)$Stage[chlB_area_fresh_valid],
                                 GrowthForm=meta(fresh_spec_all_train)$GrowthForm[chlB_area_fresh_valid],
                                 measured=meta(fresh_spec_all_train)$chlB_area[chlB_area_fresh_valid],
                                 val_pred=chlB_area_fresh$validation$pred[,,ncomp_chlB_area_fresh])
ggplot(chlB_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl b from fresh-leaf spectra")

car_area_fresh<-plsr(meta(fresh_spec_all_train)$car_area~as.matrix(fresh_spec_all_train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_area_fresh <- selectNcomp(car_area_fresh, method = "onesigma", plot = FALSE)
car_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$car_area))
car_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[car_area_fresh_valid],
                                Species=meta(fresh_spec_all_train)$Species[car_area_fresh_valid],
                                Project=meta(fresh_spec_all_train)$Project[car_area_fresh_valid],
                                Stage=meta(fresh_spec_all_train)$Stage[car_area_fresh_valid],
                                GrowthForm=meta(fresh_spec_all_train)$GrowthForm[car_area_fresh_valid],
                                measured=meta(fresh_spec_all_train)$car_area[car_area_fresh_valid],
                                val_pred=car_area_fresh$validation$pred[,,ncomp_car_area_fresh])
ggplot(car_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting total car from fresh-leaf spectra")

solubles_area_fresh<-plsr(meta(fresh_spec_all_train)$solubles_area~as.matrix(fresh_spec_all_train),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_area_fresh <- selectNcomp(solubles_area_fresh, method = "onesigma", plot = FALSE)
solubles_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$solubles_area))
solubles_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[solubles_area_fresh_valid],
                                     Species=meta(fresh_spec_all_train)$Species[solubles_area_fresh_valid],
                                     Project=meta(fresh_spec_all_train)$Project[solubles_area_fresh_valid],
                                     Stage=meta(fresh_spec_all_train)$Stage[solubles_area_fresh_valid],
                                     GrowthForm=meta(fresh_spec_all_train)$GrowthForm[solubles_area_fresh_valid],
                                     measured=meta(fresh_spec_all_train)$solubles_area[solubles_area_fresh_valid],
                                     val_pred=solubles_area_fresh$validation$pred[,,ncomp_solubles_area_fresh])
ggplot(solubles_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting solubles from fresh-leaf spectra")

hemicellulose_area_fresh<-plsr(meta(fresh_spec_all_train)$hemicellulose_area~as.matrix(fresh_spec_all_train),
                               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_area_fresh <- selectNcomp(hemicellulose_area_fresh, method = "onesigma", plot = FALSE)
hemicellulose_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$hemicellulose_area))
hemicellulose_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[hemicellulose_area_fresh_valid],
                                          Species=meta(fresh_spec_all_train)$Species[hemicellulose_area_fresh_valid],
                                          Project=meta(fresh_spec_all_train)$Project[hemicellulose_area_fresh_valid],
                                          Stage=meta(fresh_spec_all_train)$Stage[hemicellulose_area_fresh_valid],
                                          GrowthForm=meta(fresh_spec_all_train)$GrowthForm[hemicellulose_area_fresh_valid],
                                          measured=meta(fresh_spec_all_train)$hemicellulose_area[hemicellulose_area_fresh_valid],
                                          val_pred=hemicellulose_area_fresh$validation$pred[,,ncomp_hemicellulose_area_fresh])
ggplot(hemicellulose_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting hemicellulose from fresh-leaf spectra")

cellulose_area_fresh<-plsr(meta(fresh_spec_all_train)$cellulose_area~as.matrix(fresh_spec_all_train),
                           ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_area_fresh <- selectNcomp(cellulose_area_fresh, method = "onesigma", plot = FALSE)
cellulose_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$cellulose_area))
cellulose_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[cellulose_area_fresh_valid],
                                      Species=meta(fresh_spec_all_train)$Species[cellulose_area_fresh_valid],
                                      Project=meta(fresh_spec_all_train)$Project[cellulose_area_fresh_valid],
                                      Stage=meta(fresh_spec_all_train)$Stage[cellulose_area_fresh_valid],
                                      GrowthForm=meta(fresh_spec_all_train)$GrowthForm[cellulose_area_fresh_valid],
                                      measured=meta(fresh_spec_all_train)$cellulose_area[cellulose_area_fresh_valid],
                                      val_pred=cellulose_area_fresh$validation$pred[,,ncomp_cellulose_area_fresh])
ggplot(cellulose_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting cellulose from fresh-leaf spectra")

lignin_area_fresh<-plsr(meta(fresh_spec_all_train)$lignin_area~as.matrix(fresh_spec_all_train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_area_fresh <- selectNcomp(lignin_area_fresh, method = "onesigma", plot = FALSE)
lignin_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$lignin_area))
lignin_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[lignin_area_fresh_valid],
                                   Species=meta(fresh_spec_all_train)$Species[lignin_area_fresh_valid],
                                   Project=meta(fresh_spec_all_train)$Project[lignin_area_fresh_valid],
                                   Stage=meta(fresh_spec_all_train)$Stage[lignin_area_fresh_valid],
                                   GrowthForm=meta(fresh_spec_all_train)$GrowthForm[lignin_area_fresh_valid],
                                   measured=meta(fresh_spec_all_train)$lignin_area[lignin_area_fresh_valid],
                                   val_pred=lignin_area_fresh$validation$pred[,,ncomp_lignin_area_fresh])
ggplot(lignin_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting lignin from fresh-leaf spectra")

Al_area_fresh<-plsr(meta(fresh_spec_all_train)$Al_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_area_fresh <- selectNcomp(Al_area_fresh, method = "onesigma", plot = FALSE)
Al_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Al_area))
Al_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Al_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Al_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Al_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Al_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Al_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Al_area[Al_area_fresh_valid],
                               val_pred=Al_area_fresh$validation$pred[,,ncomp_Al_area_fresh])
ggplot(Al_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Al from fresh-leaf spectra")

Ca_area_fresh<-plsr(meta(fresh_spec_all_train)$Ca_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_area_fresh <- selectNcomp(Ca_area_fresh, method = "onesigma", plot = FALSE)
Ca_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Ca_area))
Ca_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Ca_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Ca_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Ca_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Ca_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Ca_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Ca_area[Ca_area_fresh_valid],
                               val_pred=Ca_area_fresh$validation$pred[,,ncomp_Ca_area_fresh])
ggplot(Ca_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Ca from fresh-leaf spectra")

Cu_area_fresh<-plsr(meta(fresh_spec_all_train)$Cu_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_area_fresh <- selectNcomp(Cu_area_fresh, method = "onesigma", plot = FALSE)
Cu_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Cu_area))
Cu_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Cu_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Cu_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Cu_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Cu_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Cu_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Cu_area[Cu_area_fresh_valid],
                               val_pred=Cu_area_fresh$validation$pred[,,ncomp_Cu_area_fresh])
ggplot(Cu_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Cu from fresh-leaf spectra")

Fe_area_fresh<-plsr(meta(fresh_spec_all_train)$Fe_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_area_fresh <- selectNcomp(Fe_area_fresh, method = "onesigma", plot = FALSE)
Fe_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Fe_area))
Fe_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Fe_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Fe_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Fe_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Fe_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Fe_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Fe_area[Fe_area_fresh_valid],
                               val_pred=Fe_area_fresh$validation$pred[,,ncomp_Fe_area_fresh])
ggplot(Fe_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Fe from fresh-leaf spectra")

K_area_fresh<-plsr(meta(fresh_spec_all_train)$K_area~as.matrix(fresh_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_area_fresh <- selectNcomp(K_area_fresh, method = "onesigma", plot = FALSE)
K_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$K_area))
K_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[K_area_fresh_valid],
                              Species=meta(fresh_spec_all_train)$Species[K_area_fresh_valid],
                              Project=meta(fresh_spec_all_train)$Project[K_area_fresh_valid],
                              Stage=meta(fresh_spec_all_train)$Stage[K_area_fresh_valid],
                              GrowthForm=meta(fresh_spec_all_train)$GrowthForm[K_area_fresh_valid],
                              measured=meta(fresh_spec_all_train)$K_area[K_area_fresh_valid],
                              val_pred=K_area_fresh$validation$pred[,,ncomp_K_area_fresh])
ggplot(K_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting K from fresh-leaf spectra")

Mg_area_fresh<-plsr(meta(fresh_spec_all_train)$Mg_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_area_fresh <- selectNcomp(Mg_area_fresh, method = "onesigma", plot = FALSE)
Mg_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Mg_area))
Mg_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Mg_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Mg_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Mg_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Mg_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Mg_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Mg_area[Mg_area_fresh_valid],
                               val_pred=Mg_area_fresh$validation$pred[,,ncomp_Mg_area_fresh])
ggplot(Mg_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mg from fresh-leaf spectra")

Mn_area_fresh<-plsr(meta(fresh_spec_all_train)$Mn_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_area_fresh <- selectNcomp(Mn_area_fresh, method = "onesigma", plot = FALSE)
Mn_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Mn_area))
Mn_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Mn_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Mn_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Mn_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Mn_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Mn_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Mn_area[Mn_area_fresh_valid],
                               val_pred=Mn_area_fresh$validation$pred[,,ncomp_Mn_area_fresh])
ggplot(Mn_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mn from fresh-leaf spectra")

Na_area_fresh<-plsr(meta(fresh_spec_all_train)$Na_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_area_fresh <- selectNcomp(Na_area_fresh, method = "onesigma", plot = FALSE)
Na_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Na_area))
Na_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Na_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Na_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Na_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Na_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Na_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Na_area[Na_area_fresh_valid],
                               val_pred=Na_area_fresh$validation$pred[,,ncomp_Na_area_fresh])
ggplot(Na_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Na from fresh-leaf spectra")

P_area_fresh<-plsr(meta(fresh_spec_all_train)$P_area~as.matrix(fresh_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_area_fresh <- selectNcomp(P_area_fresh, method = "onesigma", plot = FALSE)
P_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$P_area))
P_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[P_area_fresh_valid],
                              Species=meta(fresh_spec_all_train)$Species[P_area_fresh_valid],
                              Project=meta(fresh_spec_all_train)$Project[P_area_fresh_valid],
                              Stage=meta(fresh_spec_all_train)$Stage[P_area_fresh_valid],
                              GrowthForm=meta(fresh_spec_all_train)$GrowthForm[P_area_fresh_valid],
                              measured=meta(fresh_spec_all_train)$P_area[P_area_fresh_valid],
                              val_pred=P_area_fresh$validation$pred[,,ncomp_P_area_fresh])
ggplot(P_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting P from fresh-leaf spectra")

Zn_area_fresh<-plsr(meta(fresh_spec_all_train)$Zn_area~as.matrix(fresh_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_area_fresh <- selectNcomp(Zn_area_fresh, method = "onesigma", plot = FALSE)
Zn_area_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$Zn_area))
Zn_area_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[Zn_area_fresh_valid],
                               Species=meta(fresh_spec_all_train)$Species[Zn_area_fresh_valid],
                               Project=meta(fresh_spec_all_train)$Project[Zn_area_fresh_valid],
                               Stage=meta(fresh_spec_all_train)$Stage[Zn_area_fresh_valid],
                               GrowthForm=meta(fresh_spec_all_train)$GrowthForm[Zn_area_fresh_valid],
                               measured=meta(fresh_spec_all_train)$Zn_area[Zn_area_fresh_valid],
                               val_pred=Zn_area_fresh$validation$pred[,,ncomp_Zn_area_fresh])
ggplot(Zn_area_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Zn from fresh-leaf spectra")

###############################################
## jackknife tests + prediction of validation data
## fresh leaves

solubles_area_jack_coefs_fresh<-list()
hemicellulose_area_jack_coefs_fresh<-list()
cellulose_area_jack_coefs_fresh<-list()
lignin_area_jack_coefs_fresh<-list()
perC_area_jack_coefs_fresh<-list()
perN_area_jack_coefs_fresh<-list()
chlA_area_jack_coefs_fresh<-list()
chlB_area_jack_coefs_fresh<-list()
car_area_jack_coefs_fresh<-list()
Al_area_jack_coefs_fresh<-list()
Ca_area_jack_coefs_fresh<-list()
Cu_area_jack_coefs_fresh<-list()
Fe_area_jack_coefs_fresh<-list()
K_area_jack_coefs_fresh<-list()
Mg_area_jack_coefs_fresh<-list()
Mn_area_jack_coefs_fresh<-list()
Na_area_jack_coefs_fresh<-list()
P_area_jack_coefs_fresh<-list()
Zn_area_jack_coefs_fresh<-list()

solubles_area_jack_stats_fresh<-list()
hemicellulose_area_jack_stats_fresh<-list()
cellulose_area_jack_stats_fresh<-list()
lignin_area_jack_stats_fresh<-list()
perC_area_jack_stats_fresh<-list()
perN_area_jack_stats_fresh<-list()
chlA_area_jack_stats_fresh<-list()
chlB_area_jack_stats_fresh<-list()
car_area_jack_stats_fresh<-list()
Al_area_jack_stats_fresh<-list()
Ca_area_jack_stats_fresh<-list()
Cu_area_jack_stats_fresh<-list()
Fe_area_jack_stats_fresh<-list()
K_area_jack_stats_fresh<-list()
Mg_area_jack_stats_fresh<-list()
Mn_area_jack_stats_fresh<-list()
Na_area_jack_stats_fresh<-list()
P_area_jack_stats_fresh<-list()
Zn_area_jack_stats_fresh<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_fresh<-nrow(fresh_spec_all_train)
  train_jack_fresh<-sample(1:n_cal_spec_fresh,floor(0.7*n_cal_spec_fresh))
  test_jack_fresh<-setdiff(1:n_cal_spec_fresh,train_jack_fresh)
  
  calib_jack_fresh<-fresh_spec_all_train[train_jack_fresh]
  val_jack_fresh<-fresh_spec_all_train[test_jack_fresh]
  
  solubles_area_fresh_jack<-plsr(meta(calib_jack_fresh)$solubles_area~as.matrix(calib_jack_fresh),
                                 ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_area_fresh_jack<-plsr(meta(calib_jack_fresh)$hemicellulose_area~as.matrix(calib_jack_fresh),
                                      ncomp=30,method = "oscorespls",validation="none")
  cellulose_area_fresh_jack<-plsr(meta(calib_jack_fresh)$cellulose_area~as.matrix(calib_jack_fresh),
                                  ncomp=30,method = "oscorespls",validation="none")
  lignin_area_fresh_jack<-plsr(meta(calib_jack_fresh)$lignin_area~as.matrix(calib_jack_fresh),
                               ncomp=30,method = "oscorespls",validation="none")
  perC_area_fresh_jack<-plsr(meta(calib_jack_fresh)$C_area~as.matrix(calib_jack_fresh),
                             ncomp=30,method = "oscorespls",validation="none")
  perN_area_fresh_jack<-plsr(meta(calib_jack_fresh)$N_area~as.matrix(calib_jack_fresh),
                             ncomp=30,method = "oscorespls",validation="none")
  chlA_area_fresh_jack<-plsr(meta(calib_jack_fresh)$chlA_area~as.matrix(calib_jack_fresh),
                             ncomp=30,method = "oscorespls",validation="none")
  chlB_area_fresh_jack<-plsr(meta(calib_jack_fresh)$chlB_area~as.matrix(calib_jack_fresh),
                             ncomp=30,method = "oscorespls",validation="none")
  car_area_fresh_jack<-plsr(meta(calib_jack_fresh)$car_area~as.matrix(calib_jack_fresh),
                            ncomp=30,method = "oscorespls",validation="none")
  Al_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Al_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  Ca_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Ca_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  Cu_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Cu_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  Fe_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Fe_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  K_area_fresh_jack<-plsr(meta(calib_jack_fresh)$K_area~as.matrix(calib_jack_fresh),
                          ncomp=30,method = "oscorespls",validation="none")
  Mg_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Mg_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  Mn_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Mn_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  Na_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Na_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  P_area_fresh_jack<-plsr(meta(calib_jack_fresh)$P_area~as.matrix(calib_jack_fresh),
                          ncomp=30,method = "oscorespls",validation="none")
  Zn_area_fresh_jack<-plsr(meta(calib_jack_fresh)$Zn_area~as.matrix(calib_jack_fresh),
                           ncomp=30,method = "oscorespls",validation="none")
  
  solubles_area_jack_val_pred_fresh<-as.vector(predict(solubles_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_solubles_area_fresh)[,,1])
  solubles_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$solubles_area~solubles_area_jack_val_pred_fresh)
  solubles_area_jack_stats_fresh[[i]]<-c(R2=summary(solubles_area_jack_val_fit_fresh)$r.squared,
                                         RMSE=RMSD(meta(val_jack_fresh)$solubles_area,solubles_area_jack_val_pred_fresh),
                                         perRMSE=percentRMSD(meta(val_jack_fresh)$solubles_area,solubles_area_jack_val_pred_fresh,0.025,0.975),
                                         bias=mean(solubles_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$solubles_area,na.rm=T))
  
  hemicellulose_area_jack_val_pred_fresh<-as.vector(predict(hemicellulose_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_hemicellulose_area_fresh)[,,1])
  hemicellulose_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$hemicellulose_area~hemicellulose_area_jack_val_pred_fresh)
  hemicellulose_area_jack_stats_fresh[[i]]<-c(R2=summary(hemicellulose_area_jack_val_fit_fresh)$r.squared,
                                              RMSE=RMSD(meta(val_jack_fresh)$hemicellulose_area,hemicellulose_area_jack_val_pred_fresh),
                                              perRMSE=percentRMSD(meta(val_jack_fresh)$hemicellulose_area,hemicellulose_area_jack_val_pred_fresh,0.025,0.975),
                                              bias=mean(hemicellulose_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$hemicellulose_area,na.rm=T))
  
  cellulose_area_jack_val_pred_fresh<-as.vector(predict(cellulose_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_cellulose_area_fresh)[,,1])
  cellulose_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$cellulose_area~cellulose_area_jack_val_pred_fresh)
  cellulose_area_jack_stats_fresh[[i]]<-c(R2=summary(cellulose_area_jack_val_fit_fresh)$r.squared,
                                          RMSE=RMSD(meta(val_jack_fresh)$cellulose_area,cellulose_area_jack_val_pred_fresh),
                                          perRMSE=percentRMSD(meta(val_jack_fresh)$cellulose_area,cellulose_area_jack_val_pred_fresh,0.025,0.975),
                                          bias=mean(cellulose_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$cellulose_area,na.rm=T))
  
  lignin_area_jack_val_pred_fresh<-as.vector(predict(lignin_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_lignin_area_fresh)[,,1])
  lignin_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$lignin_area~lignin_area_jack_val_pred_fresh)
  lignin_area_jack_stats_fresh[[i]]<-c(R2=summary(lignin_area_jack_val_fit_fresh)$r.squared,
                                       RMSE=RMSD(meta(val_jack_fresh)$lignin_area,lignin_area_jack_val_pred_fresh),
                                       perRMSE=percentRMSD(meta(val_jack_fresh)$lignin_area,lignin_area_jack_val_pred_fresh,0.025,0.975),
                                       bias=mean(lignin_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$lignin_area,na.rm=T))
  
  perC_area_jack_val_pred_fresh<-as.vector(predict(perC_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_perC_area_fresh)[,,1])
  perC_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$C_area~perC_area_jack_val_pred_fresh)
  perC_area_jack_stats_fresh[[i]]<-c(R2=summary(perC_area_jack_val_fit_fresh)$r.squared,
                                     RMSE=RMSD(meta(val_jack_fresh)$C_area,perC_area_jack_val_pred_fresh),
                                     perRMSE=percentRMSD(meta(val_jack_fresh)$C_area,perC_area_jack_val_pred_fresh,0.025,0.975),
                                     bias=mean(perC_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$C_area,na.rm=T))
  
  perN_area_jack_val_pred_fresh<-as.vector(predict(perN_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_perN_area_fresh)[,,1])
  perN_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$N_area~perN_area_jack_val_pred_fresh)
  perN_area_jack_stats_fresh[[i]]<-c(R2=summary(perN_area_jack_val_fit_fresh)$r.squared,
                                     RMSE=RMSD(meta(val_jack_fresh)$N_area,perN_area_jack_val_pred_fresh),
                                     perRMSE=percentRMSD(meta(val_jack_fresh)$N_area,perN_area_jack_val_pred_fresh,0.025,0.975),
                                     bias=mean(perN_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$N_area,na.rm=T))
  
  chlA_area_jack_val_pred_fresh<-as.vector(predict(chlA_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_chlA_area_fresh)[,,1])
  chlA_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$chlA_area~chlA_area_jack_val_pred_fresh)
  chlA_area_jack_stats_fresh[[i]]<-c(R2=summary(chlA_area_jack_val_fit_fresh)$r.squared,
                                     RMSE=RMSD(meta(val_jack_fresh)$chlA_area,chlA_area_jack_val_pred_fresh),
                                     perRMSE=percentRMSD(meta(val_jack_fresh)$chlA_area,chlA_area_jack_val_pred_fresh,0.025,0.975),
                                     bias=mean(chlA_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$chlA_area,na.rm=T))
  
  chlB_area_jack_val_pred_fresh<-as.vector(predict(chlB_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_chlB_area_fresh)[,,1])
  chlB_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$chlB_area~chlB_area_jack_val_pred_fresh)
  chlB_area_jack_stats_fresh[[i]]<-c(R2=summary(chlB_area_jack_val_fit_fresh)$r.squared,
                                     RMSE=RMSD(meta(val_jack_fresh)$chlB_area,chlB_area_jack_val_pred_fresh),
                                     perRMSE=percentRMSD(meta(val_jack_fresh)$chlB_area,chlB_area_jack_val_pred_fresh,0.025,0.975),
                                     bias=mean(chlB_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$chlB_area,na.rm=T))
  
  car_area_jack_val_pred_fresh<-as.vector(predict(car_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_car_area_fresh)[,,1])
  car_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$car_area~car_area_jack_val_pred_fresh)
  car_area_jack_stats_fresh[[i]]<-c(R2=summary(car_area_jack_val_fit_fresh)$r.squared,
                                    RMSE=RMSD(meta(val_jack_fresh)$car_area,car_area_jack_val_pred_fresh),
                                    perRMSE=percentRMSD(meta(val_jack_fresh)$car_area,car_area_jack_val_pred_fresh,0.025,0.975),
                                    bias=mean(car_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$car_area,na.rm=T))
  
  Al_area_jack_val_pred_fresh<-as.vector(predict(Al_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Al_area_fresh)[,,1])
  Al_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Al_area~Al_area_jack_val_pred_fresh)
  Al_area_jack_stats_fresh[[i]]<-c(R2=summary(Al_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Al_area,Al_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Al_area,Al_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Al_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Al_area,na.rm=T))
  
  Ca_area_jack_val_pred_fresh<-as.vector(predict(Ca_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Ca_area_fresh)[,,1])
  Ca_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Ca_area~Ca_area_jack_val_pred_fresh)
  Ca_area_jack_stats_fresh[[i]]<-c(R2=summary(Ca_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Ca_area,Ca_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Ca_area,Ca_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Ca_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Ca_area,na.rm=T))
  
  Cu_area_jack_val_pred_fresh<-as.vector(predict(Cu_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Cu_area_fresh)[,,1])
  Cu_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Cu_area~Cu_area_jack_val_pred_fresh)
  Cu_area_jack_stats_fresh[[i]]<-c(R2=summary(Cu_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Cu_area,Cu_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Cu_area,Cu_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Cu_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Cu_area,na.rm=T))
  
  Fe_area_jack_val_pred_fresh<-as.vector(predict(Fe_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Fe_area_fresh)[,,1])
  Fe_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Fe_area~Fe_area_jack_val_pred_fresh)
  Fe_area_jack_stats_fresh[[i]]<-c(R2=summary(Fe_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Fe_area,Fe_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Fe_area,Fe_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Fe_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Fe_area,na.rm=T))
  
  K_area_jack_val_pred_fresh<-as.vector(predict(K_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_K_area_fresh)[,,1])
  K_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$K_area~K_area_jack_val_pred_fresh)
  K_area_jack_stats_fresh[[i]]<-c(R2=summary(K_area_jack_val_fit_fresh)$r.squared,
                                  RMSE=RMSD(meta(val_jack_fresh)$K_area,K_area_jack_val_pred_fresh),
                                  perRMSE=percentRMSD(meta(val_jack_fresh)$K_area,K_area_jack_val_pred_fresh,0.025,0.975),
                                  bias=mean(K_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$K_area,na.rm=T))
  
  Mg_area_jack_val_pred_fresh<-as.vector(predict(Mg_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Mg_area_fresh)[,,1])
  Mg_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Mg_area~Mg_area_jack_val_pred_fresh)
  Mg_area_jack_stats_fresh[[i]]<-c(R2=summary(Mg_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Mg_area,Mg_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Mg_area,Mg_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Mg_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Mg_area,na.rm=T))
  
  Mn_area_jack_val_pred_fresh<-as.vector(predict(Mn_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Mn_area_fresh)[,,1])
  Mn_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Mn_area~Mn_area_jack_val_pred_fresh)
  Mn_area_jack_stats_fresh[[i]]<-c(R2=summary(Mn_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Mn_area,Mn_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Mn_area,Mn_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Mn_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Mn_area,na.rm=T))
  
  Na_area_jack_val_pred_fresh<-as.vector(predict(Na_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Na_area_fresh)[,,1])
  Na_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Na_area~Na_area_jack_val_pred_fresh)
  Na_area_jack_stats_fresh[[i]]<-c(R2=summary(Na_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Na_area,Na_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Na_area,Na_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Na_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Na_area,na.rm=T))
  
  P_area_jack_val_pred_fresh<-as.vector(predict(P_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_P_area_fresh)[,,1])
  P_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$P_area~P_area_jack_val_pred_fresh)
  P_area_jack_stats_fresh[[i]]<-c(R2=summary(P_area_jack_val_fit_fresh)$r.squared,
                                  RMSE=RMSD(meta(val_jack_fresh)$P_area,P_area_jack_val_pred_fresh),
                                  perRMSE=percentRMSD(meta(val_jack_fresh)$P_area,P_area_jack_val_pred_fresh,0.025,0.975),
                                  bias=mean(P_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$P_area,na.rm=T))
  
  Zn_area_jack_val_pred_fresh<-as.vector(predict(Zn_area_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_Zn_area_fresh)[,,1])
  Zn_area_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$Zn_area~Zn_area_jack_val_pred_fresh)
  Zn_area_jack_stats_fresh[[i]]<-c(R2=summary(Zn_area_jack_val_fit_fresh)$r.squared,
                                   RMSE=RMSD(meta(val_jack_fresh)$Zn_area,Zn_area_jack_val_pred_fresh),
                                   perRMSE=percentRMSD(meta(val_jack_fresh)$Zn_area,Zn_area_jack_val_pred_fresh,0.025,0.975),
                                   bias=mean(Zn_area_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$Zn_area,na.rm=T))
  
  solubles_area_jack_coefs_fresh[[i]]<-as.vector(coef(solubles_area_fresh_jack,ncomp=ncomp_solubles_area_fresh,intercept=TRUE))
  hemicellulose_area_jack_coefs_fresh[[i]]<-as.vector(coef(hemicellulose_area_fresh_jack,ncomp=ncomp_hemicellulose_area_fresh,intercept=TRUE))
  cellulose_area_jack_coefs_fresh[[i]]<-as.vector(coef(cellulose_area_fresh_jack,ncomp=ncomp_cellulose_area_fresh,intercept=TRUE))
  lignin_area_jack_coefs_fresh[[i]]<-as.vector(coef(lignin_area_fresh_jack,ncomp=ncomp_lignin_area_fresh,intercept=TRUE))
  perC_area_jack_coefs_fresh[[i]]<-as.vector(coef(perC_area_fresh_jack,ncomp=ncomp_perC_area_fresh,intercept=TRUE))
  perN_area_jack_coefs_fresh[[i]]<-as.vector(coef(perN_area_fresh_jack,ncomp=ncomp_perN_area_fresh,intercept=TRUE))
  chlA_area_jack_coefs_fresh[[i]]<-as.vector(coef(chlA_area_fresh_jack,ncomp=ncomp_chlA_area_fresh,intercept=TRUE))
  chlB_area_jack_coefs_fresh[[i]]<-as.vector(coef(chlB_area_fresh_jack,ncomp=ncomp_chlB_area_fresh,intercept=TRUE))
  car_area_jack_coefs_fresh[[i]]<-as.vector(coef(car_area_fresh_jack,ncomp=ncomp_car_area_fresh,intercept=TRUE))
  Al_area_jack_coefs_fresh[[i]]<-as.vector(coef(Al_area_fresh_jack,ncomp=ncomp_Al_area_fresh,intercept=TRUE))
  Ca_area_jack_coefs_fresh[[i]]<-as.vector(coef(Ca_area_fresh_jack,ncomp=ncomp_Ca_area_fresh,intercept=TRUE))
  Cu_area_jack_coefs_fresh[[i]]<-as.vector(coef(Cu_area_fresh_jack,ncomp=ncomp_Cu_area_fresh,intercept=TRUE))
  Fe_area_jack_coefs_fresh[[i]]<-as.vector(coef(Fe_area_fresh_jack,ncomp=ncomp_Fe_area_fresh,intercept=TRUE))
  K_area_jack_coefs_fresh[[i]]<-as.vector(coef(K_area_fresh_jack,ncomp=ncomp_K_area_fresh,intercept=TRUE))
  Mg_area_jack_coefs_fresh[[i]]<-as.vector(coef(Mg_area_fresh_jack,ncomp=ncomp_Mg_area_fresh,intercept=TRUE))
  Mn_area_jack_coefs_fresh[[i]]<-as.vector(coef(Mn_area_fresh_jack,ncomp=ncomp_Mn_area_fresh,intercept=TRUE))
  Na_area_jack_coefs_fresh[[i]]<-as.vector(coef(Na_area_fresh_jack,ncomp=ncomp_Na_area_fresh,intercept=TRUE))
  P_area_jack_coefs_fresh[[i]]<-as.vector(coef(P_area_fresh_jack,ncomp=ncomp_P_area_fresh,intercept=TRUE))
  Zn_area_jack_coefs_fresh[[i]]<-as.vector(coef(Zn_area_fresh_jack,ncomp=ncomp_Zn_area_fresh,intercept=TRUE))
}

solubles_area_jack_pred_fresh<-apply.coefs(solubles_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
solubles_area_jack_stat_fresh<-t(apply(solubles_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_area_jack_df_fresh<-data.frame(pred_mean=solubles_area_jack_stat_fresh[,1],
                                        pred_low=solubles_area_jack_stat_fresh[,2],
                                        pred_high=solubles_area_jack_stat_fresh[,3],
                                        Measured=meta(fresh_spec_all_test)$solubles_area,
                                        ncomp=ncomp_solubles_area_fresh,
                                        Project=meta(fresh_spec_all_test)$Project,
                                        GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                        Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                        ID=meta(fresh_spec_all_test)$ID)

hemicellulose_area_jack_pred_fresh<-apply.coefs(hemicellulose_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
hemicellulose_area_jack_stat_fresh<-t(apply(hemicellulose_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_area_jack_df_fresh<-data.frame(pred_mean=hemicellulose_area_jack_stat_fresh[,1],
                                             pred_low=hemicellulose_area_jack_stat_fresh[,2],
                                             pred_high=hemicellulose_area_jack_stat_fresh[,3],
                                             Measured=meta(fresh_spec_all_test)$hemicellulose_area,
                                             ncomp=ncomp_hemicellulose_area_fresh,
                                             Project=meta(fresh_spec_all_test)$Project,
                                             GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                             Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                             ID=meta(fresh_spec_all_test)$ID)

cellulose_area_jack_pred_fresh<-apply.coefs(cellulose_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
cellulose_area_jack_stat_fresh<-t(apply(cellulose_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_area_jack_df_fresh<-data.frame(pred_mean=cellulose_area_jack_stat_fresh[,1],
                                         pred_low=cellulose_area_jack_stat_fresh[,2],
                                         pred_high=cellulose_area_jack_stat_fresh[,3],
                                         Measured=meta(fresh_spec_all_test)$cellulose_area,
                                         ncomp=ncomp_cellulose_area_fresh,
                                         Project=meta(fresh_spec_all_test)$Project,
                                         GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                         Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                         ID=meta(fresh_spec_all_test)$ID)

lignin_area_jack_pred_fresh<-apply.coefs(lignin_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
lignin_area_jack_stat_fresh<-t(apply(lignin_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_area_jack_df_fresh<-data.frame(pred_mean=lignin_area_jack_stat_fresh[,1],
                                      pred_low=lignin_area_jack_stat_fresh[,2],
                                      pred_high=lignin_area_jack_stat_fresh[,3],
                                      Measured=meta(fresh_spec_all_test)$lignin_area,
                                      ncomp=ncomp_lignin_area_fresh,
                                      Project=meta(fresh_spec_all_test)$Project,
                                      GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                      Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                      ID=meta(fresh_spec_all_test)$ID)

perC_area_jack_pred_fresh<-apply.coefs(perC_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
perC_area_jack_stat_fresh<-t(apply(perC_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_area_jack_df_fresh<-data.frame(pred_mean=perC_area_jack_stat_fresh[,1],
                                    pred_low=perC_area_jack_stat_fresh[,2],
                                    pred_high=perC_area_jack_stat_fresh[,3],
                                    Measured=meta(fresh_spec_all_test)$C_area,
                                    ncomp=ncomp_perC_area_fresh,
                                    Project=meta(fresh_spec_all_test)$Project,
                                    GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                    Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                    ID=meta(fresh_spec_all_test)$ID)

perN_area_jack_pred_fresh<-apply.coefs(perN_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
perN_area_jack_stat_fresh<-t(apply(perN_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_area_jack_df_fresh<-data.frame(pred_mean=perN_area_jack_stat_fresh[,1],
                                    pred_low=perN_area_jack_stat_fresh[,2],
                                    pred_high=perN_area_jack_stat_fresh[,3],
                                    Measured=meta(fresh_spec_all_test)$N_area,
                                    ncomp=ncomp_perN_area_fresh,
                                    Project=meta(fresh_spec_all_test)$Project,
                                    GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                    Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                    ID=meta(fresh_spec_all_test)$ID)

chlA_area_jack_pred_fresh<-apply.coefs(chlA_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
chlA_area_jack_stat_fresh<-t(apply(chlA_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_area_jack_df_fresh<-data.frame(pred_mean=chlA_area_jack_stat_fresh[,1],
                                    pred_low=chlA_area_jack_stat_fresh[,2],
                                    pred_high=chlA_area_jack_stat_fresh[,3],
                                    Measured=meta(fresh_spec_all_test)$chlA_area,
                                    ncomp=ncomp_chlA_area_fresh,
                                    Project=meta(fresh_spec_all_test)$Project,
                                    GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                    Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                    ID=meta(fresh_spec_all_test)$ID)

chlB_area_jack_pred_fresh<-apply.coefs(chlB_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
chlB_area_jack_stat_fresh<-t(apply(chlB_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_area_jack_df_fresh<-data.frame(pred_mean=chlB_area_jack_stat_fresh[,1],
                                    pred_low=chlB_area_jack_stat_fresh[,2],
                                    pred_high=chlB_area_jack_stat_fresh[,3],
                                    Measured=meta(fresh_spec_all_test)$chlB_area,
                                    ncomp=ncomp_chlB_area_fresh,
                                    Project=meta(fresh_spec_all_test)$Project,
                                    GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                    Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                    ID=meta(fresh_spec_all_test)$ID)

car_area_jack_pred_fresh<-apply.coefs(car_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
car_area_jack_stat_fresh<-t(apply(car_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_area_jack_df_fresh<-data.frame(pred_mean=car_area_jack_stat_fresh[,1],
                                   pred_low=car_area_jack_stat_fresh[,2],
                                   pred_high=car_area_jack_stat_fresh[,3],
                                   Measured=meta(fresh_spec_all_test)$car_area,
                                   ncomp=ncomp_car_area_fresh,
                                   Project=meta(fresh_spec_all_test)$Project,
                                   GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                   Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                   ID=meta(fresh_spec_all_test)$ID)

Al_area_jack_pred_fresh<-apply.coefs(Al_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Al_area_jack_stat_fresh<-t(apply(Al_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_area_jack_df_fresh<-data.frame(pred_mean=Al_area_jack_stat_fresh[,1],
                                  pred_low=Al_area_jack_stat_fresh[,2],
                                  pred_high=Al_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Al_area,
                                  ncomp=ncomp_Al_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

Ca_area_jack_pred_fresh<-apply.coefs(Ca_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Ca_area_jack_stat_fresh<-t(apply(Ca_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_area_jack_df_fresh<-data.frame(pred_mean=Ca_area_jack_stat_fresh[,1],
                                  pred_low=Ca_area_jack_stat_fresh[,2],
                                  pred_high=Ca_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Ca_area,
                                  ncomp=ncomp_Ca_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

Cu_area_jack_pred_fresh<-apply.coefs(Cu_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Cu_area_jack_stat_fresh<-t(apply(Cu_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_area_jack_df_fresh<-data.frame(pred_mean=Cu_area_jack_stat_fresh[,1],
                                  pred_low=Cu_area_jack_stat_fresh[,2],
                                  pred_high=Cu_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Cu_area,
                                  ncomp=ncomp_Cu_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

Fe_area_jack_pred_fresh<-apply.coefs(Fe_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Fe_area_jack_stat_fresh<-t(apply(Fe_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_area_jack_df_fresh<-data.frame(pred_mean=Fe_area_jack_stat_fresh[,1],
                                  pred_low=Fe_area_jack_stat_fresh[,2],
                                  pred_high=Fe_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Fe_area,
                                  ncomp=ncomp_Fe_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

K_area_jack_pred_fresh<-apply.coefs(K_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
K_area_jack_stat_fresh<-t(apply(K_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_area_jack_df_fresh<-data.frame(pred_mean=K_area_jack_stat_fresh[,1],
                                 pred_low=K_area_jack_stat_fresh[,2],
                                 pred_high=K_area_jack_stat_fresh[,3],
                                 Measured=meta(fresh_spec_all_test)$K_area,
                                 ncomp=ncomp_K_area_fresh,
                                 Project=meta(fresh_spec_all_test)$Project,
                                 GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                 Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                 ID=meta(fresh_spec_all_test)$ID)

Mg_area_jack_pred_fresh<-apply.coefs(Mg_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Mg_area_jack_stat_fresh<-t(apply(Mg_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_area_jack_df_fresh<-data.frame(pred_mean=Mg_area_jack_stat_fresh[,1],
                                  pred_low=Mg_area_jack_stat_fresh[,2],
                                  pred_high=Mg_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Mg_area,
                                  ncomp=ncomp_Mg_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

Mn_area_jack_pred_fresh<-apply.coefs(Mn_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Mn_area_jack_stat_fresh<-t(apply(Mn_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_area_jack_df_fresh<-data.frame(pred_mean=Mn_area_jack_stat_fresh[,1],
                                  pred_low=Mn_area_jack_stat_fresh[,2],
                                  pred_high=Mn_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Mn_area,
                                  ncomp=ncomp_Mn_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

Na_area_jack_pred_fresh<-apply.coefs(Na_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Na_area_jack_stat_fresh<-t(apply(Na_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_area_jack_df_fresh<-data.frame(pred_mean=Na_area_jack_stat_fresh[,1],
                                  pred_low=Na_area_jack_stat_fresh[,2],
                                  pred_high=Na_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Na_area,
                                  ncomp=ncomp_Na_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

P_area_jack_pred_fresh<-apply.coefs(P_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
P_area_jack_stat_fresh<-t(apply(P_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_area_jack_df_fresh<-data.frame(pred_mean=P_area_jack_stat_fresh[,1],
                                 pred_low=P_area_jack_stat_fresh[,2],
                                 pred_high=P_area_jack_stat_fresh[,3],
                                 Measured=meta(fresh_spec_all_test)$P_area,
                                 ncomp=ncomp_P_area_fresh,
                                 Project=meta(fresh_spec_all_test)$Project,
                                 GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                 Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                 ID=meta(fresh_spec_all_test)$ID)

Zn_area_jack_pred_fresh<-apply.coefs(Zn_area_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
Zn_area_jack_stat_fresh<-t(apply(Zn_area_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_area_jack_df_fresh<-data.frame(pred_mean=Zn_area_jack_stat_fresh[,1],
                                  pred_low=Zn_area_jack_stat_fresh[,2],
                                  pred_high=Zn_area_jack_stat_fresh[,3],
                                  Measured=meta(fresh_spec_all_test)$Zn_area,
                                  ncomp=ncomp_Zn_area_fresh,
                                  Project=meta(fresh_spec_all_test)$Project,
                                  GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                                  Discoloration=meta(fresh_spec_all_test)$Discoloration,
                                  ID=meta(fresh_spec_all_test)$ID)

##########################################################
## save output

fresh_area_jack_coef_list<-list(C=perC_area_jack_coefs_fresh,
                                N=perN_area_jack_coefs_fresh,
                                sol=solubles_area_jack_coefs_fresh,
                                hemi=hemicellulose_area_jack_coefs_fresh,
                                cell=cellulose_area_jack_coefs_fresh,
                                lign=lignin_area_jack_coefs_fresh,
                                chlA=chlA_area_jack_coefs_fresh,
                                chlB=chlB_area_jack_coefs_fresh,
                                car=car_area_jack_coefs_fresh,
                                Al=Al_area_jack_coefs_fresh,
                                Ca=Ca_area_jack_coefs_fresh,
                                Cu=Cu_area_jack_coefs_fresh,
                                Fe=Fe_area_jack_coefs_fresh,
                                K=K_area_jack_coefs_fresh,
                                Mg=Mg_area_jack_coefs_fresh,
                                Mn=Mn_area_jack_coefs_fresh,
                                Na=Na_area_jack_coefs_fresh,
                                P=P_area_jack_coefs_fresh,
                                Zn=Zn_area_jack_coefs_fresh)
saveRDS(fresh_area_jack_coef_list,"SavedResults/fresh_area_jack_coefs_list.rds")

fresh_area_jack_df_list<-list(C=perC_area_jack_df_fresh,
                              N=perN_area_jack_df_fresh,
                              sol=solubles_area_jack_df_fresh,
                              hemi=hemicellulose_area_jack_df_fresh,
                              cell=cellulose_area_jack_df_fresh,
                              lign=lignin_area_jack_df_fresh,
                              chlA=chlA_area_jack_df_fresh,
                              chlB=chlB_area_jack_df_fresh,
                              car=car_area_jack_df_fresh,
                              Al=Al_area_jack_df_fresh,
                              Ca=Ca_area_jack_df_fresh,
                              Cu=Cu_area_jack_df_fresh,
                              Fe=Fe_area_jack_df_fresh,
                              K=K_area_jack_df_fresh,
                              Mg=Mg_area_jack_df_fresh,
                              Mn=Mn_area_jack_df_fresh,
                              Na=Na_area_jack_df_fresh,
                              P=P_area_jack_df_fresh,
                              Zn=Zn_area_jack_df_fresh)
saveRDS(fresh_area_jack_df_list,"SavedResults/fresh_area_jack_df_list.rds")
