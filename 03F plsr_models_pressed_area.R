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

pressed_spec_all_train<-readRDS("ProcessedSpectralData/pressed_spec_all_train.rds")
pressed_spec_all_test<-readRDS("ProcessedSpectralData/pressed_spec_all_test.rds")

################################################
## building calibration models

perC_area_pressed<-plsr(meta(pressed_spec_all_train)$C_area~as.matrix(pressed_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perC_area_pressed <- selectNcomp(perC_area_pressed, method = "onesigma", plot = FALSE)
perC_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$C_area))
perC_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[perC_area_pressed_valid],
                                 Species=meta(pressed_spec_all_train)$Species[perC_area_pressed_valid],
                                 Project=meta(pressed_spec_all_train)$Project[perC_area_pressed_valid],
                                 Stage=meta(pressed_spec_all_train)$Stage[perC_area_pressed_valid],
                                 GrowthForm=meta(pressed_spec_all_train)$GrowthForm[perC_area_pressed_valid],
                                 measured=meta(pressed_spec_all_train)$C_area[perC_area_pressed_valid],
                                 val_pred=perC_area_pressed$validation$pred[,,ncomp_perC_area_pressed])
ggplot(perC_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %C from pressed-leaf spectra")

perN_area_pressed<-plsr(meta(pressed_spec_all_train)$N_area~as.matrix(pressed_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_perN_area_pressed <- selectNcomp(perN_area_pressed, method = "onesigma", plot = FALSE)
perN_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$N_area))
perN_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[perN_area_pressed_valid],
                                 Species=meta(pressed_spec_all_train)$Species[perN_area_pressed_valid],
                                 Project=meta(pressed_spec_all_train)$Project[perN_area_pressed_valid],
                                 Stage=meta(pressed_spec_all_train)$Stage[perN_area_pressed_valid],
                                 GrowthForm=meta(pressed_spec_all_train)$GrowthForm[perN_area_pressed_valid],
                                 measured=meta(pressed_spec_all_train)$N_area[perN_area_pressed_valid],
                                 val_pred=perN_area_pressed$validation$pred[,,ncomp_perN_area_pressed])
ggplot(perN_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting %N from pressed-leaf spectra")+guides(color=F)

chlA_area_pressed<-plsr(meta(pressed_spec_all_train)$chlA_area~as.matrix(pressed_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_area_pressed <- selectNcomp(chlA_area_pressed, method = "onesigma", plot = FALSE)
chlA_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$chlA_area))
chlA_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[chlA_area_pressed_valid],
                                 Species=meta(pressed_spec_all_train)$Species[chlA_area_pressed_valid],
                                 Project=meta(pressed_spec_all_train)$Project[chlA_area_pressed_valid],
                                 Stage=meta(pressed_spec_all_train)$Stage[chlA_area_pressed_valid],
                                 GrowthForm=meta(pressed_spec_all_train)$GrowthForm[chlA_area_pressed_valid],
                                 measured=meta(pressed_spec_all_train)$chlA_area[chlA_area_pressed_valid],
                                 val_pred=chlA_area_pressed$validation$pred[,,ncomp_chlA_area_pressed])
ggplot(chlA_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl a from pressed-leaf spectra")

chlB_area_pressed<-plsr(meta(pressed_spec_all_train)$chlB_area~as.matrix(pressed_spec_all_train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_area_pressed <- selectNcomp(chlB_area_pressed, method = "onesigma", plot = FALSE)
chlB_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$chlB_area))
chlB_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[chlB_area_pressed_valid],
                                 Species=meta(pressed_spec_all_train)$Species[chlB_area_pressed_valid],
                                 Project=meta(pressed_spec_all_train)$Project[chlB_area_pressed_valid],
                                 Stage=meta(pressed_spec_all_train)$Stage[chlB_area_pressed_valid],
                                 GrowthForm=meta(pressed_spec_all_train)$GrowthForm[chlB_area_pressed_valid],
                                 measured=meta(pressed_spec_all_train)$chlB_area[chlB_area_pressed_valid],
                                 val_pred=chlB_area_pressed$validation$pred[,,ncomp_chlB_area_pressed])
ggplot(chlB_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Chl b from pressed-leaf spectra")

car_area_pressed<-plsr(meta(pressed_spec_all_train)$car_area~as.matrix(pressed_spec_all_train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_area_pressed <- selectNcomp(car_area_pressed, method = "onesigma", plot = FALSE)
car_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$car_area))
car_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[car_area_pressed_valid],
                                Species=meta(pressed_spec_all_train)$Species[car_area_pressed_valid],
                                Project=meta(pressed_spec_all_train)$Project[car_area_pressed_valid],
                                Stage=meta(pressed_spec_all_train)$Stage[car_area_pressed_valid],
                                GrowthForm=meta(pressed_spec_all_train)$GrowthForm[car_area_pressed_valid],
                                measured=meta(pressed_spec_all_train)$car_area[car_area_pressed_valid],
                                val_pred=car_area_pressed$validation$pred[,,ncomp_car_area_pressed])
ggplot(car_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting total car from pressed-leaf spectra")

solubles_area_pressed<-plsr(meta(pressed_spec_all_train)$solubles_area~as.matrix(pressed_spec_all_train),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_area_pressed <- selectNcomp(solubles_area_pressed, method = "onesigma", plot = FALSE)
solubles_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$solubles_area))
solubles_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[solubles_area_pressed_valid],
                                     Species=meta(pressed_spec_all_train)$Species[solubles_area_pressed_valid],
                                     Project=meta(pressed_spec_all_train)$Project[solubles_area_pressed_valid],
                                     Stage=meta(pressed_spec_all_train)$Stage[solubles_area_pressed_valid],
                                     GrowthForm=meta(pressed_spec_all_train)$GrowthForm[solubles_area_pressed_valid],
                                     measured=meta(pressed_spec_all_train)$solubles_area[solubles_area_pressed_valid],
                                     val_pred=solubles_area_pressed$validation$pred[,,ncomp_solubles_area_pressed])
ggplot(solubles_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting solubles from pressed-leaf spectra")

hemicellulose_area_pressed<-plsr(meta(pressed_spec_all_train)$hemicellulose_area~as.matrix(pressed_spec_all_train),
                               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_area_pressed <- selectNcomp(hemicellulose_area_pressed, method = "onesigma", plot = FALSE)
hemicellulose_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$hemicellulose_area))
hemicellulose_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[hemicellulose_area_pressed_valid],
                                          Species=meta(pressed_spec_all_train)$Species[hemicellulose_area_pressed_valid],
                                          Project=meta(pressed_spec_all_train)$Project[hemicellulose_area_pressed_valid],
                                          Stage=meta(pressed_spec_all_train)$Stage[hemicellulose_area_pressed_valid],
                                          GrowthForm=meta(pressed_spec_all_train)$GrowthForm[hemicellulose_area_pressed_valid],
                                          measured=meta(pressed_spec_all_train)$hemicellulose_area[hemicellulose_area_pressed_valid],
                                          val_pred=hemicellulose_area_pressed$validation$pred[,,ncomp_hemicellulose_area_pressed])
ggplot(hemicellulose_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting hemicellulose from pressed-leaf spectra")

cellulose_area_pressed<-plsr(meta(pressed_spec_all_train)$cellulose_area~as.matrix(pressed_spec_all_train),
                           ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_area_pressed <- selectNcomp(cellulose_area_pressed, method = "onesigma", plot = FALSE)
cellulose_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$cellulose_area))
cellulose_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[cellulose_area_pressed_valid],
                                      Species=meta(pressed_spec_all_train)$Species[cellulose_area_pressed_valid],
                                      Project=meta(pressed_spec_all_train)$Project[cellulose_area_pressed_valid],
                                      Stage=meta(pressed_spec_all_train)$Stage[cellulose_area_pressed_valid],
                                      GrowthForm=meta(pressed_spec_all_train)$GrowthForm[cellulose_area_pressed_valid],
                                      measured=meta(pressed_spec_all_train)$cellulose_area[cellulose_area_pressed_valid],
                                      val_pred=cellulose_area_pressed$validation$pred[,,ncomp_cellulose_area_pressed])
ggplot(cellulose_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting cellulose from pressed-leaf spectra")

lignin_area_pressed<-plsr(meta(pressed_spec_all_train)$lignin_area~as.matrix(pressed_spec_all_train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_area_pressed <- selectNcomp(lignin_area_pressed, method = "onesigma", plot = FALSE)
lignin_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$lignin_area))
lignin_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[lignin_area_pressed_valid],
                                   Species=meta(pressed_spec_all_train)$Species[lignin_area_pressed_valid],
                                   Project=meta(pressed_spec_all_train)$Project[lignin_area_pressed_valid],
                                   Stage=meta(pressed_spec_all_train)$Stage[lignin_area_pressed_valid],
                                   GrowthForm=meta(pressed_spec_all_train)$GrowthForm[lignin_area_pressed_valid],
                                   measured=meta(pressed_spec_all_train)$lignin_area[lignin_area_pressed_valid],
                                   val_pred=lignin_area_pressed$validation$pred[,,ncomp_lignin_area_pressed])
ggplot(lignin_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting lignin from pressed-leaf spectra")

Al_area_pressed<-plsr(meta(pressed_spec_all_train)$Al_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_area_pressed <- selectNcomp(Al_area_pressed, method = "onesigma", plot = FALSE)
Al_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Al_area))
Al_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Al_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Al_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Al_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Al_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Al_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Al_area[Al_area_pressed_valid],
                               val_pred=Al_area_pressed$validation$pred[,,ncomp_Al_area_pressed])
ggplot(Al_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Al from pressed-leaf spectra")

Ca_area_pressed<-plsr(meta(pressed_spec_all_train)$Ca_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_area_pressed <- selectNcomp(Ca_area_pressed, method = "onesigma", plot = FALSE)
Ca_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Ca_area))
Ca_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Ca_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Ca_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Ca_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Ca_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Ca_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Ca_area[Ca_area_pressed_valid],
                               val_pred=Ca_area_pressed$validation$pred[,,ncomp_Ca_area_pressed])
ggplot(Ca_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Ca from pressed-leaf spectra")

Cu_area_pressed<-plsr(meta(pressed_spec_all_train)$Cu_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_area_pressed <- selectNcomp(Cu_area_pressed, method = "onesigma", plot = FALSE)
Cu_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Cu_area))
Cu_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Cu_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Cu_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Cu_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Cu_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Cu_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Cu_area[Cu_area_pressed_valid],
                               val_pred=Cu_area_pressed$validation$pred[,,ncomp_Cu_area_pressed])
ggplot(Cu_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Cu from pressed-leaf spectra")

Fe_area_pressed<-plsr(meta(pressed_spec_all_train)$Fe_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_area_pressed <- selectNcomp(Fe_area_pressed, method = "onesigma", plot = FALSE)
Fe_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Fe_area))
Fe_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Fe_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Fe_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Fe_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Fe_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Fe_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Fe_area[Fe_area_pressed_valid],
                               val_pred=Fe_area_pressed$validation$pred[,,ncomp_Fe_area_pressed])
ggplot(Fe_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Fe from pressed-leaf spectra")

K_area_pressed<-plsr(meta(pressed_spec_all_train)$K_area~as.matrix(pressed_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_area_pressed <- selectNcomp(K_area_pressed, method = "onesigma", plot = FALSE)
K_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$K_area))
K_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[K_area_pressed_valid],
                              Species=meta(pressed_spec_all_train)$Species[K_area_pressed_valid],
                              Project=meta(pressed_spec_all_train)$Project[K_area_pressed_valid],
                              Stage=meta(pressed_spec_all_train)$Stage[K_area_pressed_valid],
                              GrowthForm=meta(pressed_spec_all_train)$GrowthForm[K_area_pressed_valid],
                              measured=meta(pressed_spec_all_train)$K_area[K_area_pressed_valid],
                              val_pred=K_area_pressed$validation$pred[,,ncomp_K_area_pressed])
ggplot(K_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting K from pressed-leaf spectra")

Mg_area_pressed<-plsr(meta(pressed_spec_all_train)$Mg_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_area_pressed <- selectNcomp(Mg_area_pressed, method = "onesigma", plot = FALSE)
Mg_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Mg_area))
Mg_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Mg_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Mg_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Mg_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Mg_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Mg_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Mg_area[Mg_area_pressed_valid],
                               val_pred=Mg_area_pressed$validation$pred[,,ncomp_Mg_area_pressed])
ggplot(Mg_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mg from pressed-leaf spectra")

Mn_area_pressed<-plsr(meta(pressed_spec_all_train)$Mn_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_area_pressed <- selectNcomp(Mn_area_pressed, method = "onesigma", plot = FALSE)
Mn_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Mn_area))
Mn_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Mn_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Mn_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Mn_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Mn_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Mn_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Mn_area[Mn_area_pressed_valid],
                               val_pred=Mn_area_pressed$validation$pred[,,ncomp_Mn_area_pressed])
ggplot(Mn_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Mn from pressed-leaf spectra")

Na_area_pressed<-plsr(meta(pressed_spec_all_train)$Na_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_area_pressed <- selectNcomp(Na_area_pressed, method = "onesigma", plot = FALSE)
Na_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Na_area))
Na_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Na_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Na_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Na_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Na_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Na_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Na_area[Na_area_pressed_valid],
                               val_pred=Na_area_pressed$validation$pred[,,ncomp_Na_area_pressed])
ggplot(Na_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Na from pressed-leaf spectra")

P_area_pressed<-plsr(meta(pressed_spec_all_train)$P_area~as.matrix(pressed_spec_all_train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_area_pressed <- selectNcomp(P_area_pressed, method = "onesigma", plot = FALSE)
P_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$P_area))
P_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[P_area_pressed_valid],
                              Species=meta(pressed_spec_all_train)$Species[P_area_pressed_valid],
                              Project=meta(pressed_spec_all_train)$Project[P_area_pressed_valid],
                              Stage=meta(pressed_spec_all_train)$Stage[P_area_pressed_valid],
                              GrowthForm=meta(pressed_spec_all_train)$GrowthForm[P_area_pressed_valid],
                              measured=meta(pressed_spec_all_train)$P_area[P_area_pressed_valid],
                              val_pred=P_area_pressed$validation$pred[,,ncomp_P_area_pressed])
ggplot(P_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting P from pressed-leaf spectra")

Zn_area_pressed<-plsr(meta(pressed_spec_all_train)$Zn_area~as.matrix(pressed_spec_all_train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_area_pressed <- selectNcomp(Zn_area_pressed, method = "onesigma", plot = FALSE)
Zn_area_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$Zn_area))
Zn_area_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[Zn_area_pressed_valid],
                               Species=meta(pressed_spec_all_train)$Species[Zn_area_pressed_valid],
                               Project=meta(pressed_spec_all_train)$Project[Zn_area_pressed_valid],
                               Stage=meta(pressed_spec_all_train)$Stage[Zn_area_pressed_valid],
                               GrowthForm=meta(pressed_spec_all_train)$GrowthForm[Zn_area_pressed_valid],
                               measured=meta(pressed_spec_all_train)$Zn_area[Zn_area_pressed_valid],
                               val_pred=Zn_area_pressed$validation$pred[,,ncomp_Zn_area_pressed])
ggplot(Zn_area_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting Zn from pressed-leaf spectra")

###############################################
## jackknife tests + prediction of validation data
## pressed leaves

solubles_area_jack_coefs_pressed<-list()
hemicellulose_area_jack_coefs_pressed<-list()
cellulose_area_jack_coefs_pressed<-list()
lignin_area_jack_coefs_pressed<-list()
perC_area_jack_coefs_pressed<-list()
perN_area_jack_coefs_pressed<-list()
chlA_area_jack_coefs_pressed<-list()
chlB_area_jack_coefs_pressed<-list()
car_area_jack_coefs_pressed<-list()
Al_area_jack_coefs_pressed<-list()
Ca_area_jack_coefs_pressed<-list()
Cu_area_jack_coefs_pressed<-list()
Fe_area_jack_coefs_pressed<-list()
K_area_jack_coefs_pressed<-list()
Mg_area_jack_coefs_pressed<-list()
Mn_area_jack_coefs_pressed<-list()
Na_area_jack_coefs_pressed<-list()
P_area_jack_coefs_pressed<-list()
Zn_area_jack_coefs_pressed<-list()

solubles_area_jack_stats_pressed<-list()
hemicellulose_area_jack_stats_pressed<-list()
cellulose_area_jack_stats_pressed<-list()
lignin_area_jack_stats_pressed<-list()
perC_area_jack_stats_pressed<-list()
perN_area_jack_stats_pressed<-list()
chlA_area_jack_stats_pressed<-list()
chlB_area_jack_stats_pressed<-list()
car_area_jack_stats_pressed<-list()
Al_area_jack_stats_pressed<-list()
Ca_area_jack_stats_pressed<-list()
Cu_area_jack_stats_pressed<-list()
Fe_area_jack_stats_pressed<-list()
K_area_jack_stats_pressed<-list()
Mg_area_jack_stats_pressed<-list()
Mn_area_jack_stats_pressed<-list()
Na_area_jack_stats_pressed<-list()
P_area_jack_stats_pressed<-list()
Zn_area_jack_stats_pressed<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_pressed<-nrow(pressed_spec_all_train)
  train_jack_pressed<-sample(1:n_cal_spec_pressed,floor(0.7*n_cal_spec_pressed))
  test_jack_pressed<-setdiff(1:n_cal_spec_pressed,train_jack_pressed)
  
  calib_jack_pressed<-pressed_spec_all_train[train_jack_pressed]
  val_jack_pressed<-pressed_spec_all_train[test_jack_pressed]
  
  solubles_area_pressed_jack<-plsr(meta(calib_jack_pressed)$solubles_area~as.matrix(calib_jack_pressed),
                                 ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_area_pressed_jack<-plsr(meta(calib_jack_pressed)$hemicellulose_area~as.matrix(calib_jack_pressed),
                                      ncomp=30,method = "oscorespls",validation="none")
  cellulose_area_pressed_jack<-plsr(meta(calib_jack_pressed)$cellulose_area~as.matrix(calib_jack_pressed),
                                  ncomp=30,method = "oscorespls",validation="none")
  lignin_area_pressed_jack<-plsr(meta(calib_jack_pressed)$lignin_area~as.matrix(calib_jack_pressed),
                               ncomp=30,method = "oscorespls",validation="none")
  perC_area_pressed_jack<-plsr(meta(calib_jack_pressed)$C_area~as.matrix(calib_jack_pressed),
                             ncomp=30,method = "oscorespls",validation="none")
  perN_area_pressed_jack<-plsr(meta(calib_jack_pressed)$N_area~as.matrix(calib_jack_pressed),
                             ncomp=30,method = "oscorespls",validation="none")
  chlA_area_pressed_jack<-plsr(meta(calib_jack_pressed)$chlA_area~as.matrix(calib_jack_pressed),
                             ncomp=30,method = "oscorespls",validation="none")
  chlB_area_pressed_jack<-plsr(meta(calib_jack_pressed)$chlB_area~as.matrix(calib_jack_pressed),
                             ncomp=30,method = "oscorespls",validation="none")
  car_area_pressed_jack<-plsr(meta(calib_jack_pressed)$car_area~as.matrix(calib_jack_pressed),
                            ncomp=30,method = "oscorespls",validation="none")
  Al_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Al_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  Ca_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Ca_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  Cu_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Cu_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  Fe_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Fe_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  K_area_pressed_jack<-plsr(meta(calib_jack_pressed)$K_area~as.matrix(calib_jack_pressed),
                          ncomp=30,method = "oscorespls",validation="none")
  Mg_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Mg_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  Mn_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Mn_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  Na_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Na_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  P_area_pressed_jack<-plsr(meta(calib_jack_pressed)$P_area~as.matrix(calib_jack_pressed),
                          ncomp=30,method = "oscorespls",validation="none")
  Zn_area_pressed_jack<-plsr(meta(calib_jack_pressed)$Zn_area~as.matrix(calib_jack_pressed),
                           ncomp=30,method = "oscorespls",validation="none")
  
  solubles_area_jack_val_pred_pressed<-as.vector(predict(solubles_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_solubles_area_pressed)[,,1])
  solubles_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$solubles_area~solubles_area_jack_val_pred_pressed)
  solubles_area_jack_stats_pressed[[i]]<-c(R2=summary(solubles_area_jack_val_fit_pressed)$r.squared,
                                         RMSE=RMSD(meta(val_jack_pressed)$solubles_area,solubles_area_jack_val_pred_pressed),
                                         perRMSE=percentRMSD(meta(val_jack_pressed)$solubles_area,solubles_area_jack_val_pred_pressed,0.025,0.975),
                                         bias=mean(solubles_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$solubles_area,na.rm=T))
  
  hemicellulose_area_jack_val_pred_pressed<-as.vector(predict(hemicellulose_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_hemicellulose_area_pressed)[,,1])
  hemicellulose_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$hemicellulose_area~hemicellulose_area_jack_val_pred_pressed)
  hemicellulose_area_jack_stats_pressed[[i]]<-c(R2=summary(hemicellulose_area_jack_val_fit_pressed)$r.squared,
                                              RMSE=RMSD(meta(val_jack_pressed)$hemicellulose_area,hemicellulose_area_jack_val_pred_pressed),
                                              perRMSE=percentRMSD(meta(val_jack_pressed)$hemicellulose_area,hemicellulose_area_jack_val_pred_pressed,0.025,0.975),
                                              bias=mean(hemicellulose_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$hemicellulose_area,na.rm=T))
  
  cellulose_area_jack_val_pred_pressed<-as.vector(predict(cellulose_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_cellulose_area_pressed)[,,1])
  cellulose_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$cellulose_area~cellulose_area_jack_val_pred_pressed)
  cellulose_area_jack_stats_pressed[[i]]<-c(R2=summary(cellulose_area_jack_val_fit_pressed)$r.squared,
                                          RMSE=RMSD(meta(val_jack_pressed)$cellulose_area,cellulose_area_jack_val_pred_pressed),
                                          perRMSE=percentRMSD(meta(val_jack_pressed)$cellulose_area,cellulose_area_jack_val_pred_pressed,0.025,0.975),
                                          bias=mean(cellulose_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$cellulose_area,na.rm=T))
  
  lignin_area_jack_val_pred_pressed<-as.vector(predict(lignin_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_lignin_area_pressed)[,,1])
  lignin_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$lignin_area~lignin_area_jack_val_pred_pressed)
  lignin_area_jack_stats_pressed[[i]]<-c(R2=summary(lignin_area_jack_val_fit_pressed)$r.squared,
                                       RMSE=RMSD(meta(val_jack_pressed)$lignin_area,lignin_area_jack_val_pred_pressed),
                                       perRMSE=percentRMSD(meta(val_jack_pressed)$lignin_area,lignin_area_jack_val_pred_pressed,0.025,0.975),
                                       bias=mean(lignin_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$lignin_area,na.rm=T))
  
  perC_area_jack_val_pred_pressed<-as.vector(predict(perC_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_perC_area_pressed)[,,1])
  perC_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$C_area~perC_area_jack_val_pred_pressed)
  perC_area_jack_stats_pressed[[i]]<-c(R2=summary(perC_area_jack_val_fit_pressed)$r.squared,
                                     RMSE=RMSD(meta(val_jack_pressed)$C_area,perC_area_jack_val_pred_pressed),
                                     perRMSE=percentRMSD(meta(val_jack_pressed)$C_area,perC_area_jack_val_pred_pressed,0.025,0.975),
                                     bias=mean(perC_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$C_area,na.rm=T))
  
  perN_area_jack_val_pred_pressed<-as.vector(predict(perN_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_perN_area_pressed)[,,1])
  perN_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$N_area~perN_area_jack_val_pred_pressed)
  perN_area_jack_stats_pressed[[i]]<-c(R2=summary(perN_area_jack_val_fit_pressed)$r.squared,
                                     RMSE=RMSD(meta(val_jack_pressed)$N_area,perN_area_jack_val_pred_pressed),
                                     perRMSE=percentRMSD(meta(val_jack_pressed)$N_area,perN_area_jack_val_pred_pressed,0.025,0.975),
                                     bias=mean(perN_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$N_area,na.rm=T))
  
  chlA_area_jack_val_pred_pressed<-as.vector(predict(chlA_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_chlA_area_pressed)[,,1])
  chlA_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$chlA_area~chlA_area_jack_val_pred_pressed)
  chlA_area_jack_stats_pressed[[i]]<-c(R2=summary(chlA_area_jack_val_fit_pressed)$r.squared,
                                     RMSE=RMSD(meta(val_jack_pressed)$chlA_area,chlA_area_jack_val_pred_pressed),
                                     perRMSE=percentRMSD(meta(val_jack_pressed)$chlA_area,chlA_area_jack_val_pred_pressed,0.025,0.975),
                                     bias=mean(chlA_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$chlA_area,na.rm=T))
  
  chlB_area_jack_val_pred_pressed<-as.vector(predict(chlB_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_chlB_area_pressed)[,,1])
  chlB_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$chlB_area~chlB_area_jack_val_pred_pressed)
  chlB_area_jack_stats_pressed[[i]]<-c(R2=summary(chlB_area_jack_val_fit_pressed)$r.squared,
                                     RMSE=RMSD(meta(val_jack_pressed)$chlB_area,chlB_area_jack_val_pred_pressed),
                                     perRMSE=percentRMSD(meta(val_jack_pressed)$chlB_area,chlB_area_jack_val_pred_pressed,0.025,0.975),
                                     bias=mean(chlB_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$chlB_area,na.rm=T))
  
  car_area_jack_val_pred_pressed<-as.vector(predict(car_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_car_area_pressed)[,,1])
  car_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$car_area~car_area_jack_val_pred_pressed)
  car_area_jack_stats_pressed[[i]]<-c(R2=summary(car_area_jack_val_fit_pressed)$r.squared,
                                    RMSE=RMSD(meta(val_jack_pressed)$car_area,car_area_jack_val_pred_pressed),
                                    perRMSE=percentRMSD(meta(val_jack_pressed)$car_area,car_area_jack_val_pred_pressed,0.025,0.975),
                                    bias=mean(car_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$car_area,na.rm=T))
  
  Al_area_jack_val_pred_pressed<-as.vector(predict(Al_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Al_area_pressed)[,,1])
  Al_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Al_area~Al_area_jack_val_pred_pressed)
  Al_area_jack_stats_pressed[[i]]<-c(R2=summary(Al_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Al_area,Al_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Al_area,Al_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Al_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Al_area,na.rm=T))
  
  Ca_area_jack_val_pred_pressed<-as.vector(predict(Ca_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Ca_area_pressed)[,,1])
  Ca_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Ca_area~Ca_area_jack_val_pred_pressed)
  Ca_area_jack_stats_pressed[[i]]<-c(R2=summary(Ca_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Ca_area,Ca_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Ca_area,Ca_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Ca_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Ca_area,na.rm=T))
  
  Cu_area_jack_val_pred_pressed<-as.vector(predict(Cu_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Cu_area_pressed)[,,1])
  Cu_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Cu_area~Cu_area_jack_val_pred_pressed)
  Cu_area_jack_stats_pressed[[i]]<-c(R2=summary(Cu_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Cu_area,Cu_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Cu_area,Cu_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Cu_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Cu_area,na.rm=T))
  
  Fe_area_jack_val_pred_pressed<-as.vector(predict(Fe_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Fe_area_pressed)[,,1])
  Fe_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Fe_area~Fe_area_jack_val_pred_pressed)
  Fe_area_jack_stats_pressed[[i]]<-c(R2=summary(Fe_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Fe_area,Fe_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Fe_area,Fe_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Fe_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Fe_area,na.rm=T))
  
  K_area_jack_val_pred_pressed<-as.vector(predict(K_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_K_area_pressed)[,,1])
  K_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$K_area~K_area_jack_val_pred_pressed)
  K_area_jack_stats_pressed[[i]]<-c(R2=summary(K_area_jack_val_fit_pressed)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed)$K_area,K_area_jack_val_pred_pressed),
                                  perRMSE=percentRMSD(meta(val_jack_pressed)$K_area,K_area_jack_val_pred_pressed,0.025,0.975),
                                  bias=mean(K_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$K_area,na.rm=T))
  
  Mg_area_jack_val_pred_pressed<-as.vector(predict(Mg_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Mg_area_pressed)[,,1])
  Mg_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Mg_area~Mg_area_jack_val_pred_pressed)
  Mg_area_jack_stats_pressed[[i]]<-c(R2=summary(Mg_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Mg_area,Mg_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Mg_area,Mg_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Mg_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Mg_area,na.rm=T))
  
  Mn_area_jack_val_pred_pressed<-as.vector(predict(Mn_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Mn_area_pressed)[,,1])
  Mn_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Mn_area~Mn_area_jack_val_pred_pressed)
  Mn_area_jack_stats_pressed[[i]]<-c(R2=summary(Mn_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Mn_area,Mn_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Mn_area,Mn_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Mn_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Mn_area,na.rm=T))
  
  Na_area_jack_val_pred_pressed<-as.vector(predict(Na_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Na_area_pressed)[,,1])
  Na_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Na_area~Na_area_jack_val_pred_pressed)
  Na_area_jack_stats_pressed[[i]]<-c(R2=summary(Na_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Na_area,Na_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Na_area,Na_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Na_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Na_area,na.rm=T))
  
  P_area_jack_val_pred_pressed<-as.vector(predict(P_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_P_area_pressed)[,,1])
  P_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$P_area~P_area_jack_val_pred_pressed)
  P_area_jack_stats_pressed[[i]]<-c(R2=summary(P_area_jack_val_fit_pressed)$r.squared,
                                  RMSE=RMSD(meta(val_jack_pressed)$P_area,P_area_jack_val_pred_pressed),
                                  perRMSE=percentRMSD(meta(val_jack_pressed)$P_area,P_area_jack_val_pred_pressed,0.025,0.975),
                                  bias=mean(P_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$P_area,na.rm=T))
  
  Zn_area_jack_val_pred_pressed<-as.vector(predict(Zn_area_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_Zn_area_pressed)[,,1])
  Zn_area_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$Zn_area~Zn_area_jack_val_pred_pressed)
  Zn_area_jack_stats_pressed[[i]]<-c(R2=summary(Zn_area_jack_val_fit_pressed)$r.squared,
                                   RMSE=RMSD(meta(val_jack_pressed)$Zn_area,Zn_area_jack_val_pred_pressed),
                                   perRMSE=percentRMSD(meta(val_jack_pressed)$Zn_area,Zn_area_jack_val_pred_pressed,0.025,0.975),
                                   bias=mean(Zn_area_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$Zn_area,na.rm=T))
  
  solubles_area_jack_coefs_pressed[[i]]<-as.vector(coef(solubles_area_pressed_jack,ncomp=ncomp_solubles_area_pressed,intercept=TRUE))
  hemicellulose_area_jack_coefs_pressed[[i]]<-as.vector(coef(hemicellulose_area_pressed_jack,ncomp=ncomp_hemicellulose_area_pressed,intercept=TRUE))
  cellulose_area_jack_coefs_pressed[[i]]<-as.vector(coef(cellulose_area_pressed_jack,ncomp=ncomp_cellulose_area_pressed,intercept=TRUE))
  lignin_area_jack_coefs_pressed[[i]]<-as.vector(coef(lignin_area_pressed_jack,ncomp=ncomp_lignin_area_pressed,intercept=TRUE))
  perC_area_jack_coefs_pressed[[i]]<-as.vector(coef(perC_area_pressed_jack,ncomp=ncomp_perC_area_pressed,intercept=TRUE))
  perN_area_jack_coefs_pressed[[i]]<-as.vector(coef(perN_area_pressed_jack,ncomp=ncomp_perN_area_pressed,intercept=TRUE))
  chlA_area_jack_coefs_pressed[[i]]<-as.vector(coef(chlA_area_pressed_jack,ncomp=ncomp_chlA_area_pressed,intercept=TRUE))
  chlB_area_jack_coefs_pressed[[i]]<-as.vector(coef(chlB_area_pressed_jack,ncomp=ncomp_chlB_area_pressed,intercept=TRUE))
  car_area_jack_coefs_pressed[[i]]<-as.vector(coef(car_area_pressed_jack,ncomp=ncomp_car_area_pressed,intercept=TRUE))
  Al_area_jack_coefs_pressed[[i]]<-as.vector(coef(Al_area_pressed_jack,ncomp=ncomp_Al_area_pressed,intercept=TRUE))
  Ca_area_jack_coefs_pressed[[i]]<-as.vector(coef(Ca_area_pressed_jack,ncomp=ncomp_Ca_area_pressed,intercept=TRUE))
  Cu_area_jack_coefs_pressed[[i]]<-as.vector(coef(Cu_area_pressed_jack,ncomp=ncomp_Cu_area_pressed,intercept=TRUE))
  Fe_area_jack_coefs_pressed[[i]]<-as.vector(coef(Fe_area_pressed_jack,ncomp=ncomp_Fe_area_pressed,intercept=TRUE))
  K_area_jack_coefs_pressed[[i]]<-as.vector(coef(K_area_pressed_jack,ncomp=ncomp_K_area_pressed,intercept=TRUE))
  Mg_area_jack_coefs_pressed[[i]]<-as.vector(coef(Mg_area_pressed_jack,ncomp=ncomp_Mg_area_pressed,intercept=TRUE))
  Mn_area_jack_coefs_pressed[[i]]<-as.vector(coef(Mn_area_pressed_jack,ncomp=ncomp_Mn_area_pressed,intercept=TRUE))
  Na_area_jack_coefs_pressed[[i]]<-as.vector(coef(Na_area_pressed_jack,ncomp=ncomp_Na_area_pressed,intercept=TRUE))
  P_area_jack_coefs_pressed[[i]]<-as.vector(coef(P_area_pressed_jack,ncomp=ncomp_P_area_pressed,intercept=TRUE))
  Zn_area_jack_coefs_pressed[[i]]<-as.vector(coef(Zn_area_pressed_jack,ncomp=ncomp_Zn_area_pressed,intercept=TRUE))
}

solubles_area_jack_pred_pressed<-apply.coefs(solubles_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
solubles_area_jack_stat_pressed<-t(apply(solubles_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_area_jack_df_pressed<-data.frame(pred_mean=solubles_area_jack_stat_pressed[,1],
                                        pred_low=solubles_area_jack_stat_pressed[,2],
                                        pred_high=solubles_area_jack_stat_pressed[,3],
                                        Measured=meta(pressed_spec_all_test)$solubles_area,
                                        ncomp=ncomp_solubles_area_pressed,
                                        Project=meta(pressed_spec_all_test)$Project,
                                        GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                        Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                        ID=meta(pressed_spec_all_test)$ID)

hemicellulose_area_jack_pred_pressed<-apply.coefs(hemicellulose_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
hemicellulose_area_jack_stat_pressed<-t(apply(hemicellulose_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_area_jack_df_pressed<-data.frame(pred_mean=hemicellulose_area_jack_stat_pressed[,1],
                                             pred_low=hemicellulose_area_jack_stat_pressed[,2],
                                             pred_high=hemicellulose_area_jack_stat_pressed[,3],
                                             Measured=meta(pressed_spec_all_test)$hemicellulose_area,
                                             ncomp=ncomp_hemicellulose_area_pressed,
                                             Project=meta(pressed_spec_all_test)$Project,
                                             GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                             Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                             ID=meta(pressed_spec_all_test)$ID)

cellulose_area_jack_pred_pressed<-apply.coefs(cellulose_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
cellulose_area_jack_stat_pressed<-t(apply(cellulose_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_area_jack_df_pressed<-data.frame(pred_mean=cellulose_area_jack_stat_pressed[,1],
                                         pred_low=cellulose_area_jack_stat_pressed[,2],
                                         pred_high=cellulose_area_jack_stat_pressed[,3],
                                         Measured=meta(pressed_spec_all_test)$cellulose_area,
                                         ncomp=ncomp_cellulose_area_pressed,
                                         Project=meta(pressed_spec_all_test)$Project,
                                         GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                         Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                         ID=meta(pressed_spec_all_test)$ID)

lignin_area_jack_pred_pressed<-apply.coefs(lignin_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
lignin_area_jack_stat_pressed<-t(apply(lignin_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_area_jack_df_pressed<-data.frame(pred_mean=lignin_area_jack_stat_pressed[,1],
                                      pred_low=lignin_area_jack_stat_pressed[,2],
                                      pred_high=lignin_area_jack_stat_pressed[,3],
                                      Measured=meta(pressed_spec_all_test)$lignin_area,
                                      ncomp=ncomp_lignin_area_pressed,
                                      Project=meta(pressed_spec_all_test)$Project,
                                      GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                      Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                      ID=meta(pressed_spec_all_test)$ID)

perC_area_jack_pred_pressed<-apply.coefs(perC_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
perC_area_jack_stat_pressed<-t(apply(perC_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_area_jack_df_pressed<-data.frame(pred_mean=perC_area_jack_stat_pressed[,1],
                                    pred_low=perC_area_jack_stat_pressed[,2],
                                    pred_high=perC_area_jack_stat_pressed[,3],
                                    Measured=meta(pressed_spec_all_test)$C_area,
                                    ncomp=ncomp_perC_area_pressed,
                                    Project=meta(pressed_spec_all_test)$Project,
                                    GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                    Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                    ID=meta(pressed_spec_all_test)$ID)

perN_area_jack_pred_pressed<-apply.coefs(perN_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
perN_area_jack_stat_pressed<-t(apply(perN_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_area_jack_df_pressed<-data.frame(pred_mean=perN_area_jack_stat_pressed[,1],
                                    pred_low=perN_area_jack_stat_pressed[,2],
                                    pred_high=perN_area_jack_stat_pressed[,3],
                                    Measured=meta(pressed_spec_all_test)$N_area,
                                    ncomp=ncomp_perN_area_pressed,
                                    Project=meta(pressed_spec_all_test)$Project,
                                    GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                    Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                    ID=meta(pressed_spec_all_test)$ID)

chlA_area_jack_pred_pressed<-apply.coefs(chlA_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
chlA_area_jack_stat_pressed<-t(apply(chlA_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_area_jack_df_pressed<-data.frame(pred_mean=chlA_area_jack_stat_pressed[,1],
                                    pred_low=chlA_area_jack_stat_pressed[,2],
                                    pred_high=chlA_area_jack_stat_pressed[,3],
                                    Measured=meta(pressed_spec_all_test)$chlA_area,
                                    ncomp=ncomp_chlA_area_pressed,
                                    Project=meta(pressed_spec_all_test)$Project,
                                    GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                    Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                    ID=meta(pressed_spec_all_test)$ID)

chlB_area_jack_pred_pressed<-apply.coefs(chlB_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
chlB_area_jack_stat_pressed<-t(apply(chlB_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_area_jack_df_pressed<-data.frame(pred_mean=chlB_area_jack_stat_pressed[,1],
                                    pred_low=chlB_area_jack_stat_pressed[,2],
                                    pred_high=chlB_area_jack_stat_pressed[,3],
                                    Measured=meta(pressed_spec_all_test)$chlB_area,
                                    ncomp=ncomp_chlB_area_pressed,
                                    Project=meta(pressed_spec_all_test)$Project,
                                    GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                    Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                    ID=meta(pressed_spec_all_test)$ID)

car_area_jack_pred_pressed<-apply.coefs(car_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
car_area_jack_stat_pressed<-t(apply(car_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_area_jack_df_pressed<-data.frame(pred_mean=car_area_jack_stat_pressed[,1],
                                   pred_low=car_area_jack_stat_pressed[,2],
                                   pred_high=car_area_jack_stat_pressed[,3],
                                   Measured=meta(pressed_spec_all_test)$car_area,
                                   ncomp=ncomp_car_area_pressed,
                                   Project=meta(pressed_spec_all_test)$Project,
                                   GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                   Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                   ID=meta(pressed_spec_all_test)$ID)

Al_area_jack_pred_pressed<-apply.coefs(Al_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Al_area_jack_stat_pressed<-t(apply(Al_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_area_jack_df_pressed<-data.frame(pred_mean=Al_area_jack_stat_pressed[,1],
                                  pred_low=Al_area_jack_stat_pressed[,2],
                                  pred_high=Al_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Al_area,
                                  ncomp=ncomp_Al_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

Ca_area_jack_pred_pressed<-apply.coefs(Ca_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Ca_area_jack_stat_pressed<-t(apply(Ca_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_area_jack_df_pressed<-data.frame(pred_mean=Ca_area_jack_stat_pressed[,1],
                                  pred_low=Ca_area_jack_stat_pressed[,2],
                                  pred_high=Ca_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Ca_area,
                                  ncomp=ncomp_Ca_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

Cu_area_jack_pred_pressed<-apply.coefs(Cu_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Cu_area_jack_stat_pressed<-t(apply(Cu_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_area_jack_df_pressed<-data.frame(pred_mean=Cu_area_jack_stat_pressed[,1],
                                  pred_low=Cu_area_jack_stat_pressed[,2],
                                  pred_high=Cu_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Cu_area,
                                  ncomp=ncomp_Cu_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

Fe_area_jack_pred_pressed<-apply.coefs(Fe_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Fe_area_jack_stat_pressed<-t(apply(Fe_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_area_jack_df_pressed<-data.frame(pred_mean=Fe_area_jack_stat_pressed[,1],
                                  pred_low=Fe_area_jack_stat_pressed[,2],
                                  pred_high=Fe_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Fe_area,
                                  ncomp=ncomp_Fe_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

K_area_jack_pred_pressed<-apply.coefs(K_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
K_area_jack_stat_pressed<-t(apply(K_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_area_jack_df_pressed<-data.frame(pred_mean=K_area_jack_stat_pressed[,1],
                                 pred_low=K_area_jack_stat_pressed[,2],
                                 pred_high=K_area_jack_stat_pressed[,3],
                                 Measured=meta(pressed_spec_all_test)$K_area,
                                 ncomp=ncomp_K_area_pressed,
                                 Project=meta(pressed_spec_all_test)$Project,
                                 GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                 ID=meta(pressed_spec_all_test)$ID)

Mg_area_jack_pred_pressed<-apply.coefs(Mg_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Mg_area_jack_stat_pressed<-t(apply(Mg_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_area_jack_df_pressed<-data.frame(pred_mean=Mg_area_jack_stat_pressed[,1],
                                  pred_low=Mg_area_jack_stat_pressed[,2],
                                  pred_high=Mg_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Mg_area,
                                  ncomp=ncomp_Mg_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

Mn_area_jack_pred_pressed<-apply.coefs(Mn_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Mn_area_jack_stat_pressed<-t(apply(Mn_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_area_jack_df_pressed<-data.frame(pred_mean=Mn_area_jack_stat_pressed[,1],
                                  pred_low=Mn_area_jack_stat_pressed[,2],
                                  pred_high=Mn_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Mn_area,
                                  ncomp=ncomp_Mn_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

Na_area_jack_pred_pressed<-apply.coefs(Na_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Na_area_jack_stat_pressed<-t(apply(Na_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_area_jack_df_pressed<-data.frame(pred_mean=Na_area_jack_stat_pressed[,1],
                                  pred_low=Na_area_jack_stat_pressed[,2],
                                  pred_high=Na_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Na_area,
                                  ncomp=ncomp_Na_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

P_area_jack_pred_pressed<-apply.coefs(P_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
P_area_jack_stat_pressed<-t(apply(P_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_area_jack_df_pressed<-data.frame(pred_mean=P_area_jack_stat_pressed[,1],
                                 pred_low=P_area_jack_stat_pressed[,2],
                                 pred_high=P_area_jack_stat_pressed[,3],
                                 Measured=meta(pressed_spec_all_test)$P_area,
                                 ncomp=ncomp_P_area_pressed,
                                 Project=meta(pressed_spec_all_test)$Project,
                                 GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                 Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                 ID=meta(pressed_spec_all_test)$ID)

Zn_area_jack_pred_pressed<-apply.coefs(Zn_area_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
Zn_area_jack_stat_pressed<-t(apply(Zn_area_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_area_jack_df_pressed<-data.frame(pred_mean=Zn_area_jack_stat_pressed[,1],
                                  pred_low=Zn_area_jack_stat_pressed[,2],
                                  pred_high=Zn_area_jack_stat_pressed[,3],
                                  Measured=meta(pressed_spec_all_test)$Zn_area,
                                  ncomp=ncomp_Zn_area_pressed,
                                  Project=meta(pressed_spec_all_test)$Project,
                                  GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                                  Discoloration=meta(pressed_spec_all_test)$Discoloration,
                                  ID=meta(pressed_spec_all_test)$ID)

##########################################################
## save output

pressed_area_jack_coef_list<-list(C=perC_area_jack_coefs_pressed,
                                N=perN_area_jack_coefs_pressed,
                                sol=solubles_area_jack_coefs_pressed,
                                hemi=hemicellulose_area_jack_coefs_pressed,
                                cell=cellulose_area_jack_coefs_pressed,
                                lign=lignin_area_jack_coefs_pressed,
                                chlA=chlA_area_jack_coefs_pressed,
                                chlB=chlB_area_jack_coefs_pressed,
                                car=car_area_jack_coefs_pressed,
                                Al=Al_area_jack_coefs_pressed,
                                Ca=Ca_area_jack_coefs_pressed,
                                Cu=Cu_area_jack_coefs_pressed,
                                Fe=Fe_area_jack_coefs_pressed,
                                K=K_area_jack_coefs_pressed,
                                Mg=Mg_area_jack_coefs_pressed,
                                Mn=Mn_area_jack_coefs_pressed,
                                Na=Na_area_jack_coefs_pressed,
                                P=P_area_jack_coefs_pressed,
                                Zn=Zn_area_jack_coefs_pressed)
saveRDS(pressed_area_jack_coef_list,"SavedResults/pressed_area_jack_coefs_list.rds")

pressed_area_jack_df_list<-list(C=perC_area_jack_df_pressed,
                              N=perN_area_jack_df_pressed,
                              sol=solubles_area_jack_df_pressed,
                              hemi=hemicellulose_area_jack_df_pressed,
                              cell=cellulose_area_jack_df_pressed,
                              lign=lignin_area_jack_df_pressed,
                              chlA=chlA_area_jack_df_pressed,
                              chlB=chlB_area_jack_df_pressed,
                              car=car_area_jack_df_pressed,
                              Al=Al_area_jack_df_pressed,
                              Ca=Ca_area_jack_df_pressed,
                              Cu=Cu_area_jack_df_pressed,
                              Fe=Fe_area_jack_df_pressed,
                              K=K_area_jack_df_pressed,
                              Mg=Mg_area_jack_df_pressed,
                              Mn=Mn_area_jack_df_pressed,
                              Na=Na_area_jack_df_pressed,
                              P=P_area_jack_df_pressed,
                              Zn=Zn_area_jack_df_pressed)
saveRDS(pressed_area_jack_df_list,"SavedResults/pressed_area_jack_df_list.rds")
