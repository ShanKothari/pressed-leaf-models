setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")

library(spectrolab)
library(ggplot2)
library(RColorBrewer)
library(pls)
library(reshape2)
library(ggpubr)
source("Scripts/pressed-leaf-models/00 useful_functions.R")

###################################################
## the purpose of this section is to add on the corrected
## values of EWT to the fresh-leaf spectra
## without creating a new training/testing data split

## script 03 has been modified so that the trait data
## incorporates the corrected values of EWT in the file
## fresh_spec_all.rds
## we'll use those values to create a new column in the
## training and testing splits

## read file
fresh_spec_all<-readRDS("ProcessedSpectralData/fresh_spec_all.rds")
pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_all.rds")
ground_spec_all<-readRDS("ProcessedSpectralData/ground_spec_all.rds")

## optionally, we could 'fill in' EWT (fresh leaves) for Dessain
## with EWT_rehydrated since no fresh mass was recorded
## however, for the correction won't!
## we include only the other projects

# meta(fresh_spec_all)$EWT[meta(fresh_spec_all)$Project=="Dessain"]<-
# meta(fresh_spec_all)$EWT_rehydrated[meta(fresh_spec_all)$Project=="Dessain"]

## read original train/test split
fresh_spec_all_train<-readRDS("ProcessedSpectralData/fresh_spec_all_train.rds")
fresh_spec_all_test<-readRDS("ProcessedSpectralData/fresh_spec_all_test.rds")

pressed_spec_all_train<-readRDS("ProcessedSpectralData/pressed_spec_all_train.rds")
pressed_spec_all_test<-readRDS("ProcessedSpectralData/pressed_spec_all_test.rds")

ground_spec_all_train<-readRDS("ProcessedSpectralData/ground_spec_all_train.rds")
ground_spec_all_test<-readRDS("ProcessedSpectralData/ground_spec_all_test.rds")

## in that original split, EWT_rehydrated is simply called EWT
## so we add a new column for true EWT called EWT_actual
meta(fresh_spec_all_train)$EWT_actual<-meta(fresh_spec_all)$EWT[match(meta(fresh_spec_all_train)$ID,
                                                                      meta(fresh_spec_all)$ID)]
meta(fresh_spec_all_test)$EWT_actual<-meta(fresh_spec_all)$EWT[match(meta(fresh_spec_all_test)$ID,
                                                                     meta(fresh_spec_all)$ID)]

meta(pressed_spec_all_train)$EWT_actual<-meta(pressed_spec_all)$EWT[match(meta(pressed_spec_all_train)$ID,
                                                                          meta(pressed_spec_all)$ID)]
meta(pressed_spec_all_test)$EWT_actual<-meta(pressed_spec_all)$EWT[match(meta(pressed_spec_all_test)$ID,
                                                                         meta(pressed_spec_all)$ID)]

meta(ground_spec_all_train)$EWT_actual<-meta(ground_spec_all)$EWT[match(meta(ground_spec_all_train)$ID,
                                                                        meta(ground_spec_all)$ID)]
meta(ground_spec_all_test)$EWT_actual<-meta(ground_spec_all)$EWT[match(meta(ground_spec_all_test)$ID,
                                                                       meta(ground_spec_all)$ID)]

##################################
## here we develop new models for fresh EWT using the corrected data
## and fresh-leaf spectra

EWT_fresh<-plsr(meta(fresh_spec_all_train)$EWT_actual~as.matrix(fresh_spec_all_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_fresh <- selectNcomp(EWT_fresh, method = "onesigma", plot = FALSE)
EWT_fresh_valid <- which(!is.na(meta(fresh_spec_all_train)$EWT_actual))
EWT_fresh_pred<-data.frame(ID=meta(fresh_spec_all_train)$ID[EWT_fresh_valid],
                           Species=meta(fresh_spec_all_train)$Species[EWT_fresh_valid],
                           Project=meta(fresh_spec_all_train)$Project[EWT_fresh_valid],
                           Stage=meta(fresh_spec_all_train)$Stage[EWT_fresh_valid],
                           GrowthForm=meta(fresh_spec_all_train)$GrowthForm[EWT_fresh_valid],
                           measured=meta(fresh_spec_all_train)$EWT_actual[EWT_fresh_valid],
                           val_pred=EWT_fresh$validation$pred[,,ncomp_EWT_fresh])

ggplot(EWT_fresh_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting EWT from fresh-leaf spectra")

EWT_jack_coefs_fresh<-list()
EWT_jack_stats_fresh<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_fresh<-nrow(fresh_spec_all_train)
  train_jack_fresh<-sample(1:n_cal_spec_fresh,floor(0.7*n_cal_spec_fresh))
  test_jack_fresh<-setdiff(1:n_cal_spec_fresh,train_jack_fresh)
  
  calib_jack_fresh<-fresh_spec_all_train[train_jack_fresh]
  val_jack_fresh<-fresh_spec_all_train[test_jack_fresh]
  
  EWT_fresh_jack<-plsr(meta(calib_jack_fresh)$EWT~as.matrix(calib_jack_fresh),
                       ncomp=30,method = "oscorespls",validation="none")

  EWT_jack_val_pred_fresh<-as.vector(predict(EWT_fresh_jack,newdata=as.matrix(val_jack_fresh),ncomp=ncomp_EWT_fresh)[,,1])
  EWT_jack_val_fit_fresh<-lm(meta(val_jack_fresh)$EWT~EWT_jack_val_pred_fresh)
  EWT_jack_stats_fresh[[i]]<-c(R2=summary(EWT_jack_val_fit_fresh)$r.squared,
                               RMSE=RMSD(meta(val_jack_fresh)$EWT,EWT_jack_val_pred_fresh),
                               perRMSE=percentRMSD(meta(val_jack_fresh)$EWT,EWT_jack_val_pred_fresh,0.025,0.975),
                               bias=mean(EWT_jack_val_pred_fresh,na.rm=T)-mean(meta(val_jack_fresh)$EWT,na.rm=T))
  
  EWT_jack_coefs_fresh[[i]]<-as.vector(coef(EWT_fresh_jack,ncomp=ncomp_EWT_fresh,intercept=TRUE))
}

EWT_jack_pred_fresh<-apply.coefs(EWT_jack_coefs_fresh,as.matrix(fresh_spec_all_test))
EWT_jack_stat_fresh<-t(apply(EWT_jack_pred_fresh,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_jack_df_fresh<-data.frame(pred_mean=EWT_jack_stat_fresh[,1],
                              pred_low=EWT_jack_stat_fresh[,2],
                              pred_high=EWT_jack_stat_fresh[,3],
                              Measured=meta(fresh_spec_all_test)$EWT_actual,
                              ncomp=ncomp_EWT_fresh,
                              Project=meta(fresh_spec_all_test)$Project,
                              GrowthForm=meta(fresh_spec_all_test)$GrowthForm,
                              Discoloration=meta(fresh_spec_all_test)$Discoloration,
                              ID=meta(fresh_spec_all_test)$ID)

############################################
## here we develop new models for fresh EWT using the corrected data
## and pressed-leaf spectra

EWT_pressed<-plsr(meta(pressed_spec_all_train)$EWT_actual~as.matrix(pressed_spec_all_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_pressed <- selectNcomp(EWT_pressed, method = "onesigma", plot = FALSE)
EWT_pressed_valid <- which(!is.na(meta(pressed_spec_all_train)$EWT_actual))
EWT_pressed_pred<-data.frame(ID=meta(pressed_spec_all_train)$ID[EWT_pressed_valid],
                             Species=meta(pressed_spec_all_train)$Species[EWT_pressed_valid],
                             Project=meta(pressed_spec_all_train)$Project[EWT_pressed_valid],
                             Stage=meta(pressed_spec_all_train)$Stage[EWT_pressed_valid],
                             GrowthForm=meta(pressed_spec_all_train)$GrowthForm[EWT_pressed_valid],
                             measured=meta(pressed_spec_all_train)$EWT_actual[EWT_pressed_valid],
                             val_pred=EWT_pressed$validation$pred[,,ncomp_EWT_pressed])

ggplot(EWT_pressed_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting EWT from pressed-leaf spectra")

EWT_jack_coefs_pressed<-list()
EWT_jack_stats_pressed<-list()

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_pressed<-nrow(pressed_spec_all_train)
  train_jack_pressed<-sample(1:n_cal_spec_pressed,floor(0.7*n_cal_spec_pressed))
  test_jack_pressed<-setdiff(1:n_cal_spec_pressed,train_jack_pressed)
  
  calib_jack_pressed<-pressed_spec_all_train[train_jack_pressed]
  val_jack_pressed<-pressed_spec_all_train[test_jack_pressed]
  
  EWT_pressed_jack<-plsr(meta(calib_jack_pressed)$EWT~as.matrix(calib_jack_pressed),
                       ncomp=30,method = "oscorespls",validation="none")
  
  EWT_jack_val_pred_pressed<-as.vector(predict(EWT_pressed_jack,newdata=as.matrix(val_jack_pressed),ncomp=ncomp_EWT_pressed)[,,1])
  EWT_jack_val_fit_pressed<-lm(meta(val_jack_pressed)$EWT~EWT_jack_val_pred_pressed)
  EWT_jack_stats_pressed[[i]]<-c(R2=summary(EWT_jack_val_fit_pressed)$r.squared,
                               RMSE=RMSD(meta(val_jack_pressed)$EWT,EWT_jack_val_pred_pressed),
                               perRMSE=percentRMSD(meta(val_jack_pressed)$EWT,EWT_jack_val_pred_pressed,0.025,0.975),
                               bias=mean(EWT_jack_val_pred_pressed,na.rm=T)-mean(meta(val_jack_pressed)$EWT,na.rm=T))
  
  EWT_jack_coefs_pressed[[i]]<-as.vector(coef(EWT_pressed_jack,ncomp=ncomp_EWT_pressed,intercept=TRUE))
}

EWT_jack_pred_pressed<-apply.coefs(EWT_jack_coefs_pressed,as.matrix(pressed_spec_all_test))
EWT_jack_stat_pressed<-t(apply(EWT_jack_pred_pressed,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_jack_df_pressed<-data.frame(pred_mean=EWT_jack_stat_pressed[,1],
                              pred_low=EWT_jack_stat_pressed[,2],
                              pred_high=EWT_jack_stat_pressed[,3],
                              Measured=meta(pressed_spec_all_test)$EWT_actual,
                              ncomp=ncomp_EWT_pressed,
                              Project=meta(pressed_spec_all_test)$Project,
                              GrowthForm=meta(pressed_spec_all_test)$GrowthForm,
                              Discoloration=meta(pressed_spec_all_test)$Discoloration,
                              ID=meta(pressed_spec_all_test)$ID)

#######################################
## here we develop new models for fresh EWT using the corrected data
## and ground-leaf spectra

EWT_jack_coefs_ground<-list()
EWT_jack_stats_ground<-list()

EWT_ground<-plsr(meta(ground_spec_all_train)$EWT_actual~as.matrix(ground_spec_all_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_ground <- selectNcomp(EWT_ground, method = "onesigma", plot = FALSE)
EWT_ground_valid <- which(!is.na(meta(ground_spec_all_train)$EWT_actual))
EWT_ground_pred<-data.frame(ID=meta(ground_spec_all_train)$ID[EWT_ground_valid],
                            Species=meta(ground_spec_all_train)$Species[EWT_ground_valid],
                            Project=meta(ground_spec_all_train)$Project[EWT_ground_valid],
                            Stage=meta(ground_spec_all_train)$Stage[EWT_ground_valid],
                            GrowthForm=meta(ground_spec_all_train)$GrowthForm[EWT_ground_valid],
                            measured=meta(ground_spec_all_train)$EWT_actual[EWT_ground_valid],
                            val_pred=EWT_ground$validation$pred[,,ncomp_EWT_ground])

ggplot(EWT_ground_pred,aes(y=measured,x=val_pred,color=GrowthForm))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Predicted",y="Measured")+
  ggtitle("Predicting EWT from ground-leaf spectra")

for(i in 1:nreps){
  print(i)
  
  n_cal_spec_ground<-nrow(ground_spec_all_train)
  train_jack_ground<-sample(1:n_cal_spec_ground,floor(0.7*n_cal_spec_ground))
  test_jack_ground<-setdiff(1:n_cal_spec_ground,train_jack_ground)
  
  calib_jack_ground<-ground_spec_all_train[train_jack_ground]
  val_jack_ground<-ground_spec_all_train[test_jack_ground]
  
  EWT_ground_jack<-plsr(meta(calib_jack_ground)$EWT~as.matrix(calib_jack_ground),
                       ncomp=30,method = "oscorespls",validation="none")
  
  EWT_jack_val_pred_ground<-as.vector(predict(EWT_ground_jack,newdata=as.matrix(val_jack_ground),ncomp=ncomp_EWT_ground)[,,1])
  EWT_jack_val_fit_ground<-lm(meta(val_jack_ground)$EWT~EWT_jack_val_pred_ground)
  EWT_jack_stats_ground[[i]]<-c(R2=summary(EWT_jack_val_fit_ground)$r.squared,
                               RMSE=RMSD(meta(val_jack_ground)$EWT,EWT_jack_val_pred_ground),
                               perRMSE=percentRMSD(meta(val_jack_ground)$EWT,EWT_jack_val_pred_ground,0.025,0.975),
                               bias=mean(EWT_jack_val_pred_ground,na.rm=T)-mean(meta(val_jack_ground)$EWT,na.rm=T))
  
  EWT_jack_coefs_ground[[i]]<-as.vector(coef(EWT_ground_jack,ncomp=ncomp_EWT_ground,intercept=TRUE))
}

EWT_jack_pred_ground<-apply.coefs(EWT_jack_coefs_ground,as.matrix(ground_spec_all_test))
EWT_jack_stat_ground<-t(apply(EWT_jack_pred_ground,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_jack_df_ground<-data.frame(pred_mean=EWT_jack_stat_ground[,1],
                              pred_low=EWT_jack_stat_ground[,2],
                              pred_high=EWT_jack_stat_ground[,3],
                              Measured=meta(ground_spec_all_test)$EWT_actual,
                              ncomp=ncomp_EWT_ground,
                              Project=meta(ground_spec_all_test)$Project,
                              GrowthForm=meta(ground_spec_all_test)$GrowthForm,
                              Discoloration=meta(ground_spec_all_test)$Discoloration,
                              ID=meta(ground_spec_all_test)$ID)

############################################
## save outputs

saveRDS(EWT_jack_coefs_fresh,"SavedResults/EWT_corrected_jack_coefs_fresh.rds")
saveRDS(EWT_jack_df_fresh,"SavedResults/EWT_corrected_jack_df_fresh.rds")
saveRDS(EWT_jack_stats_fresh,"SavedResults/EWT_corrected_jack_stats_fresh.rds")

saveRDS(EWT_jack_coefs_pressed,"SavedResults/EWT_corrected_jack_coefs_pressed.rds")
saveRDS(EWT_jack_df_pressed,"SavedResults/EWT_corrected_jack_df_pressed.rds")
saveRDS(EWT_jack_stats_pressed,"SavedResults/EWT_corrected_jack_stats_pressed.rds")

saveRDS(EWT_jack_coefs_ground,"SavedResults/EWT_corrected_jack_coefs_ground.rds")
saveRDS(EWT_jack_df_ground,"SavedResults/EWT_corrected_jack_df_ground.rds")
saveRDS(EWT_jack_stats_ground,"SavedResults/EWT_corrected_jack_stats_ground.rds")

## save coefficients
write.coefs<-function(obj,path,filename){
  coef_mat<-matrix(unlist(obj),nrow=length(obj),byrow=T)
  colnames(coef_mat)<-c("intercept",400:2400)
  write.csv(coef_mat,
            paste(path,filename,".csv",sep=""),
            row.names=F)
}

write.coefs(obj=EWT_jack_coefs_fresh,
            path="ModelCoefficients/EWTCorrectedModels/",
            filename="EWTActual_fresh")

write.coefs(obj=EWT_jack_coefs_pressed,
            path="ModelCoefficients/EWTCorrectedModels/",
            filename="EWTActual_pressed")

write.coefs(obj=EWT_jack_coefs_ground,
            path="ModelCoefficients/EWTCorrectedModels/",
            filename="EWTActual_ground")

###############################################
## plotting internal validation output

focal_palette=palette(brewer.pal(8,name="Set2")[c(3,4,5,6,8,1,2)])

all.EWT<-c(EWT_jack_df_pressed$Measured,
           EWT_jack_df_fresh$pred.mean,
           EWT_jack_df_pressed$pred.mean,
           EWT_jack_df_ground$pred.mean)
EWT_upper<-max(all.EWT,na.rm=T)+0.02
EWT_lower<-min(all.EWT,na.rm=T)-0.02

EWT_actual_fresh_plot<-ggplot(EWT_jack_df_fresh,
                                aes(y=Measured,x=pred_mean,color=GrowthForm))+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color="none")+
  scale_color_manual(values=focal_palette)+
  ggtitle("Fresh-leaf spectra")

EWT_actual_pressed_plot<-ggplot(EWT_jack_df_pressed,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color="none")+
  scale_color_manual(values=focal_palette)+
  ggtitle("Pressed-leaf spectra")

EWT_actual_ground_plot<-ggplot(EWT_jack_df_ground,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)",
       color="Growth form")+
  scale_color_manual(values=focal_palette)+
  ggtitle("Ground-leaf spectra")

pdf("Manuscript/EWT_corrected_val_plot.pdf",width=16,height=7)
(EWT_actual_fresh_plot + EWT_actual_pressed_plot + EWT_actual_ground_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

summary(lm(Measured~pred_mean,data=EWT_jack_df_fresh))
with(EWT_jack_df_fresh,
     RMSD(measured = Measured,predicted = pred_mean))
with(EWT_jack_df_fresh,
     percentRMSD(measured = Measured,predicted = pred_mean,
                 min=0.025,max=0.975))

########################################
## external validation

pressed_spec_MN_RWC_agg<-read.csv("ProcessedSpectralData/pressed_spec_MN_RWC_avg.csv")
colnames(pressed_spec_MN_RWC_agg)<-gsub("X","",colnames(pressed_spec_MN_RWC_agg))
pressed_spec_MN_RWC_spec<-as_spectra(pressed_spec_MN_RWC_agg,
                                     name_idx = 1,
                                     meta_idxs = 2:19)
meta(pressed_spec_MN_RWC_spec)$EWT<-as.numeric(meta(pressed_spec_MN_RWC_spec)$EWT)

EWT_ext_pressed<-apply.coefs(EWT_jack_coefs_pressed,
                             val.spec = pressed_spec_MN_RWC_spec)
EWT_ext_pressed_stat<-t(apply(EWT_ext_pressed,1,
                              function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_ext_pressed_pred_df<-data.frame(Measured=meta(pressed_spec_MN_RWC_spec)$EWT,
                                    pred_mean=EWT_ext_pressed_stat[,1],
                                    pred_low=EWT_ext_pressed_stat[,2],
                                    pred_high=EWT_ext_pressed_stat[,3],
                                    FunctionalGroup=meta(pressed_spec_MN_RWC_spec)$FunctionalGroup)
EWT_all<-with(EWT_ext_pressed_pred_df,c(pred_low[!is.na(Measured)],
                                        pred_high[!is.na(Measured)],
                                        Measured))
EWT_upper<-max(EWT_all,na.rm=T)+0.03
EWT_lower<-min(EWT_all,na.rm=T)-0.03

EWT_ind_val<-ggplot(data=EWT_ext_pressed_pred_df,
                    aes(x=pred_mean,y=Measured,color=FunctionalGroup))+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  scale_color_manual(values=focal_palette)+
  guides(color=guide_legend(title="Functional group"))

## to dos:
## output the data for upload to EcoSIS