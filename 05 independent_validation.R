setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(ggplot2)
library(dplyr)
library(signal)
library(patchwork)
library(stringr)

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

#####################################
## process data

pressed_spec_MN<-read_spectra("PressedSpectra/JCBLab","sed",exclude_if_matches = c("BAD","WR_"))
pressed_spec_MN<-match_sensors(pressed_spec_MN,splice_at=983)

pressed_spec_MN_spl<-strsplit(names(pressed_spec_MN),split = "_")
pressed_spec_MN_ID<-lapply(pressed_spec_MN_spl,function(x) x[-length(x)])
meta(pressed_spec_MN)$Plot<-unlist(lapply(pressed_spec_MN_ID,function(x) x[[1]]))
meta(pressed_spec_MN)$Species<-unlist(lapply(pressed_spec_MN_ID,function(x) x[[2]]))
meta(pressed_spec_MN)$Location<-unlist(lapply(pressed_spec_MN_ID,function(x) x[[3]]))
meta(pressed_spec_MN)$Type<-unlist(lapply(pressed_spec_MN_ID,function(x) x[[length(x)]]))
meta(pressed_spec_MN)$ID<-toupper(unlist(lapply(pressed_spec_MN_ID,
                                                function(x) paste(x[-length(x)],collapse="_"))))
splookup<-read.csv("Traits/JCBdata/splookup.csv")
meta(pressed_spec_MN)$functional.group<-splookup$Group[match(meta(pressed_spec_MN)$Species,splookup$SpeciesCode)]
## remove ACNE due to discrepancy about including the rachis or not
pressed_spec_MN<-pressed_spec_MN[-which(meta(pressed_spec_MN)$Species=="ACNE")]

## only spectra for RWC
pressed_spec_MN_RWC<-pressed_spec_MN[-grep("PIG",names(pressed_spec_MN))]

## aggregate spectra with the same ID
pressed_spec_MN_agg<-aggregate(pressed_spec_MN,by=meta(pressed_spec_MN)$ID,
                               FUN=try_keep_txt(mean))
pressed_spec_MN_RWC_agg<-aggregate(pressed_spec_MN_RWC,by=meta(pressed_spec_MN_RWC)$ID,
                               FUN=try_keep_txt(mean))

pressed_spec_MN_agg <- pressed_spec_MN_agg[,400:2400]
pressed_spec_MN_RWC_agg <- pressed_spec_MN_RWC_agg[,400:2400]

pressed_spec_MN_agg<-pressed_spec_MN_agg[-which(meta(pressed_spec_MN_agg)$ID=="31_BAPLE_1")]
pressed_spec_MN_RWC_agg<-pressed_spec_MN_RWC_agg[-which(meta(pressed_spec_MN_RWC_agg)$ID=="31_BAPLE_1")]

##################################
## read traits

freshmass<-read.csv("Traits/JCBdata/Relative Water Content -- Trait-Spectra Models 2018 - Sheet1.csv")
freshmass$ID<-toupper(apply(freshmass,1,function(x) paste(x[1:4],collapse="_")))
freshmass$ID<-sub("_$","",freshmass$ID)
freshmass$Wet.mass.petiole.rachis..g.[is.na(freshmass$Wet.mass.petiole.rachis..g.)]<-0
freshmass$lamina<-freshmass$Total.wet.mass..g.-freshmass$Wet.mass.petiole.rachis..g.

drymass<-read.csv("Traits/JCBdata/Dry Mass -- Trait-Spectra Models 2018 - Sheet1.csv")
## adjust this to include more species
drymass_RWC<-drymass[drymass$Purpose=="RWC",]
drymass_RWC$ID<-toupper(apply(drymass_RWC,1,function(x) paste(x[1:4],collapse="_")))
drymass_RWC$ID<-sub("_$","",drymass_RWC$ID)
drymass_RWC$Dry.mass.petiole.rachis..g.[is.na(drymass_RWC$Dry.mass.petiole.rachis..g.)]<-0
drymass_RWC$lamina<-drymass_RWC$Total.dry.mass..g.-drymass_RWC$Dry.mass.petiole.rachis..g.

meta(pressed_spec_MN_RWC_agg)$fresh_mass<-freshmass$lamina[match(meta(pressed_spec_MN_RWC_agg)$ID,freshmass$ID)]
## some missing matches for dry mass?
meta(pressed_spec_MN_RWC_agg)$dry_mass<-drymass_RWC$lamina[match(meta(pressed_spec_MN_RWC_agg)$ID,drymass_RWC$ID)]
meta(pressed_spec_MN_RWC_agg)$LDMC<-with(meta(pressed_spec_MN_RWC_agg),dry_mass/fresh_mass)*1000

## attach dry mass also to full (not just RWC) data

area_all<-read.csv("Traits/JCBdata/LeafArea/area_all.csv")
meta(pressed_spec_MN_RWC_agg)$area<-area_all$area[match(toupper(meta(pressed_spec_MN_RWC_agg)$ID),toupper(area_all$TrueID))]
meta(pressed_spec_MN_RWC_agg)$LMA<-with(meta(pressed_spec_MN_RWC_agg),dry_mass/area)*10
meta(pressed_spec_MN_RWC_agg)$EWT<-with(meta(pressed_spec_MN_RWC_agg),(1/(LDMC/1000)-1)*LMA)/10

CN<-read.csv("Traits/JCBdata/ElementalAnalysis/SummaryTables/fullsummarytable.csv")
CN$SampleID<-gsub(pattern=" ",replacement="_",x = CN$SampleID)
meta(pressed_spec_MN_agg)$C<-CN$Carbon[match(meta(pressed_spec_MN_agg)$ID,CN$SampleID)]
meta(pressed_spec_MN_agg)$N<-CN$Nitrogen[match(meta(pressed_spec_MN_agg)$ID,CN$SampleID)]
meta(pressed_spec_MN_agg)$EA_run<-CN$Run[match(meta(pressed_spec_MN_agg)$ID,CN$SampleID)]

#####################################
## read coefficients

pressed_jack_coef_list<-readRDS("SavedResults/pressed_jack_coefs_list.rds")
pressed_1300_jack_coef_list<-readRDS("SavedResults/pressed_1300_jack_coefs_list.rds")

LDMC_preds<-apply.coefs(pressed_jack_coef_list$LDMC,as.matrix(pressed_spec_MN_RWC_agg))
meta(pressed_spec_MN_RWC_agg)$LDMC_pred<-rowMeans(LDMC_preds)
meta(pressed_spec_MN_RWC_agg)$LDMC_lower<-apply(LDMC_preds,1,quantile,probs=0.025)
meta(pressed_spec_MN_RWC_agg)$LDMC_upper<-apply(LDMC_preds,1,quantile,probs=0.975)

LDMC_preds_1300<-apply.coefs(pressed_1300_jack_coef_list$LDMC,as.matrix(pressed_spec_MN_RWC_agg[,1300:2400]))
meta(pressed_spec_MN_RWC_agg)$LDMC_pred_1300<-rowMeans(LDMC_preds_1300)
meta(pressed_spec_MN_RWC_agg)$LDMC_1300_lower<-apply(LDMC_preds_1300,1,quantile,probs=0.025)
meta(pressed_spec_MN_RWC_agg)$LDMC_1300_upper<-apply(LDMC_preds_1300,1,quantile,probs=0.975)

LMA_preds<-apply.coefs(pressed_jack_coef_list$LMA,as.matrix(pressed_spec_MN_RWC_agg))
meta(pressed_spec_MN_RWC_agg)$LMA_pred<-rowMeans(LMA_preds)
meta(pressed_spec_MN_RWC_agg)$LMA_lower<-apply(LMA_preds,1,quantile,probs=0.025)
meta(pressed_spec_MN_RWC_agg)$LMA_upper<-apply(LMA_preds,1,quantile,probs=0.975)

LMA_preds_1300<-apply.coefs(pressed_1300_jack_coef_list$LMA,as.matrix(pressed_spec_MN_RWC_agg[,1300:2400]))
meta(pressed_spec_MN_RWC_agg)$LMA_pred_1300<-rowMeans(LMA_preds_1300)
meta(pressed_spec_MN_RWC_agg)$LMA_1300_lower<-apply(LMA_preds_1300,1,quantile,probs=0.025)
meta(pressed_spec_MN_RWC_agg)$LMA_1300_upper<-apply(LMA_preds_1300,1,quantile,probs=0.975)

EWT_preds<-apply.coefs(pressed_jack_coef_list$EWT,as.matrix(pressed_spec_MN_RWC_agg))
meta(pressed_spec_MN_RWC_agg)$EWT_pred<-rowMeans(EWT_preds)
meta(pressed_spec_MN_RWC_agg)$EWT_lower<-apply(EWT_preds,1,quantile,probs=0.025)
meta(pressed_spec_MN_RWC_agg)$EWT_upper<-apply(EWT_preds,1,quantile,probs=0.975)

EWT_preds_1300<-apply.coefs(pressed_1300_jack_coef_list$EWT,as.matrix(pressed_spec_MN_RWC_agg[,1300:2400]))
meta(pressed_spec_MN_RWC_agg)$EWT_pred_1300<-rowMeans(EWT_preds_1300)
meta(pressed_spec_MN_RWC_agg)$EWT_1300_lower<-apply(EWT_preds_1300,1,quantile,probs=0.025)
meta(pressed_spec_MN_RWC_agg)$EWT_1300_upper<-apply(EWT_preds_1300,1,quantile,probs=0.975)

C_preds<-apply.coefs(pressed_jack_coef_list$perC,as.matrix(pressed_spec_MN_agg))
meta(pressed_spec_MN_agg)$C_pred<-rowMeans(C_preds)
meta(pressed_spec_MN_agg)$C_lower<-apply(C_preds,1,quantile,probs=0.025)
meta(pressed_spec_MN_agg)$C_upper<-apply(C_preds,1,quantile,probs=0.975)

N_preds<-apply.coefs(pressed_jack_coef_list$perN,as.matrix(pressed_spec_MN_agg))
meta(pressed_spec_MN_agg)$N_pred<-rowMeans(N_preds)
meta(pressed_spec_MN_agg)$N_lower<-apply(N_preds,1,quantile,probs=0.025)
meta(pressed_spec_MN_agg)$N_upper<-apply(N_preds,1,quantile,probs=0.975)

C_preds_1300<-apply.coefs(pressed_1300_jack_coef_list$perC,as.matrix(pressed_spec_MN_agg[,1300:2400]))
meta(pressed_spec_MN_agg)$C_pred_1300<-rowMeans(C_preds_1300)
meta(pressed_spec_MN_agg)$C_1300_lower<-apply(C_preds_1300,1,quantile,probs=0.025)
meta(pressed_spec_MN_agg)$C_1300_upper<-apply(C_preds_1300,1,quantile,probs=0.975)

N_preds_1300<-apply.coefs(pressed_1300_jack_coef_list$perN,as.matrix(pressed_spec_MN_agg[,1300:2400]))
meta(pressed_spec_MN_agg)$N_pred_1300<-rowMeans(N_preds_1300)
meta(pressed_spec_MN_agg)$N_1300_lower<-apply(N_preds_1300,1,quantile,probs=0.025)
meta(pressed_spec_MN_agg)$N_1300_upper<-apply(N_preds_1300,1,quantile,probs=0.975)

N_indval_plot<-ggplot(meta(pressed_spec_MN_agg),
       aes(x=N_pred,y=N,color=functional.group))+
  geom_errorbarh(aes(y=N,xmin=N_lower,xmax=N_upper),
               color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(-1,5),ylim=c(-1,5))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Predicted N (%)",y="Measured N (%)")+
  guides(color=guide_legend(title="Functional group"))

N_1300_indval_plot<-ggplot(meta(pressed_spec_MN_agg),
       aes(x=N_pred_1300,y=N,color=functional.group))+
  geom_errorbarh(aes(y=N,xmin=N_1300_lower,xmax=N_1300_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(-1,5),ylim=c(-1,5))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Predicted N (%)",y="Measured N (%)")+
  guides(color=guide_legend(title="Functional group"))

LMA_indval_plot<-ggplot(meta(pressed_spec_MN_RWC_agg),
       aes(x=LMA_pred,y=LMA,color=functional.group))+
  geom_errorbarh(aes(y=LMA,xmin=LMA_lower,xmax=LMA_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.25))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=guide_legend(title="Functional group"))

LMA_1300_indval_plot<-ggplot(meta(pressed_spec_MN_RWC_agg),
                             aes(x=LMA_pred_1300,y=LMA,color=functional.group))+
  geom_errorbarh(aes(y=LMA,xmin=LMA_1300_lower,xmax=LMA_1300_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.25))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=guide_legend(title="Functional group"))

EWT_indval_plot<-ggplot(meta(pressed_spec_MN_RWC_agg),
       aes(x=EWT_pred,y=EWT,color=functional.group))+
  geom_errorbarh(aes(y=EWT,xmin=EWT_lower,xmax=EWT_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(0,0.04),ylim=c(0,0.04))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Predicted EWT (cm)",y="Measured EWT (cm)")+
  guides(color=guide_legend(title="Functional group"))

EWT_1300_indval_plot<-ggplot(meta(pressed_spec_MN_RWC_agg),
       aes(x=EWT_pred_1300,y=EWT,color=functional.group))+
  geom_errorbarh(aes(y=EWT,xmin=EWT_1300_lower,xmax=EWT_1300_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(0,0.04),ylim=c(0,0.04))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Predicted EWT (cm)",y="Measured EWT (cm)")+
  guides(color=guide_legend(title="Functional group"))

C_indval_plot<-ggplot(meta(pressed_spec_MN_agg),
       aes(x=C_pred,y=C,color=functional.group))+
  geom_errorbarh(aes(y=C,xmin=C_lower,xmax=C_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(35,60),ylim=c(35,60))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Predicted C (%)",y="Measured C (%)")+
  guides(color=guide_legend(title="Functional group"))

C_1300_indval_plot<-ggplot(meta(pressed_spec_MN_agg),
       aes(x=C_pred_1300,y=C,color=functional.group))+
  geom_errorbarh(aes(y=C,xmin=C_1300_lower,xmax=C_1300_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(35,60),ylim=c(35,60))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Predicted C (%)",y="Measured C (%)")+
  guides(color=guide_legend(title="Functional group"))

LDMC_indval_plot<-ggplot(meta(pressed_spec_MN_RWC_agg),
       aes(x=LDMC_pred,y=LDMC,color=functional.group))+
  geom_errorbarh(aes(y=LDMC,xmin=LDMC_lower,xmax=LDMC_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(100,600),ylim=c(100,600))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression("Predicted LDMC (mg g"^-1*")"),
       y=expression("Measured LDMC (mg g"^-1*")"))+
  guides(color=guide_legend(title="Functional group"))

LDMC_1300_indval_plot<-ggplot(meta(pressed_spec_MN_RWC_agg),
       aes(x=LDMC_pred_1300,y=LDMC,color=functional.group))+
  geom_errorbarh(aes(y=LDMC,xmin=LDMC_1300_lower,xmax=LDMC_1300_upper),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  coord_cartesian(xlim=c(100,700),ylim=c(100,700))+
  geom_abline(slope=1,intercept=0,size=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression("Predicted LDMC (mg g"^-1*")"),
       y=expression("Measured LDMC (mg g"^-1*")"))+
  guides(color=guide_legend(title="Functional group"))

pdf("Manuscript/Fig5.pdf",width=9,height=12)
((LMA_indval_plot/LDMC_indval_plot/EWT_indval_plot)|
  (N_indval_plot/C_indval_plot/plot_spacer()))+
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')
dev.off()

pdf("Manuscript/FigS2.pdf",width=9,height=12)
((LMA_1300_indval_plot/LDMC_1300_indval_plot/EWT_1300_indval_plot)|
    (N_1300_indval_plot/C_1300_indval_plot/plot_spacer()))+
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')
dev.off()