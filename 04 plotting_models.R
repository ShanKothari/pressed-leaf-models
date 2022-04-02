setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(reshape2)

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

#############################################
## VIP plots

VIP_fresh<-readRDS("SavedResults/VIP_fresh.rds")
VIP_pressed<-readRDS("SavedResults/VIP_pressed.rds")
VIP_ground<-readRDS("SavedResults/VIP_ground.rds")

focal_palette=palette(brewer.pal(8,name="Set2")[c(3,4,5,6,8,1,2)])

VIP_fresh_long<-melt(VIP_fresh,id.vars = "wavelength")
VIP_pressed_long<-melt(VIP_pressed,id.vars = "wavelength")
VIP_ground_long<-melt(VIP_ground,id.vars = "wavelength")

## fresh

VIP_fiber_fresh_plot<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("sol","hemi","cell","lign"),],
                             aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Fresh")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_pigments_fresh_plot<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("chlA","chlB","car"),],
                                aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_other_fresh_plot<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("C","N","LMA","LDMC","EWT"),],
                             aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_ICP1_fresh_plot<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("Al","Ca","Cu","Fe","K"),],
                            aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_ICP2_fresh_plot<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("Mg","Mn","Na","P","Zn"),],
                            aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

pdf("Manuscript/FigS16.pdf",width=8,height=14)
VIP_fiber_fresh_plot+VIP_pigments_fresh_plot+
  VIP_other_fresh_plot+VIP_ICP1_fresh_plot+
  VIP_ICP2_fresh_plot+plot_layout(ncol = 1)
dev.off()

## pressed

VIP_fiber_pressed_plot<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("sol","hemi","cell","lign"),],
                             aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Pressed")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_pigments_pressed_plot<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("chlA","chlB","car"),],
                                aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_other_pressed_plot<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("C","N","LMA","LDMC","EWT"),],
                             aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_ICP1_pressed_plot<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("Al","Ca","Cu","Fe","K"),],
                            aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_ICP2_pressed_plot<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("Mg","Mn","Na","P","Zn"),],
                            aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

pdf("Manuscript/FigS17.pdf",width=8,height=14)
VIP_fiber_pressed_plot+VIP_pigments_pressed_plot+
  VIP_other_pressed_plot+VIP_ICP1_pressed_plot+
  VIP_ICP2_pressed_plot+plot_layout(ncol = 1)
dev.off()

## ground

VIP_fiber_ground_plot<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("sol","hemi","cell","lign"),],
                               aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Ground")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_pigments_ground_plot<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("chlA","chlB","car"),],
                                  aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_other_ground_plot<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("C","N","LMA","LDMC","EWT"),],
                               aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_ICP1_ground_plot<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("Al","Ca","Cu","Fe","K"),],
                              aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

VIP_ICP2_ground_plot<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("Mg","Mn","Na","P","Zn"),],
                              aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))

pdf("Manuscript/FigS18.pdf",width=8,height=14)
VIP_fiber_ground_plot+VIP_pigments_ground_plot+
  VIP_other_ground_plot+VIP_ICP1_ground_plot+
  VIP_ICP2_ground_plot+plot_layout(ncol = 1)
dev.off()

## all together

VIP_ms1_fresh<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("LMA","EWT","cell","chlA"),],
                      aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette[1:4])+
  ggtitle("Fresh")+
  coord_cartesian(ylim=c(0,3.5))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  geom_hline(yintercept=0.8,linetype="dashed",size=1.5)+
  guides(color=F)

VIP_ms2_fresh<-ggplot(VIP_fresh_long[VIP_fresh_long$variable %in% c("N","K","Mn"),],
                      aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette[5:7])+
  coord_cartesian(ylim=c(0,3))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  geom_hline(yintercept=0.8,linetype="dashed",size=1.5)+
  guides(color=F)

VIP_ms1_pressed<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("LMA","EWT","cell","chlA"),],
                        aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette[1:4])+
  ggtitle("Pressed")+
  coord_cartesian(ylim=c(0,3.5))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  geom_hline(yintercept=0.8,linetype="dashed",size=1.5)+
  guides(color=F)

VIP_ms2_pressed<-ggplot(VIP_pressed_long[VIP_pressed_long$variable %in% c("N","K","Mn"),],
                        aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette[5:7])+
  coord_cartesian(ylim=c(0,3))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  geom_hline(yintercept=0.8,linetype="dashed",size=1.5)+
  guides(color=F)

VIP_ms1_ground<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("LMA","EWT","cell","chlA"),],
                        aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette[1:4],
                     labels=c("LMA","EWT","Cell","Chl a"))+
  ggtitle("Ground")+
  coord_cartesian(ylim=c(0,3.5))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  geom_hline(yintercept=0.8,linetype="dashed",size=1.5)

VIP_ms2_ground<-ggplot(VIP_ground_long[VIP_ground_long$variable %in% c("N","K","Mn"),],
                        aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  scale_color_manual(values=focal_palette[5:7])+
  coord_cartesian(ylim=c(0,3))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  geom_hline(yintercept=0.8,linetype="dashed",size=1.5)

pdf("Manuscript/Fig4.pdf",width=15,height=7)
(VIP_ms1_fresh/VIP_ms2_fresh)|
  (VIP_ms1_pressed/VIP_ms2_pressed)|
  (VIP_ms1_ground/VIP_ms2_ground)
dev.off()

#############################################
## read saved data

fresh_jack_df_list<-readRDS("SavedResults/fresh_jack_df_list.rds")
fresh_jack_coef_list<-readRDS("SavedResults/fresh_jack_coefs_list.rds")

pressed_jack_df_list<-readRDS("SavedResults/pressed_jack_df_list.rds")
pressed_jack_coef_list<-readRDS("SavedResults/pressed_jack_coefs_list.rds")

pressed_1300_jack_df_list<-readRDS("SavedResults/pressed_1300_jack_df_list.rds")
pressed_1300_jack_coef_list<-readRDS("SavedResults/pressed_1300_jack_coefs_list.rds")

ground_jack_df_list<-readRDS("SavedResults/ground_jack_df_list.rds")
ground_jack_coef_list<-readRDS("SavedResults/ground_jack_coefs_list.rds")

#######################################
## axis limits for validation plots

all.solubles<-c(fresh_jack_df_list$sol$Measured,
                fresh_jack_df_list$sol$pred_mean,
                pressed_jack_df_list$sol$pred_mean,
                ground_jack_df_list$sol$pred_mean)
solubles_upper<-max(all.solubles,na.rm=T)+3
solubles_lower<-min(all.solubles,na.rm=T)-3

all.hemicellulose<-c(fresh_jack_df_list$hemi$Measured,
                fresh_jack_df_list$hemi$pred_mean,
                pressed_jack_df_list$hemi$pred_mean,
                ground_jack_df_list$hemi$pred_mean)
hemicellulose_upper<-max(all.hemicellulose,na.rm=T)+3
hemicellulose_lower<-min(all.hemicellulose,na.rm=T)-3

all.cellulose<-c(fresh_jack_df_list$cell$Measured,
                     fresh_jack_df_list$cell$pred_mean,
                     pressed_jack_df_list$cell$pred_mean,
                     ground_jack_df_list$cell$pred_mean)
cellulose_upper<-max(all.cellulose,na.rm=T)+2
cellulose_lower<-min(all.cellulose,na.rm=T)-2

all.lignin<-c(fresh_jack_df_list$lign$Measured,
                     fresh_jack_df_list$lign$pred_mean,
                     pressed_jack_df_list$lign$pred_mean,
                     ground_jack_df_list$lign$pred_mean)
lignin_upper<-max(all.lignin,na.rm=T)+2
lignin_lower<-min(all.lignin,na.rm=T)-2

all.perN<-c(fresh_jack_df_list$N$Measured,
            fresh_jack_df_list$N$pred_mean,
            pressed_jack_df_list$N$pred_mean,
            ground_jack_df_list$N$pred_mean)
perN_upper<-max(all.perN,na.rm=T)+0.2
perN_lower<-min(all.perN,na.rm=T)-0.2

all.perC<-c(fresh_jack_df_list$C$Measured,
                     fresh_jack_df_list$C$pred_mean,
                     pressed_jack_df_list$C$pred_mean,
                     ground_jack_df_list$C$pred_mean)
perC_upper<-max(all.perC,na.rm=T)+2
perC_lower<-min(all.perC,na.rm=T)-2

all.LMA<-c(fresh_jack_df_list$LMA$Measured,
            fresh_jack_df_list$LMA$pred_mean,
            pressed_jack_df_list$LMA$pred_mean,
            ground_jack_df_list$LMA$pred_mean)
LMA_upper<-max(all.LMA,na.rm=T)+0.02
LMA_lower<-min(all.LMA,na.rm=T)-0.02

all.LDMC<-c(fresh_jack_df_list$LDMC$Measured,
            fresh_jack_df_list$LDMC$pred_mean,
            pressed_jack_df_list$LDMC$pred_mean,
            ground_jack_df_list$LDMC$pred_mean)
LDMC_upper<-max(all.LDMC,na.rm=T)+30
LDMC_lower<-min(all.LDMC,na.rm=T)-30

all.EWT<-c(fresh_jack_df_list$EWT$Measured,
            fresh_jack_df_list$EWT$pred_mean,
            pressed_jack_df_list$EWT$pred_mean,
            ground_jack_df_list$EWT$pred_mean)
EWT_upper<-max(all.EWT,na.rm=T)+0.002
EWT_lower<-min(all.EWT,na.rm=T)-0.002

all.chlA<-c(fresh_jack_df_list$chlA$Measured,
            fresh_jack_df_list$chlA$pred_mean,
            pressed_jack_df_list$chlA$pred_mean,
            ground_jack_df_list$chlA$pred_mean)
chlA_upper<-max(all.chlA,na.rm=T)+1
chlA_lower<-min(all.chlA,na.rm=T)-1

all.chlB<-c(fresh_jack_df_list$chlB$Measured,
            fresh_jack_df_list$chlB$pred_mean,
            pressed_jack_df_list$chlB$pred_mean,
            ground_jack_df_list$chlB$pred_mean)
chlB_upper<-max(all.chlB,na.rm=T)+0.4
chlB_lower<-min(all.chlB,na.rm=T)-0.4

all.car<-c(fresh_jack_df_list$car$Measured,
            fresh_jack_df_list$car$pred_mean,
            pressed_jack_df_list$car$pred_mean,
            ground_jack_df_list$car$pred_mean)
car_upper<-max(all.car,na.rm=T)+0.2
car_lower<-min(all.car,na.rm=T)-0.2

all.Al<-c(fresh_jack_df_list$Al$Measured,
           fresh_jack_df_list$Al$pred_mean,
           pressed_jack_df_list$Al$pred_mean,
           ground_jack_df_list$Al$pred_mean)
Al_upper<-max(all.Al,na.rm=T)+0.02
Al_lower<-min(all.Al,na.rm=T)-0.02

all.Ca<-c(fresh_jack_df_list$Ca$Measured,
           fresh_jack_df_list$Ca$pred_mean,
           pressed_jack_df_list$Ca$pred_mean,
           ground_jack_df_list$Ca$pred_mean)
Ca_upper<-max(all.Ca,na.rm=T)+2
Ca_lower<-min(all.Ca,na.rm=T)-2

all.Cu<-c(fresh_jack_df_list$Cu$Measured,
           fresh_jack_df_list$Cu$pred_mean,
           pressed_jack_df_list$Cu$pred_mean,
           ground_jack_df_list$Cu$pred_mean)
Cu_upper<-max(all.Cu,na.rm=T)+0.003
Cu_lower<-min(all.Cu,na.rm=T)-0.003

all.Fe<-c(fresh_jack_df_list$Fe$Measured,
           fresh_jack_df_list$Fe$pred_mean,
           pressed_jack_df_list$Fe$pred_mean,
           ground_jack_df_list$Fe$pred_mean)
Fe_upper<-max(all.Fe,na.rm=T)+0.03
Fe_lower<-min(all.Fe,na.rm=T)-0.03

all.K<-c(fresh_jack_df_list$K$Measured,
           fresh_jack_df_list$K$pred_mean,
           pressed_jack_df_list$K$pred_mean,
           ground_jack_df_list$K$pred_mean)
K_upper<-max(all.K,na.rm=T)+3
K_lower<-min(all.K,na.rm=T)-3

all.Mg<-c(fresh_jack_df_list$Mg$Measured,
           fresh_jack_df_list$Mg$pred_mean,
           pressed_jack_df_list$Mg$pred_mean,
           ground_jack_df_list$Mg$pred_mean)
Mg_upper<-max(all.Mg,na.rm=T)+0.5
Mg_lower<-min(all.Mg,na.rm=T)-0.5

all.Mn<-c(fresh_jack_df_list$Mn$Measured,
           fresh_jack_df_list$Mn$pred_mean,
           pressed_jack_df_list$Mn$pred_mean,
           ground_jack_df_list$Mn$pred_mean)
Mn_upper<-max(all.Mn,na.rm=T)+0.1
Mn_lower<-min(all.Mn,na.rm=T)-0.1

all.Na<-c(fresh_jack_df_list$Na$Measured,
           fresh_jack_df_list$Na$pred_mean,
           pressed_jack_df_list$Na$pred_mean,
           ground_jack_df_list$Na$pred_mean)
Na_upper<-max(all.Na,na.rm=T)+0.5
Na_lower<-min(all.Na,na.rm=T)-0.5

all.P<-c(fresh_jack_df_list$P$Measured,
           fresh_jack_df_list$P$pred_mean,
           pressed_jack_df_list$P$pred_mean,
           ground_jack_df_list$P$pred_mean)
P_upper<-max(all.P,na.rm=T)+0.5
P_lower<-min(all.P,na.rm=T)-0.5

all.Zn<-c(fresh_jack_df_list$Zn$Measured,
           fresh_jack_df_list$Zn$pred_mean,
           pressed_jack_df_list$Zn$pred_mean,
           ground_jack_df_list$Zn$pred_mean)
Zn_upper<-max(all.Zn,na.rm=T)+0.1
Zn_lower<-min(all.Zn,na.rm=T)-0.1

#######################################
## plotting

solubles_fresh_val_plot<-ggplot(fresh_jack_df_list$sol,
                                aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  ggtitle("Fresh-leaf spectra")+guides(color=F)

hemicellulose_fresh_val_plot<-ggplot(fresh_jack_df_list$hemi,
                                     aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

cellulose_fresh_val_plot<-ggplot(fresh_jack_df_list$cell,
                                 aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)

lignin_fresh_val_plot<-ggplot(fresh_jack_df_list$lign,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)

perC_fresh_val_plot<-ggplot(fresh_jack_df_list$C,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured C (%)",x="Predicted C (%)")+
  ggtitle("Fresh-leaf spectra")+guides(color=F)

perN_fresh_val_plot<-ggplot(fresh_jack_df_list$N,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured N (%)",x="Predicted N (%)")+
  guides(color=F)

LMA_fresh_val_plot<-ggplot(fresh_jack_df_list$LMA,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Fresh-leaf spectra")+
  labs(y=expression("Measured LMA (kg m"^-2*")"),x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)

LDMC_fresh_val_plot<-ggplot(fresh_jack_df_list$LDMC,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)

EWT_fresh_val_plot<-ggplot(fresh_jack_df_list$EWT,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured EWT (mm)"),x=expression("Predicted EWT (mm)"))+
  guides(color=F)

chlA_fresh_val_plot<-ggplot(fresh_jack_df_list$chlA,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  ggtitle("Fresh-leaf spectra")+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)

chlB_fresh_val_plot<-ggplot(fresh_jack_df_list$chlB,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)

car_fresh_val_plot<-ggplot(fresh_jack_df_list$car,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),x=expression("Predicted carotenoids (mg g"^-1*")"))+
  guides(color=F)

Al_fresh_val_plot<-ggplot(fresh_jack_df_list$Al,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  ggtitle("Fresh-leaf spectra")+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)

Ca_fresh_val_plot<-ggplot(fresh_jack_df_list$Ca,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)

Cu_fresh_val_plot<-ggplot(fresh_jack_df_list$Cu,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),x=expression("Predicted Cu (mg g"^-1*")"))+
  guides(color=F)

Fe_fresh_val_plot<-ggplot(fresh_jack_df_list$Fe,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  ggtitle("Fresh-leaf spectra")+
  labs(y=expression("Measured Fe (mg g"^-1*")"),x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)

K_fresh_val_plot<-ggplot(fresh_jack_df_list$K,
                         aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured K (mg g"^-1*")"),x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mg_fresh_val_plot<-ggplot(fresh_jack_df_list$Mg,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)

Mn_fresh_val_plot<-ggplot(fresh_jack_df_list$Mn,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  ggtitle("Fresh-leaf spectra")+
  labs(y=expression("Measured Mn (mg g"^-1*")"),x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)

Na_fresh_val_plot<-ggplot(fresh_jack_df_list$Na,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured Na (mg g"^-1*")"),x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)

P_fresh_val_plot<-ggplot(fresh_jack_df_list$P,
                         aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured P (mg g"^-1*")"),x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)

Zn_fresh_val_plot<-ggplot(fresh_jack_df_list$Zn,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=20),
        plot.margin=unit(c(0,0.3,0,0),"in"),
        legend.position = c(0.85, 0.25))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)

solubles_pressed_val_plot<-ggplot(pressed_jack_df_list$sol,
                                  aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),
                  ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  ggtitle("Pressed-leaf spectra")+guides(color=F)

hemicellulose_pressed_val_plot<-ggplot(pressed_jack_df_list$hemi,
                                       aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

cellulose_pressed_val_plot<-ggplot(pressed_jack_df_list$cell,
                                   aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)

lignin_pressed_val_plot<-ggplot(pressed_jack_df_list$lign,
                                aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)

perC_pressed_val_plot<-ggplot(pressed_jack_df_list$C,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured C (%)",x="Predicted C (%)")+
  ggtitle("Pressed-leaf spectra")+guides(color=F)

perN_pressed_val_plot<-ggplot(pressed_jack_df_list$N,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured N (%)",x="Predicted N (%)")+
  guides(color=F)

LMA_pressed_val_plot<-ggplot(pressed_jack_df_list$LMA,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Pressed-leaf spectra")+
  labs(y=expression("Measured LMA (kg m"^-2*")"),x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)

LDMC_pressed_val_plot<-ggplot(pressed_jack_df_list$LDMC,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)

EWT_pressed_val_plot<-ggplot(pressed_jack_df_list$EWT,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured EWT (mm)"),x=expression("Predicted EWT (mm)"))+
  guides(color=F)

chlA_pressed_val_plot<-ggplot(pressed_jack_df_list$chlA,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Pressed-leaf spectra")+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)

chlB_pressed_val_plot<-ggplot(pressed_jack_df_list$chlB,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)

car_pressed_val_plot<-ggplot(pressed_jack_df_list$car,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),x=expression("Predicted carotenoids (mg g"^-1*")"))+
  guides(color=F)

Al_pressed_val_plot<-ggplot(pressed_jack_df_list$Al,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Pressed-leaf spectra")+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)

Ca_pressed_val_plot<-ggplot(pressed_jack_df_list$Ca,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)

Cu_pressed_val_plot<-ggplot(pressed_jack_df_list$Cu,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),x=expression("Predicted Cu (mg g"^-1*")"))+
  guides(color=F)

Fe_pressed_val_plot<-ggplot(pressed_jack_df_list$Fe,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Pressed-leaf spectra")+
  labs(y=expression("Measured Fe (mg g"^-1*")"),x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)

K_pressed_val_plot<-ggplot(pressed_jack_df_list$K,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured K (mg g"^-1*")"),x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mg_pressed_val_plot<-ggplot(pressed_jack_df_list$Mg,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)

Mn_pressed_val_plot<-ggplot(pressed_jack_df_list$Mn,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Pressed-leaf spectra")+
  labs(y=expression("Measured Mn (mg g"^-1*")"),x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)

Na_pressed_val_plot<-ggplot(pressed_jack_df_list$Na,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Na (mg g"^-1*")"),x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)

P_pressed_val_plot<-ggplot(pressed_jack_df_list$P,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured P (mg g"^-1*")"),x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)

Zn_pressed_val_plot<-ggplot(pressed_jack_df_list$Zn,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)

solubles_ground_val_plot<-ggplot(ground_jack_df_list$sol,
                                 aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  ggtitle("Ground-leaf spectra")+guides(color=F)

hemicellulose_ground_val_plot<-ggplot(ground_jack_df_list$hemi,
                                      aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

cellulose_ground_val_plot<-ggplot(ground_jack_df_list$cell,
                                  aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)

lignin_ground_val_plot<-ggplot(ground_jack_df_list$lign,
                               aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=guide_legend(title="Growth form"))

perC_ground_val_plot<-ggplot(ground_jack_df_list$C,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured C (%)",x="Predicted C (%)")+
  ggtitle("Ground-leaf spectra")+guides(color=F)

perN_ground_val_plot<-ggplot(ground_jack_df_list$N,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y="Measured N (%)",x="Predicted N (%)")+
  guides(color=F)

LMA_ground_val_plot<-ggplot(ground_jack_df_list$LMA,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Ground-leaf spectra")+
  labs(y=expression("Measured LMA (kg m"^-2*")"),x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)

LDMC_ground_val_plot<-ggplot(ground_jack_df_list$LDMC,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)

EWT_ground_val_plot<-ggplot(ground_jack_df_list$EWT,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured EWT (mm)"),x=expression("Predicted EWT (mm)"))+
  guides(color=guide_legend(title="Growth form"))

chlA_ground_val_plot<-ggplot(ground_jack_df_list$chlA,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Ground-leaf spectra")+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)

chlB_ground_val_plot<-ggplot(ground_jack_df_list$chlB,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)

car_ground_val_plot<-ggplot(ground_jack_df_list$car,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),x=expression("Predicted carotenoids (mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

Al_ground_val_plot<-ggplot(ground_jack_df_list$Al,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Ground-leaf spectra")+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)

Ca_ground_val_plot<-ggplot(ground_jack_df_list$Ca,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)

Cu_ground_val_plot<-ggplot(ground_jack_df_list$Cu,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),x=expression("Predicted Cu (mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

Fe_ground_val_plot<-ggplot(ground_jack_df_list$Fe,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Ground-leaf spectra")+
  labs(y=expression("Measured Fe (mg g"^-1*")"),x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)

K_ground_val_plot<-ggplot(ground_jack_df_list$K,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured K (mg g"^-1*")"),x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mg_ground_val_plot<-ggplot(ground_jack_df_list$Mg,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

Mn_ground_val_plot<-ggplot(ground_jack_df_list$Mn,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Ground-leaf spectra")+
  labs(y=expression("Measured Mn (mg g"^-1*")"),x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)

Na_ground_val_plot<-ggplot(ground_jack_df_list$Na,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Na (mg g"^-1*")"),x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)

P_ground_val_plot<-ggplot(ground_jack_df_list$P,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured P (mg g"^-1*")"),x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

Zn_ground_val_plot<-ggplot(ground_jack_df_list$Zn,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

##########################################################
## plotting

pdf("Manuscript/FigS6.pdf",width=15,height=15)
(LMA_fresh_val_plot + LMA_pressed_val_plot + LMA_ground_val_plot) /
  (LDMC_fresh_val_plot + LDMC_pressed_val_plot + LDMC_ground_val_plot) / 
  (EWT_fresh_val_plot + EWT_pressed_val_plot + EWT_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/FigS7.pdf",width=16,height=15)
(perC_fresh_val_plot + perC_pressed_val_plot + perC_ground_val_plot) / 
  (perN_fresh_val_plot + perN_pressed_val_plot + perN_ground_val_plot) /
  (P_fresh_val_plot + P_pressed_val_plot + P_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/FigS8.pdf",width = 16,height = 20)
(solubles_fresh_val_plot + solubles_pressed_val_plot + solubles_ground_val_plot) /
  (hemicellulose_fresh_val_plot + hemicellulose_pressed_val_plot + hemicellulose_ground_val_plot) /
  (cellulose_fresh_val_plot + cellulose_pressed_val_plot + cellulose_ground_val_plot) /
  (lignin_fresh_val_plot + lignin_pressed_val_plot + lignin_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/FigS9.pdf",width=16,height=15)
(chlA_fresh_val_plot + chlA_pressed_val_plot + chlA_ground_val_plot) /
  (chlB_fresh_val_plot + chlB_pressed_val_plot + chlB_ground_val_plot) / 
  (car_fresh_val_plot + car_pressed_val_plot + car_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/FigS10.pdf",width=16,height=15)
(Al_fresh_val_plot + Al_pressed_val_plot + Al_ground_val_plot) /
  (Ca_fresh_val_plot + Ca_pressed_val_plot + Ca_ground_val_plot) /
  (Cu_fresh_val_plot + Cu_pressed_val_plot + Cu_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/FigS11.pdf",width=16,height=15,onefile=F)
(Fe_fresh_val_plot + Fe_pressed_val_plot + Fe_ground_val_plot) /
  (K_fresh_val_plot + K_pressed_val_plot + K_ground_val_plot) /
  (Mg_fresh_val_plot + Mg_pressed_val_plot + Mg_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/FigS12.pdf",width=16,height=15,onefile=F)
(Mn_fresh_val_plot + Mn_pressed_val_plot + Mn_ground_val_plot) / 
  (Na_fresh_val_plot + Na_pressed_val_plot + Na_ground_val_plot) /
  (Zn_fresh_val_plot + Zn_pressed_val_plot + Zn_ground_val_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

## main text versions

chlA_fresh_val_plot_alt<-ggplot(fresh_jack_df_list$chlA,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)

chlA_pressed_val_plot_alt<-ggplot(pressed_jack_df_list$chlA,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)

EWT_ground_val_plot_alt<-ggplot(ground_jack_df_list$EWT,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured EWT (mm)"),x=expression("Predicted EWT (mm)"))+
  guides(color=F)

chlA_ground_val_plot_alt<-ggplot(ground_jack_df_list$chlA,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

pdf("Manuscript/Fig2.pdf",width=16,height=20)
(LMA_fresh_val_plot + LMA_pressed_val_plot + LMA_ground_val_plot) /
  (EWT_fresh_val_plot + EWT_pressed_val_plot + EWT_ground_val_plot_alt) /
  (cellulose_fresh_val_plot+cellulose_pressed_val_plot+cellulose_ground_val_plot) /
  (chlA_fresh_val_plot_alt+chlA_pressed_val_plot_alt+chlA_ground_val_plot_alt) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

# pdf("Manuscript/Fig2.pdf",width=16,height=20)
# (LMA_fresh_val_plot / EWT_fresh_val_plot / cellulose_fresh_val_plot / chlA_fresh_val_plot_alt)|
#   (LMA_pressed_val_plot / EWT_pressed_val_plot / cellulose_pressed_val_plot / chlA_pressed_val_plot_alt)|
#   (LMA_ground_val_plot / EWT_ground_val_plot_alt / cellulose_ground_val_plot / chlA_ground_val_plot_alt) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

perN_fresh_val_plot_alt<-ggplot(fresh_jack_df_list$N,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Fresh-leaf spectra")+
  labs(y="Measured N (%)",x="Predicted N (%)")+
  guides(color=F)

K_fresh_val_plot_alt<-ggplot(fresh_jack_df_list$K,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured K (mg g"^-1*")"),x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mn_fresh_val_plot_alt<-ggplot(fresh_jack_df_list$Mn,
                          aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)

perN_pressed_val_plot_alt<-ggplot(pressed_jack_df_list$N,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Pressed-leaf spectra")+
  labs(y="Measured N (%)",x="Predicted N (%)")+
  guides(color=F)

K_pressed_val_plot_alt<-ggplot(pressed_jack_df_list$K,
                               aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured K (mg g"^-1*")"),x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mn_pressed_val_plot_alt<-ggplot(pressed_jack_df_list$Mn,
                            aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)

perN_ground_val_plot_alt<-ggplot(ground_jack_df_list$N,
                             aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  ggtitle("Ground-leaf spectra")+
  labs(y="Measured N (%)",x="Predicted N (%)")+
  guides(color=F)

K_ground_val_plot_alt<-ggplot(ground_jack_df_list$K,
                              aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured K (mg g"^-1*")"),x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mn_ground_val_plot_alt<-ggplot(ground_jack_df_list$Mn,
                           aes(y=Measured,x=pred_mean,color=GrowthForm))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.85, 0.25),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=guide_legend(title="Growth form"))

pdf("Manuscript/Fig3.pdf",width=16,height=15)
(perN_fresh_val_plot_alt + perN_pressed_val_plot_alt + perN_ground_val_plot_alt) /
  (K_fresh_val_plot_alt + K_pressed_val_plot_alt + K_ground_val_plot_alt) /
  (Mn_fresh_val_plot_alt + Mn_pressed_val_plot_alt + Mn_ground_val_plot_alt) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

##################################################
## model performance

## calibration R2, RMSE, and %RMSE
## summary(lm(measured~val_pred,data=lignin_fresh_pred))
## with(lignin_fresh_pred,RMSD(measured,val_pred))
## with(lignin_fresh_pred,RMSD(measured,val_pred)/(max(measured$Measured,na.rm=T)-min(measured$Measured,na.rm=T)))

## validation R2, RMSE, and %RMSE
unlist(lapply(fresh_jack_df_list,function(x) x$ncomp[1]))
unlist(lapply(fresh_jack_df_list,function(x) summary(lm(Measured~pred_mean,data=x))$r.squared))
unlist(lapply(fresh_jack_df_list,function(x) RMSD(x$Measured,x$pred_mean)))
unlist(lapply(fresh_jack_df_list,function(x) percentRMSD(x$Measured,x$pred_mean,0.025,0.975)))*100
