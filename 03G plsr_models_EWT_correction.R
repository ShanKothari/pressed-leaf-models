setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")

library(spectrolab)
library(ggplot2)
library(RColorBrewer)

##################################
## the purpose of this script is to develop new models
## for fresh EWT using 

fresh_spec_all<-readRDS("ProcessedSpectralData/fresh_spec_all.rds")
pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_all.rds")
ground_spec_all<-readRDS("ProcessedSpectralData/ground_spec_all.rds")

## 'fill in' EWT_actual (fresh leaves) for Dessain with
## EWT (rehydrated leaves), since no fresh mass was recorded
meta(fresh_spec_all)$EWT[meta(fresh_spec_all)$Project=="Dessain"]<-
  meta(fresh_spec_all)$EWT_rehydrated[meta(fresh_spec_all)$Project=="Dessain"]
meta(pressed_spec_all)$EWT[meta(pressed_spec_all)$Project=="Dessain"]<-
  meta(pressed_spec_all)$EWT_rehydrated[meta(pressed_spec_all)$Project=="Dessain"]
meta(ground_spec_all)$EWT[meta(ground_spec_all)$Project=="Dessain"]<-
  meta(ground_spec_all)$EWT_rehydrated[meta(ground_spec_all)$Project=="Dessain"]
