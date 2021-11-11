setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(caret)

######################################################
## read data

fresh_spec_EL_agg<-readRDS("ProcessedSpectralData/fresh_spec_EL_agg.rds")
pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_EL_agg.rds")
ground_spec_EL_agg<-readRDS("ProcessedSpectralData/ground_spec_EL_agg.rds")

#####################################################
## dataset summary

fresh_meta<-meta(fresh_spec_EL_agg)
pressed_meta<-meta(pressed_spec_all)
ground_meta<-meta(ground_spec_EL_agg)

fp_merge<-merge(fresh_meta,pressed_meta,by="ID",all=T)
fp_merge$GrowthForm.x[is.na(fp_merge$GrowthForm.x)]<-fp_merge$GrowthForm.y[is.na(fp_merge$GrowthForm.x)]
fp_merge$Species.x[is.na(fp_merge$Species.x)]<-fp_merge$Species.y[is.na(fp_merge$Species.x)]
fp_merge$Project.x[is.na(fp_merge$Project.x)]<-fp_merge$Project.y[is.na(fp_merge$Project.x)]

table(fp_merge$Species.x)[order(names(table(fp_merge$Species.x)))]

######################################################
## partition data

## which data are in all three data sets?
fresh_spec_int<-fresh_spec_EL_agg[Reduce(intersect,list(names(fresh_spec_EL_agg),
                                                          names(pressed_spec_all),
                                                          names(ground_spec_EL_agg)))]

test_sample_fresh <- createDataPartition(
  y = meta(fresh_spec_int)$GrowthForm,
  p = .25,
  list = FALSE
)

fresh_spec_EL_agg_train<-fresh_spec_EL_agg[-test_sample_fresh,]
fresh_spec_EL_agg_test<-fresh_spec_EL_agg[test_sample_fresh,]
test_names<-names(fresh_spec_EL_agg_test)

test_sample_pressed<-which(names(pressed_spec_all) %in% test_names)

pressed_spec_all_train<-pressed_spec_all[-test_sample_pressed,]
pressed_spec_all_test<-pressed_spec_all[test_sample_pressed,]

test_sample_ground<-which(names(ground_spec_EL_agg) %in% test_names)

ground_spec_EL_agg_train<-ground_spec_EL_agg[-test_sample_ground,]
ground_spec_EL_agg_test<-ground_spec_EL_agg[test_sample_ground,]

#######################################
## write data

saveRDS(fresh_spec_EL_agg_train,"ProcessedSpectralData/fresh_spec_EL_agg_train.rds")
saveRDS(fresh_spec_EL_agg_test,"ProcessedSpectralData/fresh_spec_EL_agg_test.rds")
saveRDS(pressed_spec_all_train,"ProcessedSpectralData/pressed_spec_EL_agg_train.rds")
saveRDS(pressed_spec_all_test,"ProcessedSpectralData/pressed_spec_EL_agg_test.rds")
saveRDS(ground_spec_EL_agg_train,"ProcessedSpectralData/ground_spec_EL_agg_train.rds")
saveRDS(ground_spec_EL_agg_test,"ProcessedSpectralData/ground_spec_EL_agg_test.rds")
