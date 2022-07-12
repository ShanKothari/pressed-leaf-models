setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(caret)

## here we separate data between calibration and validation
## subsets so that the same samples are used for validation
## in analyses of fresh, pressed, and ground leaves
## so that comparisons of accuracy are more interpretable
## (no variation due to the caprices of random cal/val splits)

######################################################
## read data

fresh_spec_all<-readRDS("ProcessedSpectralData/fresh_spec_all.rds")
pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_all.rds")
ground_spec_all<-readRDS("ProcessedSpectralData/ground_spec_all.rds")

#####################################################
## dataset summary

fresh_meta<-meta(fresh_spec_all)
pressed_meta<-meta(pressed_spec_all)
ground_meta<-meta(ground_spec_all)

## archiving data!
# write.csv(as.data.frame(fresh_spec_all),"ProcessedSpectralData/fresh_spec_avg.csv",row.names = F)
# write.csv(as.data.frame(pressed_spec_all),"ProcessedSpectralData/pressed_spec_avg.csv",row.names = F)
# write.csv(as.data.frame(ground_spec_all),"ProcessedSpectralData/ground_spec_avg.csv",row.names = F)

## the few lines below are just to help count numbers of samples
## per project, growth form, or species
fp_merge<-merge(fresh_meta,pressed_meta,by="ID",all=T)
fp_merge$GrowthForm.x[is.na(fp_merge$GrowthForm.x)]<-fp_merge$GrowthForm.y[is.na(fp_merge$GrowthForm.x)]
fp_merge$Species.x[is.na(fp_merge$Species.x)]<-fp_merge$Species.y[is.na(fp_merge$Species.x)]
fp_merge$Project.x[is.na(fp_merge$Project.x)]<-fp_merge$Project.y[is.na(fp_merge$Project.x)]

table(fp_merge$Species.x)[order(names(table(fp_merge$Species.x)))]

######################################################
## partition data

## which data are in all three data sets?
fresh_spec_int<-fresh_spec_all[Reduce(intersect,list(names(fresh_spec_all),
                                                     names(pressed_spec_all),
                                                     names(ground_spec_all)))]

## choose validation samples from among those shared in common
test_sample_fresh <- createDataPartition(
  y = meta(fresh_spec_int)$GrowthForm,
  p = .25,
  list = FALSE
)

## do the splits
fresh_spec_all_train<-fresh_spec_all[-test_sample_fresh,]
fresh_spec_all_test<-fresh_spec_all[test_sample_fresh,]
test_names<-names(fresh_spec_all_test)

test_sample_pressed<-which(names(pressed_spec_all) %in% test_names)

pressed_spec_all_train<-pressed_spec_all[-test_sample_pressed,]
pressed_spec_all_test<-pressed_spec_all[test_sample_pressed,]

test_sample_ground<-which(names(ground_spec_all) %in% test_names)

ground_spec_all_train<-ground_spec_all[-test_sample_ground,]
ground_spec_all_test<-ground_spec_all[test_sample_ground,]

#######################################
## write data

saveRDS(fresh_spec_all_train,"ProcessedSpectralData/fresh_spec_all_train.rds")
saveRDS(fresh_spec_all_test,"ProcessedSpectralData/fresh_spec_all_test.rds")
saveRDS(pressed_spec_all_train,"ProcessedSpectralData/pressed_spec_all_train.rds")
saveRDS(pressed_spec_all_test,"ProcessedSpectralData/pressed_spec_all_test.rds")
saveRDS(ground_spec_all_train,"ProcessedSpectralData/ground_spec_all_train.rds")
saveRDS(ground_spec_all_test,"ProcessedSpectralData/ground_spec_all_test.rds")