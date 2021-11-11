setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(caret)
library(plyr)

######################################################
## read data

fresh_spec_EL_agg<-readRDS("ProcessedSpectralData/fresh_spec_EL_agg.rds")
pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_EL_agg.rds")
ground_spec_EL_agg<-readRDS("ProcessedSpectralData/ground_spec_EL_agg.rds")

#####################################################
## fresh spectra

common_sp<-names(which(table(meta(fresh_spec_EL_agg)$Species)>=20))
fresh_spec_EL_common<-fresh_spec_EL_agg[which(meta(fresh_spec_EL_agg)$Species %in% common_sp)]
meta(fresh_spec_EL_common)$Species<-droplevels(meta(fresh_spec_EL_common)$Species)
levels(meta(fresh_spec_EL_common)$Species)<-gsub("[^[:alnum:] ]",
                                                 replacement = "",
                                                 x = levels(meta(fresh_spec_EL_common)$Species))

fresh_spec_EL_common_df<-as.data.frame(fresh_spec_EL_common,metadata=F)
fresh_spec_EL_common_df$Species<-meta(fresh_spec_EL_common)$Species

fresh_spec_EL_common_train<-ddply(fresh_spec_EL_common_df,.(Species),function(x) x[sample(nrow(x),13),])
fresh_spec_EL_common_test<-fresh_spec_EL_common_df[!(fresh_spec_EL_common_df$sample_name %in% fresh_spec_EL_common_train$sample_name),]

fresh_spec_EL_common_train$sample_name<-NULL
fresh_spec_EL_common_test$sample_name<-NULL

ctrl <- trainControl(method = "repeatedcv", repeats = 10, number=5,
                     summaryFunction = multiClassSummary)

plsFit_fresh <- train(subset(fresh_spec_EL_common_train,select=-Species), fresh_spec_EL_common_train$Species,
                      method = "pls", tuneLength = 30,
                      trControl = ctrl,probMethod="softmax")

fresh_conmat<-confusionMatrix(predict(plsFit_fresh,
                                      subset(fresh_spec_EL_common_test,select=-Species)),
                              fresh_spec_EL_common_test$Species)

fcm_d<-as.data.frame(fresh_conmat$table)
fcm_d$Prediction<-unlist(lapply(as.character(fcm_d$Prediction),function(x){
  pred_split<-strsplit(x,split=" ")
  pred_paste<-paste(pred_split[[1]][1:2],collapse=" ")
  return(pred_paste)
}))

fcm_d$Reference<-unlist(lapply(as.character(fcm_d$Reference),function(x){
  pred_split<-strsplit(x,split=" ")
  pred_paste<-paste(pred_split[[1]][1:2],collapse=" ")
  return(pred_paste)
}))
fcm_st<-as.data.frame(fresh_conmat$overall)

fcm_d_p <-  ggplot(data = fcm_d, aes(x = Prediction , y =  Reference, fill = Freq))+
  geom_tile() +
  geom_text(aes(label = paste("",Freq)), color = 'dark gray', size = 8) +
  theme_light() +
  guides(fill=FALSE)+theme(text=element_text(size=15),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Fresh leaf PLS-DA")

##################################################
## pressed spectra

## select only the most common species
common_sp<-names(which(table(meta(pressed_spec_all)$Species)>=20))
pressed_spec_EL_common<-pressed_spec_all[which(meta(pressed_spec_all)$Species %in% common_sp)]
meta(pressed_spec_EL_common)$Species<-droplevels(meta(pressed_spec_EL_common)$Species)
levels(meta(pressed_spec_EL_common)$Species)<-gsub("[^[:alnum:] ]",
                                                   replacement = "",
                                                   x = levels(meta(pressed_spec_EL_common)$Species))

pressed_spec_EL_common_df<-as.data.frame(pressed_spec_EL_common,metadata=F)
pressed_spec_EL_common_df$Species<-meta(pressed_spec_EL_common)$Species

pressed_spec_EL_common_train<-ddply(pressed_spec_EL_common_df,.(Species),function(x) x[sample(nrow(x),13),])
pressed_spec_EL_common_test<-pressed_spec_EL_common_df[!(pressed_spec_EL_common_df$sample_name %in% pressed_spec_EL_common_train$sample_name),]

pressed_spec_EL_common_train$sample_name<-NULL
pressed_spec_EL_common_test$sample_name<-NULL

ctrl <- trainControl(method = "repeatedcv", repeats = 10, number=5,
                     summaryFunction = multiClassSummary)

plsFit_pressed <- train(subset(pressed_spec_EL_common_train,select=-Species), pressed_spec_EL_common_train$Species,
                        method = "pls", tuneLength = 30,
                        trControl = ctrl,probMethod="softmax")

pressed_conmat<-confusionMatrix(predict(plsFit_pressed,
                                        subset(pressed_spec_EL_common_test,select=-Species)),
                                pressed_spec_EL_common_test$Species)

pcm_d<-as.data.frame(pressed_conmat$table)
pcm_d$Prediction<-unlist(lapply(as.character(pcm_d$Prediction),function(x){
  pred_split<-strsplit(x,split=" ")
  pred_paste<-paste(pred_split[[1]][1:2],collapse=" ")
  return(pred_paste)
}))

pcm_d$Reference<-unlist(lapply(as.character(pcm_d$Reference),function(x){
  pred_split<-strsplit(x,split=" ")
  pred_paste<-paste(pred_split[[1]][1:2],collapse=" ")
  return(pred_paste)
}))
pcm_st<-as.data.frame(pressed_conmat$overall)

pcm_d_p <-  ggplot(data = pcm_d, aes(x = Prediction , y =  Reference, fill = Freq))+
  geom_tile() +
  geom_text(aes(label = paste("",Freq)), color = 'dark gray', size = 8) +
  theme_light() +
  guides(fill=FALSE)+theme(text=element_text(size=15),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Pressed leaf PLS-DA")

##################################################
## ground spectra

common_sp<-names(which(table(meta(ground_spec_EL_agg)$Species)>=20))
ground_spec_EL_common<-ground_spec_EL_agg[which(meta(ground_spec_EL_agg)$Species %in% common_sp)]
meta(ground_spec_EL_common)$Species<-droplevels(meta(ground_spec_EL_common)$Species)
levels(meta(ground_spec_EL_common)$Species)<-gsub("[^[:alnum:] ]",
                                                  replacement = "",
                                                  x = levels(meta(ground_spec_EL_common)$Species))

ground_spec_EL_common_df<-as.data.frame(ground_spec_EL_common,metadata=F)
ground_spec_EL_common_df$Species<-meta(ground_spec_EL_common)$Species

ground_spec_EL_common_train<-ddply(ground_spec_EL_common_df,.(Species),function(x) x[sample(nrow(x),13),])
ground_spec_EL_common_test<-ground_spec_EL_common_df[!(ground_spec_EL_common_df$sample_name %in% ground_spec_EL_common_train$sample_name),]

ground_spec_EL_common_train$sample_name<-NULL
ground_spec_EL_common_test$sample_name<-NULL

ctrl <- trainControl(method = "repeatedcv", repeats = 10, number=5,
                     summaryFunction = multiClassSummary)

plsFit_ground <- train(subset(ground_spec_EL_common_train,select=-Species), ground_spec_EL_common_train$Species,
                        method = "pls", tuneLength = 30,
                        trControl = ctrl,probMethod="softmax")

ground_conmat<-confusionMatrix(predict(plsFit_ground,
                                        subset(ground_spec_EL_common_test,select=-Species)),
                                ground_spec_EL_common_test$Species)

gcm_d<-as.data.frame(ground_conmat$table)
gcm_d$Prediction<-unlist(lapply(as.character(gcm_d$Prediction),function(x){
  pred_split<-strsplit(x,split=" ")
  pred_paste<-paste(pred_split[[1]][1:2],collapse=" ")
  return(pred_paste)
}))

gcm_d$Reference<-unlist(lapply(as.character(gcm_d$Reference),function(x){
  pred_split<-strsplit(x,split=" ")
  pred_paste<-paste(pred_split[[1]][1:2],collapse=" ")
  return(pred_paste)
}))
gcm_st<-as.data.frame(ground_conmat$overall)

gcm_d_p <-  ggplot(data = gcm_d, aes(x = Prediction , y =  Reference, fill = Freq))+
  geom_tile() +
  geom_text(aes(label = paste("",Freq)), color = 'dark gray', size = 8) +
  theme_light() +
  guides(fill=FALSE)+theme(text=element_text(size=15),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Ground leaf PLS-DA")

pdf("Manuscript/Fig5.pdf",width=7.5,height=7.5)
fcm_d_p
pcm_d_p
gcm_d_p
dev.off()
