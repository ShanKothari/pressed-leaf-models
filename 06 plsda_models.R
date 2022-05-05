setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(caret)
library(plyr)

##################################
## to dos

## redo the divisions into training and testing data so that
## training data are the same for each analysis?

######################################################
## read data

fresh_spec_all<-readRDS("ProcessedSpectralData/fresh_spec_all.rds")
pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_all.rds")
ground_spec_all<-readRDS("ProcessedSpectralData/ground_spec_all.rds")

#####################################################
## shorten species names

meta(fresh_spec_all)$SpShort<-factor(unlist(lapply(strsplit(as.character(meta(fresh_spec_all)$Species),split=" "),
                                               function(x) paste(substring(x[[1]], 1, 1),x[[2]],sep = ". "))))
meta(pressed_spec_all)$SpShort<-factor(unlist(lapply(strsplit(as.character(meta(pressed_spec_all)$Species),split=" "),
                                               function(x) paste(substring(x[[1]], 1, 1),x[[2]],sep = ". "))))
meta(ground_spec_all)$SpShort<-factor(unlist(lapply(strsplit(as.character(meta(ground_spec_all)$Species),split=" "),
                                               function(x) paste(substring(x[[1]], 1, 1),x[[2]],sep = ". "))))

#####################################################
## fresh spectra

common_sp<-names(which(table(meta(fresh_spec_all)$SpShort)>=20))
fresh_spec_EL_common<-fresh_spec_all[which(meta(fresh_spec_all)$SpShort %in% common_sp)]
meta(fresh_spec_EL_common)$SpShort<-droplevels(meta(fresh_spec_EL_common)$SpShort)

fresh_spec_EL_class<-meta(fresh_spec_EL_common)$SpShort

## split data into 60% training and 40% testing datasets, stratified by species

train <- createDataPartition(
  y = meta(fresh_spec_EL_common)$SpShort,
  p = .6,
  list = FALSE
)

fresh_spec_EL_common_train<-as.matrix(fresh_spec_EL_common[train])
fresh_spec_EL_common_test<-as.matrix(fresh_spec_EL_common[-train])
fresh_spec_EL_class_train<-fresh_spec_EL_class[train]
fresh_spec_EL_class_test<-fresh_spec_EL_class[-train]

## here are two ways to deal with imbalance in #samples per group
## ctrl1 downsamples so that no group has more members than the smallest
## ctrl2 upsamples by duplicating members of smaller groups at random
## another option is to use sample="smote" but that works less well here

ctrl1 <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                     summaryFunction = multiClassSummary, sampling="down")

ctrl2 <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                     summaryFunction = multiClassSummary, sampling="up")

## since upsampling can lead to overfitting, we limit the number
## of components (specified by tuneLength) for the upsampling PLS-DA
## to the optimal # for the downsampling PLS-DA

plsFit_fresh1 <- train(fresh_spec_EL_common_train, fresh_spec_EL_class_train,
                      method = "pls", tuneLength = 50,
                      trControl = ctrl1,probMethod="softmax")

plsFit_fresh2 <- train(fresh_spec_EL_common_train, fresh_spec_EL_class_train,
                      method = "pls", tuneLength = plsFit_fresh1$bestTune$ncomp,
                      trControl = ctrl2,probMethod="softmax")

fresh_conmat<-confusionMatrix(predict(plsFit_fresh2, fresh_spec_EL_common_test),fresh_spec_EL_class_test)

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

fcm_d$RefSum<-unlist(lapply(fcm_d$Reference,function(x) sum(fcm_d$Freq[fcm_d$Reference==x])))
fcm_d$RefPer<-round(fcm_d$Freq/fcm_d$RefSum*100,digits=0)

fcm_st<-as.data.frame(fresh_conmat$overall)

fcm_d_p <-  ggplot(data = fcm_d[fcm_d$RefPer>0,],
                   aes(x = Prediction , y =  Reference, fill = RefPer))+
  geom_tile() +
  geom_text(aes(label = RefPer), color = 'peachpuff1', size = 6) +
  theme_light() +
  guides(fill=FALSE)+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7",
                        limits=c(0,100))+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="italic"),
        axis.text.y = element_text(face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Fresh-leaf spectra")+
  labs(tag="A",y="True identity")

pdf("Manuscript/Fig6A.pdf",width=7.5,height=7.5)
fcm_d_p
dev.off()

##################################################
## pressed spectra

## select only the most common species
common_sp<-names(which(table(meta(pressed_spec_all)$SpShort)>=20))
pressed_spec_EL_common<-pressed_spec_all[which(meta(pressed_spec_all)$SpShort %in% common_sp)]
meta(pressed_spec_EL_common)$SpShort<-droplevels(meta(pressed_spec_EL_common)$SpShort)

pressed_spec_EL_class<-meta(pressed_spec_EL_common)$SpShort

train <- createDataPartition(
  y = meta(pressed_spec_EL_common)$SpShort,
  p = .6,
  list = FALSE
)

pressed_spec_EL_common_train<-as.matrix(pressed_spec_EL_common[train])
pressed_spec_EL_common_test<-as.matrix(pressed_spec_EL_common[-train])
pressed_spec_EL_class_train<-pressed_spec_EL_class[train]
pressed_spec_EL_class_test<-pressed_spec_EL_class[-train]

ctrl1 <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                      summaryFunction = multiClassSummary, sampling="down")

ctrl2 <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                      summaryFunction = multiClassSummary, sampling="up")

plsFit_pressed1 <- train(pressed_spec_EL_common_train, pressed_spec_EL_class_train,
                       method = "pls", tuneLength = 50,
                       trControl = ctrl1,probMethod="softmax")

plsFit_pressed2 <- train(pressed_spec_EL_common_train, pressed_spec_EL_class_train,
                       method = "pls", tuneLength = plsFit_pressed1$bestTune$ncomp,
                       trControl = ctrl2,probMethod="softmax")

pressed_conmat<-confusionMatrix(predict(plsFit_pressed2, pressed_spec_EL_common_test),pressed_spec_EL_class_test)

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

pcm_d$RefSum<-unlist(lapply(pcm_d$Reference,function(x) sum(pcm_d$Freq[pcm_d$Reference==x])))
pcm_d$RefPer<-round(pcm_d$Freq/pcm_d$RefSum*100,digits=0)

pcm_st<-as.data.frame(pressed_conmat$overall)

pcm_d_p <-  ggplot(data = pcm_d[pcm_d$RefPer>0,],
                   aes(x = Prediction , y =  Reference, fill = RefPer))+
  geom_tile() +
  geom_text(aes(label = RefPer), color = 'peachpuff1', size = 6) +
  theme_light() +
  guides(fill=FALSE)+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7",
                        limits=c(0,100))+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="italic"),
        axis.text.y = element_text(face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Pressed-leaf spectra")+
  labs(tag="B",y="True identity")

pdf("Manuscript/Fig6B.pdf",width=7.5,height=7.5)
pcm_d_p
dev.off()

##################################################
## ground spectra

common_sp<-names(which(table(meta(ground_spec_all)$SpShort)>=20))
ground_spec_EL_common<-ground_spec_all[which(meta(ground_spec_all)$SpShort %in% common_sp)]
meta(ground_spec_EL_common)$SpShort<-droplevels(meta(ground_spec_EL_common)$SpShort)

ground_spec_EL_class<-meta(ground_spec_EL_common)$SpShort

train <- createDataPartition(
  y = meta(ground_spec_EL_common)$SpShort,
  p = .6,
  list = FALSE
)

ground_spec_EL_common_train<-as.matrix(ground_spec_EL_common[train])
ground_spec_EL_common_test<-as.matrix(ground_spec_EL_common[-train])
ground_spec_EL_class_train<-ground_spec_EL_class[train]
ground_spec_EL_class_test<-ground_spec_EL_class[-train]

ctrl1 <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                      summaryFunction = multiClassSummary, sampling="down")

ctrl2 <- trainControl(method = "repeatedcv", repeats = 10, number=10,
                      summaryFunction = multiClassSummary, sampling="up")

plsFit_ground1 <- train(ground_spec_EL_common_train, ground_spec_EL_class_train,
                       method = "pls", tuneLength = 50,
                       trControl = ctrl1,probMethod="softmax")

plsFit_ground2 <- train(ground_spec_EL_common_train, ground_spec_EL_class_train,
                       method = "pls", tuneLength = plsFit_ground1$bestTune$ncomp,
                       trControl = ctrl2,probMethod="softmax")

ground_conmat<-confusionMatrix(predict(plsFit_ground2, ground_spec_EL_common_test),ground_spec_EL_class_test)

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

gcm_d$RefSum<-unlist(lapply(gcm_d$Reference,function(x) sum(gcm_d$Freq[gcm_d$Reference==x])))
gcm_d$RefPer<-round(gcm_d$Freq/gcm_d$RefSum*100,digits=0)

gcm_st<-as.data.frame(ground_conmat$overall)

gcm_d_p <-  ggplot(data = gcm_d[gcm_d$RefPer>0,],
                   aes(x = Prediction , y =  Reference, fill = RefPer))+
  geom_tile() +
  geom_text(aes(label = RefPer), color = 'peachpuff1', size = 6) +
  theme_light() +
  guides(fill=FALSE)+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7",
                        limits=c(0,100))+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="italic"),
        axis.text.y = element_text(face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Ground-leaf spectra")+
  labs(tag="C",y="True identity")

pdf("Manuscript/Fig6C.pdf",width=7.5,height=7.5)
gcm_d_p
dev.off()

obj_list<-list(fresh=list(fresh_dat=fcm_d,
                          fresh_sum=fcm_st,
                          fresh_comp=plsFit_fresh2$bestTune$ncomp),
               pressed=list(pressed_dat=pcm_d,
                            pressed_sum=pcm_st,
                            pressed_comp=plsFit_pressed2$bestTune$ncomp),
               ground=list(ground_dat=gcm_d,
                           ground_sum=gcm_st,
                           ground_comp=plsFit_ground2$bestTune$ncomp))
saveRDS(obj_list,"SavedResults/class_summ.rds")
