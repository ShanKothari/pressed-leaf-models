setwd("C:/Users/kotha020/Dropbox/TraitModels2018/HerbariumPaper/")
library(spectrolab)
library(ggplot2)
library(dplyr)
library(signal)
library(patchwork)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(geomtextpath)

########################################
## process ground spectra

ground_spec_EL<-read_spectra(path = "GroundSpectra/LaliberteLab","sed",exclude_if_matches = c("BAD","WR"))

# write.csv(as.matrix(ground_spec_EL),"ProcessedSpectralData/ground_spec_nonavg.csv")

## match sensors if desired
## not necessary, and not performed for the paper
# ground_spec_EL_match<-match_sensors(ground_spec_EL,splice_at=983)

ground_spec_EL_spl<-strsplit(x = names(ground_spec_EL), split = "_")
meta(ground_spec_EL)$ID<-unlist(lapply(ground_spec_EL_spl, function(x) x[[1]]))
meta(ground_spec_EL)$dummy<-0 # dummy variable just so metadata (ID) stays intact
## aggregate by sample
ground_spec_all_raw<-aggregate(ground_spec_EL,by=meta(ground_spec_EL)$ID,
                              FUN=try_keep_txt(mean))

## if smoothing is desired
## but this may obscure some small absorption features
## so I didn't do it for the paper
# ground_spec_all_vis<-t(apply(as.matrix(ground_spec_all_raw[,350:715]),1,sgolayfilt,p=3,n=21))
# ground_spec_all_nir<-t(apply(as.matrix(ground_spec_all_raw[,716:1390]),1,sgolayfilt,p=3,n=35))
# ground_spec_all_swir1<-t(apply(as.matrix(ground_spec_all_raw[,1391:1880]),1,sgolayfilt,p=3,n=35))
# ground_spec_all_swir2<-t(apply(as.matrix(ground_spec_all_raw[,1881:2500]),1,sgolayfilt,p=3,n=35))

# ground_spec_all<-cbind(ground_spec_all_vis,ground_spec_all_nir,ground_spec_all_swir1,ground_spec_all_swir2)
# colnames(ground_spec_all)<-350:2500
# ground_spec_all<-as_spectra(ground_spec_all)
# meta(ground_spec_all)$ID<-meta(ground_spec_all_raw)$ID
# ground_spec_all <- ground_spec_all[ ,400:2400]

# trim to 400-2400 nm
ground_spec_all <- ground_spec_all_raw[ ,400:2400]

#####################################################
## process pressed spectra

pressed_spec_EL<-read_spectra(path = "PressedSpectra/LaliberteLab","sed",exclude_if_matches = c("BAD","WR","BLACK"))

# write.csv(as.matrix(pressed_spec_EL),"ProcessedSpectralData/pressed_spec_nonavg.csv")

## for pressed spectra I did match sensors
pressed_spec_EL_match<-match_sensors(pressed_spec_EL,splice_at=983)

pressed_spec_EL_spl<-strsplit(x = names(pressed_spec_EL_match), split = "_")
meta(pressed_spec_EL_match)$ID<-as.factor(unlist(lapply(pressed_spec_EL_spl, function(x) x[[1]])))
meta(pressed_spec_EL_match)$dummy<-0
pressed_spec_all<-aggregate(pressed_spec_EL_match,by=meta(pressed_spec_EL_match)$ID,
                               FUN=try_keep_txt(mean))

## here again it is possible to smooth as one can do for ground spectra
## but again, not necessarily advisable

## trim to 400-2400 nm
pressed_spec_all <- pressed_spec_all[,400:2400]

#######################################################
## read and process fresh spectra

## most of the basic processing was done according to the script
## labeled 01 in the GitHub repository archived at
## DOI: 10.5281/zenodo.6820487
## here, we just read in the already-processed data
## (compiled across projects, and including
## some samples from the Dessain project found in the
## independent validation in that repository)

# source("FreshSpectra/process_fresh_spectra.R")
fresh_spec_df<-read.csv("FreshSpectra/fresh_spec_EL_processed.csv")
colnames(fresh_spec_df)<-gsub(pattern="X",replacement="",colnames(fresh_spec_df))
CABO_general_df<-read.csv("FreshSpectra/CABOGeneralOther_spec_processed.csv")
colnames(CABO_general_df)<-gsub(pattern="X",replacement="",colnames(CABO_general_df))

## trim ranges
fresh_spec<-as_spectra(fresh_spec_df,name_idx=1)
fresh_spec<-fresh_spec[ ,400:2400]
CABO_general<-as_spectra(CABO_general_df,name_idx=1)
CABO_general<-CABO_general[ ,400:2400]
fresh_spec_all<-spectrolab::combine(fresh_spec,CABO_general)

meta(fresh_spec_all)$ID<-c(as.character(fresh_spec_df$sample_id),CABO_general_df$sample_id)

##########################################################
## remove bad spectra

## removing yellow/red leaves, or other ones excluded
## because they deviated from the standard sampling protocol
bad_spectra<-c("11851575","13221003","22405244",
               "21854888","21854774","21854585",
               "10449089","10450505","10452131",
               "10453164","10454172","10461643",
               "10465061","10466192","10466735",
               "10468015")
fresh_spec_all<-fresh_spec_all[-which(meta(fresh_spec_all)$ID %in% bad_spectra)]
ground_spec_all<-ground_spec_all[-which(meta(ground_spec_all)$ID %in% bad_spectra)]
pressed_spec_all<-pressed_spec_all[-which(meta(pressed_spec_all)$ID %in% bad_spectra)]

## remove samples only found at one stage
unique_spectra<-union(union(setdiff(meta(fresh_spec_all)$ID,union(meta(ground_spec_all)$ID,meta(pressed_spec_all)$ID)),
                            setdiff(meta(ground_spec_all)$ID,union(meta(fresh_spec_all)$ID,meta(pressed_spec_all)$ID))),
                      setdiff(meta(pressed_spec_all)$ID,union(meta(ground_spec_all)$ID,meta(fresh_spec_all)$ID)))

fresh_spec_all<-fresh_spec_all[-which(meta(fresh_spec_all)$ID %in% unique_spectra)]
ground_spec_all<-ground_spec_all[-which(meta(ground_spec_all)$ID %in% unique_spectra)]
pressed_spec_all<-pressed_spec_all[-which(meta(pressed_spec_all)$ID %in% unique_spectra)]

#######################################################
## read in summary data

summary_BeauchampRioux<-read.csv("ELSummary/Summary_CABOData_ShanKProject - BeauchampRioux+Cabo-General.csv")
sBR_sub<-data.frame(ID=summary_BeauchampRioux$Bulk.sample.ID,
                    Species=summary_BeauchampRioux$Species,
                    Project="Beauchamp-Rioux",
                    Discoloration=summary_BeauchampRioux$Discoloration,
                    Stage=NA)

summary_Boucherville<-read.csv("ELSummary/Summary_CABOData_ShanKProject - Boucherville.csv")
sBV_sub<-data.frame(ID=summary_Boucherville$Bulk.sample.ID,
                    Species=summary_Boucherville$Species,
                    Project="Boucherville",
                    Discoloration=summary_Boucherville$Discoloration,
                    Stage=NA)

summary_Warren<-read.csv("ELSummary/Summary_CABOData_ShanKProject - SWA-Warren.csv")
sWN_sub<-data.frame(ID=summary_Warren$Bulk.sample.ID,
                    Species=summary_Warren$Species,
                    Project="Warren",
                    Discoloration=summary_Warren$Discoloration,
                    Stage=summary_Warren$Stage)

summary_Dessain<-read.csv("ELSummary/Summary_CABOData_ShanKProject - Dessain.csv")
sDN_sub<-data.frame(ID=summary_Dessain$Bulk.sample.ID,
                    Species=summary_Dessain$Species,
                    Project="Dessain",
                    Discoloration=summary_Dessain$Discoloration,
                    Stage=NA)

summary_all<-do.call(rbind,args=list(sBR_sub,sBV_sub,sWN_sub,sDN_sub))

## separate genus and species columns
sp_split<-strsplit(as.character(summary_all$Species),split=" ")
summary_all$LatinGenus<-unlist(lapply(sp_split,function(entry) entry[[1]]))
summary_all$LatinSpecies<-unlist(lapply(sp_split,function(entry) entry[[2]]))

## attach site-level latitude and longitude
## as well as the measurement date
Dessain_gps<-read.csv("ELSummary/Dessainleafspectra - gps.csv")
Dessain_date<-sapply(Dessain_gps$parentEventID,function(id) substr(id,start=1,stop=10))
Dessain_gps$date_measured<-as.POSIXlt(Dessain_date,format="%Y-%m-%d")
Dessain_gps_sub<-data.frame(ID=as.character(Dessain_gps$parentEventID),
                            Project="Dessain",
                            date_measured=Dessain_gps$date_measured,
                            latitude=Dessain_gps$decimalLatitude,
                            longitude=Dessain_gps$decimalLongitude)

spectra_summary<-read.csv("ELSummary/leaf_spectra.csv")
bulks<-read.csv("ELsummary/bulk_leaf_samples.csv")
plants<-read.csv("ELSummary/plants.csv")
sites<-read.csv("ELSummary/sites.csv")
spectra_summary$plant_id<-bulks$plant_id[match(spectra_summary$sample_id,bulks$sample_id)]

## grab lat/long from plants app if possible
spectra_summary$latitude<-plants$latitude[match(spectra_summary$plant_id,plants$plant_id)]
spectra_summary$longitude<-plants$longitude[match(spectra_summary$plant_id,plants$plant_id)]

## or otherwise sites app
sites_missing<-spectra_summary$site_id[is.na(spectra_summary$latitude)]
spectra_summary$latitude[is.na(spectra_summary$latitude)]<-sites$latitude[match(sites_missing,sites$site_id)]
spectra_summary$longitude[is.na(spectra_summary$longitude)]<-sites$longitude[match(sites_missing,sites$site_id)]

spectra_summary$date_measured<-as.POSIXlt(spectra_summary$date_measured,format="%Y-%m-%d")
spectra_summary_sub<-data.frame(ID=as.character(spectra_summary$sample_id),
                                Project=spectra_summary$project,
                                date_measured=spectra_summary$date_measured,
                                latitude=spectra_summary$latitude,
                                longitude=spectra_summary$longitude)

all_gps<-rbind(Dessain_gps_sub,spectra_summary_sub)

summary_all$latitude<-all_gps$latitude[match(summary_all$ID,all_gps$ID)]
summary_all$longitude<-all_gps$longitude[match(summary_all$ID,all_gps$ID)]
summary_all$date_measured<-all_gps$date_measured[match(summary_all$ID,all_gps$ID)]

## use VASCAN (Desmet & Brouillet 2013) for growth form and taxonomic group
## species included in multiple VASCAN growth forms are disambiguated manually
vascan<-read.csv("ELSummary/vascan.csv")
summary_all$GrowthForm<-vascan$growth_form[match(summary_all$Species,vascan$scientific_name)]
summary_all$GrowthForm[summary_all$Species=="Acer pensylvanicum Linnaeus"]<-"tree"
summary_all$GrowthForm[summary_all$Species=="Betula populifolia Marshall"]<-"tree"
summary_all$GrowthForm[summary_all$Species=="Staphylea trifolia Linnaeus"]<-"shrub"
summary_all$GrowthForm[summary_all$Species=="Prunus pensylvanica Linnaeus f."]<-"tree"
summary_all$GrowthForm[summary_all$Species=="Salix Linnaeus"]<-"shrub"
summary_all$GrowthForm[summary_all$Species=="Agonis flexuosa (Willd.) Sweet"]<-"tree"
summary_all$GrowthForm[summary_all$Species=="Amelanchier laevis Wiegand"]<-"tree"
summary_all$GrowthForm[summary_all$Species=="Calamagrostis canadensis Michaux Palisot de Beauvois"]<-"herb"

## all tree species here in the CABO dataset are broadleaf
levels(summary_all$GrowthForm)[levels(summary_all$GrowthForm)=="tree"]<-"broadleaf"

meta(fresh_spec_all)$Species<-
  summary_all$Species[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$LatinGenus<-
  summary_all$LatinGenus[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$LatinSpecies<-
  summary_all$LatinSpecies[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$Project<-
  summary_all$Project[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$Discoloration<-
  summary_all$Discoloration[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$Stage<-
  summary_all$Stage[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$GrowthForm<-
  summary_all$GrowthForm[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$Latitude<-
  summary_all$latitude[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$Longitude<-
  summary_all$longitude[match(meta(fresh_spec_all)$ID,summary_all$ID)]
meta(fresh_spec_all)$date_measured<-
  summary_all$date_measured[match(meta(fresh_spec_all)$ID,summary_all$ID)]
## this gets rid of any remaining 2017 Dessain spectra that aren't in the summary file
fresh_spec_all<-fresh_spec_all[-which(is.na(meta(fresh_spec_all)$Project))]

meta(pressed_spec_all)$Species<-
  summary_all$Species[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$LatinGenus<-
  summary_all$LatinGenus[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$LatinSpecies<-
  summary_all$LatinSpecies[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$Project<-
  summary_all$Project[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$Discoloration<-
  summary_all$Discoloration[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$Stage<-
  summary_all$Stage[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$GrowthForm<-
  summary_all$GrowthForm[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$Latitude<-
  summary_all$latitude[match(meta(pressed_spec_all)$ID,summary_all$ID)]
meta(pressed_spec_all)$Longitude<-
  summary_all$longitude[match(meta(pressed_spec_all)$ID,summary_all$ID)]
## should yield an error if all spectra IDs can be linked to a project
pressed_spec_all<-pressed_spec_all[-which(is.na(meta(pressed_spec_all)$Project))]

meta(ground_spec_all)$Species<-
  summary_all$Species[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$LatinGenus<-
  summary_all$LatinGenus[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$LatinSpecies<-
  summary_all$LatinSpecies[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$Project<-
  summary_all$Project[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$Discoloration<-
  summary_all$Discoloration[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$Stage<-
  summary_all$Stage[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$GrowthForm<-
  summary_all$GrowthForm[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$Latitude<-
  summary_all$latitude[match(meta(ground_spec_all)$ID,summary_all$ID)]
meta(ground_spec_all)$Longitude<-
  summary_all$longitude[match(meta(ground_spec_all)$ID,summary_all$ID)]
## should yield an error if all spectra IDs can be linked to a project
ground_spec_all<-ground_spec_all[-which(is.na(meta(ground_spec_all)$Project))]

## remove dummy variables
## for some reason, this doesn't actually work (must be done manually)
meta(pressed_spec_all)$dummy<-NULL
meta(ground_spec_all)$dummy<-NULL

############################################
## read in SLA / LDMC

SLA_Dessain<-read.csv("Traits/SLA/SLA_data_Aurelie_Dessain - Lab_data.csv")
SLA_Dessain<-data.frame(ID=SLA_Dessain$parentEventID,
                        SLA=SLA_Dessain$SLA_m2_kg,
                        LDMC=SLA_Dessain$LDMC_mg_g,
                        LDMC_actual=NA,
                        Thickness=NA)
SLA_Dessain$ID<-as.character(SLA_Dessain$ID)
SLA_Dessain$SLA<-as.numeric(as.character(SLA_Dessain$SLA))
SLA_Dessain$LDMC<-as.numeric(as.character(SLA_Dessain$LDMC))

SLA_other<-read.csv("Traits/SLA/leaf_area_and_water_samples.csv")
SLA_other<-data.frame(ID=SLA_other$sample_id,
                      SLA=as.numeric(as.character(SLA_other$specific_leaf_area_m2_kg)),
                      LDMC=as.numeric(as.character(SLA_other$leaf_dry_matter_content_mg_g)),
                      LDMC_actual=as.numeric(as.character(SLA_other$actual_leaf_dry_matter_content_perc))*10,
                      Thickness=NA)
SLA_other$ID<-as.character(SLA_other$ID)

SLA_all<-do.call(rbind,args=list(SLA_Dessain,SLA_other))
## flagged bad data points
SLA_all$SLA[SLA_all$ID==13404937]<-NA
SLA_all$LDMC[SLA_all$ID %in% c(10290262,10966273,13404937)]<-NA
SLA_all$LDMC_actual[SLA_all$ID %in% c(10290262,10966273,13404937)]<-NA

## LMA in kg/m2
SLA_all$LMA<-1/SLA_all$SLA
## EWT in mm
## based on either rehydrated (EWT) or field fresh (EWT_actual) leaves
SLA_all$EWT<-(1/(SLA_all$LDMC/1000)-1)*SLA_all$LMA
SLA_all$EWT_actual<-(1/(SLA_all$LDMC_actual/1000)-1)*SLA_all$LMA

meta(fresh_spec_all)$LDMC<-SLA_all$LDMC[match(meta(fresh_spec_all)$ID,SLA_all$ID)]
meta(fresh_spec_all)$LMA<-SLA_all$LMA[match(meta(fresh_spec_all)$ID,SLA_all$ID)]
meta(fresh_spec_all)$EWT<-SLA_all$EWT[match(meta(fresh_spec_all)$ID,SLA_all$ID)]
meta(fresh_spec_all)$EWT_actual<-SLA_all$EWT_actual[match(meta(fresh_spec_all)$ID,SLA_all$ID)]

meta(pressed_spec_all)$LDMC<-SLA_all$LDMC[match(meta(pressed_spec_all)$ID,SLA_all$ID)]
meta(pressed_spec_all)$LMA<-SLA_all$LMA[match(meta(pressed_spec_all)$ID,SLA_all$ID)]
meta(pressed_spec_all)$EWT<-SLA_all$EWT[match(meta(pressed_spec_all)$ID,SLA_all$ID)]
meta(pressed_spec_all)$EWT_actual<-SLA_all$EWT_actual[match(meta(pressed_spec_all)$ID,SLA_all$ID)]

meta(ground_spec_all)$LDMC<-SLA_all$LDMC[match(meta(ground_spec_all)$ID,SLA_all$ID)]
meta(ground_spec_all)$LMA<-SLA_all$LMA[match(meta(ground_spec_all)$ID,SLA_all$ID)]
meta(ground_spec_all)$EWT<-SLA_all$EWT[match(meta(ground_spec_all)$ID,SLA_all$ID)]
meta(ground_spec_all)$EWT_actual<-SLA_all$EWT_actual[match(meta(ground_spec_all)$ID,SLA_all$ID)]

## check relationship between EWT and EWT_actual
summary(lm(EWT_actual~EWT,data=meta(pressed_spec_all)))
with(meta(pressed_spec_all),quantile(EWT_actual/EWT,
                                     probs=c(0.025,0.25,0.5,0.75,0.975),
                                     na.rm=T))

## 'fill in' EWT_actual (fresh leaves) for Dessain with
## EWT (rehydrated leaves), since no fresh mass was recorded
meta(fresh_spec_all)$EWT_actual[meta(fresh_spec_all)$Project=="Dessain"]<-
  meta(fresh_spec_all)$EWT[meta(fresh_spec_all)$Project=="Dessain"]
meta(pressed_spec_all)$EWT_actual[meta(pressed_spec_all)$Project=="Dessain"]<-
  meta(pressed_spec_all)$EWT[meta(pressed_spec_all)$Project=="Dessain"]
meta(ground_spec_all)$EWT_actual[meta(ground_spec_all)$Project=="Dessain"]<-
  meta(ground_spec_all)$EWT[meta(ground_spec_all)$Project=="Dessain"]

############################################
## read in C/N

CN_BeauchampRioux<-read.csv("Traits/ElementalAnalysis/CN_data_all_projects - BeauchampRioux + CABO General.csv")

CN_Boucherville<-read.csv("Traits/ElementalAnalysis/CN_data_all_projects - Boucherville.csv")
CN_Boucherville$Sample_id<-as.factor(CN_Boucherville$Sample_id)
CN_Boucherville$Name<-as.factor(CN_Boucherville$Name)

CN_Warren<-read.csv("Traits/ElementalAnalysis/CN_data_all_projects - SWA-Warren.csv")
CN_Warren$Sample_id<-as.factor(CN_Warren$Sample_id)
CN_Warren$Name<-as.factor(CN_Warren$Name)
colnames(CN_Warren)[2]<-"No." ## reconciling column names

CN_Dessain<-read.csv("Traits/ElementalAnalysis/CN_data_all_projects - Dessain.csv")

CN_all<-do.call(rbind,args=list(CN_BeauchampRioux,CN_Boucherville,CN_Warren,CN_Dessain))

meta(fresh_spec_all)$N<-
  CN_all$N....[match(meta(fresh_spec_all)$ID,CN_all$Sample_id)]
meta(fresh_spec_all)$C<-
  CN_all$C.....[match(meta(fresh_spec_all)$ID,CN_all$Sample_id)]

meta(pressed_spec_all)$N<-
  CN_all$N....[match(meta(pressed_spec_all)$ID,CN_all$Sample_id)]
meta(pressed_spec_all)$C<-
  CN_all$C.....[match(meta(pressed_spec_all)$ID,CN_all$Sample_id)]

meta(ground_spec_all)$N<-
  CN_all$N....[match(meta(ground_spec_all)$ID,CN_all$Sample_id)]
meta(ground_spec_all)$C<-
  CN_all$C.....[match(meta(ground_spec_all)$ID,CN_all$Sample_id)]

###################################################
## read in carbon fractions

Fiber<-read.csv("Traits/CarbonFractions/carbon_fractions_bags.csv")
Fiber<-data.frame(ID=Fiber$bottle_id,
                  chem_samp=Fiber$leaf_chemistry_sample,
                  NDF=Fiber$ndf_perc,
                  ADF=Fiber$adf_perc,
                  ADL=Fiber$adl_perc,
                  solubles=Fiber$soluble_perc,
                  hemicellulose=Fiber$hemicellulose_perc,
                  cellulose=Fiber$cellulose_perc,
                  lignin=Fiber$lignin_perc)
Fiber$ID<-as.character(Fiber$ID)

meta(fresh_spec_all)$NDF<-
  Fiber$NDF[match(meta(fresh_spec_all)$ID,Fiber$ID)]
meta(fresh_spec_all)$ADF<-
  Fiber$ADF[match(meta(fresh_spec_all)$ID,Fiber$ID)]
meta(fresh_spec_all)$ADL<-
  Fiber$ADL[match(meta(fresh_spec_all)$ID,Fiber$ID)]
meta(fresh_spec_all)$solubles<-
  Fiber$solubles[match(meta(fresh_spec_all)$ID,Fiber$ID)]
meta(fresh_spec_all)$hemicellulose<-
  Fiber$hemicellulose[match(meta(fresh_spec_all)$ID,Fiber$ID)]
meta(fresh_spec_all)$cellulose<-
  Fiber$cellulose[match(meta(fresh_spec_all)$ID,Fiber$ID)]
meta(fresh_spec_all)$lignin<-
  Fiber$lignin[match(meta(fresh_spec_all)$ID,Fiber$ID)]

meta(pressed_spec_all)$NDF<-
  Fiber$NDF[match(meta(pressed_spec_all)$ID,Fiber$ID)]
meta(pressed_spec_all)$ADF<-
  Fiber$ADF[match(meta(pressed_spec_all)$ID,Fiber$ID)]
meta(pressed_spec_all)$ADL<-
  Fiber$ADL[match(meta(pressed_spec_all)$ID,Fiber$ID)]
meta(pressed_spec_all)$solubles<-
  Fiber$solubles[match(meta(pressed_spec_all)$ID,Fiber$ID)]
meta(pressed_spec_all)$hemicellulose<-
  Fiber$hemicellulose[match(meta(pressed_spec_all)$ID,Fiber$ID)]
meta(pressed_spec_all)$cellulose<-
  Fiber$cellulose[match(meta(pressed_spec_all)$ID,Fiber$ID)]
meta(pressed_spec_all)$lignin<-
  Fiber$lignin[match(meta(pressed_spec_all)$ID,Fiber$ID)]

meta(ground_spec_all)$NDF<-
  Fiber$NDF[match(meta(ground_spec_all)$ID,Fiber$ID)]
meta(ground_spec_all)$ADF<-
  Fiber$ADF[match(meta(ground_spec_all)$ID,Fiber$ID)]
meta(ground_spec_all)$ADL<-
  Fiber$ADL[match(meta(ground_spec_all)$ID,Fiber$ID)]
meta(ground_spec_all)$solubles<-
  Fiber$solubles[match(meta(ground_spec_all)$ID,Fiber$ID)]
meta(ground_spec_all)$hemicellulose<-
  Fiber$hemicellulose[match(meta(ground_spec_all)$ID,Fiber$ID)]
meta(ground_spec_all)$cellulose<-
  Fiber$cellulose[match(meta(ground_spec_all)$ID,Fiber$ID)]
meta(ground_spec_all)$lignin<-
  Fiber$lignin[match(meta(ground_spec_all)$ID,Fiber$ID)]

###############################################
## read pigments

pigments_Dessain<-read.csv("Traits/Pigments/Aurelie_pigments_valeurs_brutes - Analyses_Aurelie.csv")
pigments_Dessain<-pigments_Dessain[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

pigments_Boucherville<-read.csv("Traits/Pigments/Boucherville_pigments_valeurs_brutes - Sheet1.csv")
pigments_Boucherville<-pigments_Boucherville[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

pigments_Warren<-read.csv("Traits/Pigments/Warren_pigments_valeurs_brutes - Valeurs.csv")
pigments_Warren<-pigments_Warren[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

pigments_BeauchampRioux<-read.csv("Traits/Pigments/RBR_pigments_valeurs_brutes - valeurs.csv")
pigments_BeauchampRioux<-pigments_BeauchampRioux[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides_mg_g","Notes")]
colnames(pigments_BeauchampRioux)[colnames(pigments_BeauchampRioux)=="carotenoides_mg_g"]<-"carotenoides._mg_g"

pigments_all<-do.call(rbind,args=list(pigments_BeauchampRioux,pigments_Boucherville,pigments_Warren,pigments_Dessain))
pigments_all$chlA_mg_g<-as.numeric(as.character(pigments_all$chlA_mg_g))
pigments_all$chlB_mg_g<-as.numeric(as.character(pigments_all$chlB_mg_g))
pigments_all$carotenoides._mg_g<-as.numeric(as.character(pigments_all$carotenoides._mg_g))

## sample known to have bad pigment data
pigments_all$chlA_mg_g[pigments_all$sample_id==13404937]<-NA
pigments_all$chlB_mg_g[pigments_all$sample_id==13404937]<-NA
pigments_all$carotenoides._mg_g[pigments_all$sample_id==13404937]<-NA

## remove records who notes contain "refaire" or "refait" since these were all redone
pigments_all<-pigments_all[-which(str_detect(string = pigments_all$Notes,pattern = "efai")),]

meta(fresh_spec_all)$chlA<-
  pigments_all$chlA_mg_g[match(meta(fresh_spec_all)$ID,pigments_all$sample_id)]
meta(fresh_spec_all)$chlB<-
  pigments_all$chlB_mg_g[match(meta(fresh_spec_all)$ID,pigments_all$sample_id)]
meta(fresh_spec_all)$car<-
  pigments_all$carotenoides._mg_g[match(meta(fresh_spec_all)$ID,pigments_all$sample_id)]

meta(pressed_spec_all)$chlA<-
  pigments_all$chlA_mg_g[match(meta(pressed_spec_all)$ID,pigments_all$sample_id)]
meta(pressed_spec_all)$chlB<-
  pigments_all$chlB_mg_g[match(meta(pressed_spec_all)$ID,pigments_all$sample_id)]
meta(pressed_spec_all)$car<-
  pigments_all$carotenoides._mg_g[match(meta(pressed_spec_all)$ID,pigments_all$sample_id)]

meta(ground_spec_all)$chlA<-
  pigments_all$chlA_mg_g[match(meta(ground_spec_all)$ID,pigments_all$sample_id)]
meta(ground_spec_all)$chlB<-
  pigments_all$chlB_mg_g[match(meta(ground_spec_all)$ID,pigments_all$sample_id)]
meta(ground_spec_all)$car<-
  pigments_all$carotenoides._mg_g[match(meta(ground_spec_all)$ID,pigments_all$sample_id)]

######################################
## ICP elemental concentrations

ICP12<-read.csv("Traits/ICP/2020-03-12 1MHCl_Etienne box1,2.csv")
ICP34<-read.csv("Traits/ICP/2020-10-20 1MHCl_Etienne box 3,4.csv")
ICP5<-read.csv("Traits/ICP/2020-10-21 1MHCl_Etienne box 5.csv")
ICP67<-read.csv("Traits/ICP/2020-10-22 1MHCl_Etienne box 6,7.csv")
ICP8<-read.csv("Traits/ICP/2020-10-23 1MHCl_Etienne box 8.csv")
ICP9<-read.csv("Traits/ICP/2020-11-11 1MHCl_Etienne box 9.csv")
ICP1011<-read.csv("Traits/ICP/2020-12-16 1MHCl_Etienne box 10,11.csv")
ICP_boxes<-do.call(rbind,args=list(ICP12,ICP34,ICP5,ICP67,
                                 ICP8,ICP9,ICP1011))
num_cols<-c("Al","B","B.1","Ca","Cu","Fe","K","Mg","Mn","Na","P","Zn")
ICP_boxes[,num_cols]<-data.frame(sapply(ICP_boxes[,num_cols],
                                            function(x) as.numeric(as.character(x))))

ICP_Warren<-read.csv("Traits/ICP/icp_leaf_element_concentrations_Warren.csv")
ICP_Warren<-data.frame(Sample_id=Fiber$ID[match(ICP_Warren$leaf_chemistry_sample,Fiber$chem_samp)],
                       Scientific.name=NA,
                       Al=ICP_Warren$al_mg_g,
                       B=ICP_Warren$b_mg_g,
                       B.1=ICP_Warren$b_mg_g,
                       Ca=ICP_Warren$ca_mg_g,
                       Cu=ICP_Warren$cu_mg_g,
                       Fe=ICP_Warren$fe_mg_g,
                       K=ICP_Warren$k_mg_g,
                       Mg=ICP_Warren$mg_mg_g,
                       Mn=ICP_Warren$mn_mg_g,
                       Na=ICP_Warren$na_mg_g,
                       P=ICP_Warren$p_mg_g,
                       Zn=ICP_Warren$zn_mg_g)
ICP_all<-rbind(ICP_boxes,ICP_Warren)

## remove extreme outliers
ICP_all$Al[ICP_all$Sample_id=="2017-08-15-jbmcb-P006"]<-NA
ICP_all$Cu[ICP_all$Sample_id=="2017-06-07-ireqa-P010"]<-NA
ICP_all$Fe[ICP_all$Sample_id=="2017-05-31-jbmtb-P011"]<-NA
ICP_all$K[ICP_all$Sample_id=="2017-06-15-jbmcb-P003"]<-NA
ICP_all$Na[ICP_all$Sample_id=="2017-06-15-jbmcb-P003"]<-NA
## set samples with negative values (below detection limit) to 0
ICP_all$Al[ICP_all$Al<0]<-0
ICP_all$Na[ICP_all$Na<0]<-0

meta(fresh_spec_all)$Al<-
  ICP_all$Al[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Ca<-
  ICP_all$Ca[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Cu<-
  ICP_all$Cu[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Fe<-
  ICP_all$Fe[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$K<-
  ICP_all$K[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Mg<-
  ICP_all$Mg[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Mn<-
  ICP_all$Mn[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Na<-
  ICP_all$Na[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$P<-
  ICP_all$P[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]
meta(fresh_spec_all)$Zn<-
  ICP_all$Zn[match(meta(fresh_spec_all)$ID,ICP_all$Sample_id)]

meta(pressed_spec_all)$Al<-
  ICP_all$Al[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Ca<-
  ICP_all$Ca[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Cu<-
  ICP_all$Cu[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Fe<-
  ICP_all$Fe[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$K<-
  ICP_all$K[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Mg<-
  ICP_all$Mg[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Mn<-
  ICP_all$Mn[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Na<-
  ICP_all$Na[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$P<-
  ICP_all$P[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]
meta(pressed_spec_all)$Zn<-
  ICP_all$Zn[match(meta(pressed_spec_all)$ID,ICP_all$Sample_id)]

meta(ground_spec_all)$Al<-
  ICP_all$Al[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Ca<-
  ICP_all$Ca[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Cu<-
  ICP_all$Cu[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Fe<-
  ICP_all$Fe[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$K<-
  ICP_all$K[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Mg<-
  ICP_all$Mg[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Mn<-
  ICP_all$Mn[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Na<-
  ICP_all$Na[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$P<-
  ICP_all$P[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]
meta(ground_spec_all)$Zn<-
  ICP_all$Zn[match(meta(ground_spec_all)$ID,ICP_all$Sample_id)]

#####################################
## area based traits

## converting all to mg per cm2
meta(fresh_spec_all)$N_area<-meta(fresh_spec_all)$N*meta(fresh_spec_all)$LMA/1000
meta(fresh_spec_all)$C_area<-meta(fresh_spec_all)$C*meta(fresh_spec_all)$LMA/1000
meta(fresh_spec_all)$solubles_area<-meta(fresh_spec_all)$solubles*meta(fresh_spec_all)$LMA/1000
meta(fresh_spec_all)$hemicellulose_area<-meta(fresh_spec_all)$hemicellulose*meta(fresh_spec_all)$LMA/1000
meta(fresh_spec_all)$cellulose_area<-meta(fresh_spec_all)$cellulose*meta(fresh_spec_all)$LMA/1000
meta(fresh_spec_all)$lignin_area<-meta(fresh_spec_all)$lignin*meta(fresh_spec_all)$LMA/1000
meta(fresh_spec_all)$chlA_area<-meta(fresh_spec_all)$chlA*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$chlB_area<-meta(fresh_spec_all)$chlB*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$car_area<-meta(fresh_spec_all)$car*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Al_area<-meta(fresh_spec_all)$Al*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Ca_area<-meta(fresh_spec_all)$Ca*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Cu_area<-meta(fresh_spec_all)$Cu*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Fe_area<-meta(fresh_spec_all)$Fe*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$K_area<-meta(fresh_spec_all)$K*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Mg_area<-meta(fresh_spec_all)$Mg*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Mn_area<-meta(fresh_spec_all)$Mn*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Na_area<-meta(fresh_spec_all)$Na*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$P_area<-meta(fresh_spec_all)$P*meta(fresh_spec_all)$LMA/10000
meta(fresh_spec_all)$Zn_area<-meta(fresh_spec_all)$Zn*meta(fresh_spec_all)$LMA/10000

meta(pressed_spec_all)$N_area<-meta(pressed_spec_all)$N*meta(pressed_spec_all)$LMA/1000
meta(pressed_spec_all)$C_area<-meta(pressed_spec_all)$C*meta(pressed_spec_all)$LMA/1000
meta(pressed_spec_all)$solubles_area<-meta(pressed_spec_all)$solubles*meta(pressed_spec_all)$LMA/1000
meta(pressed_spec_all)$hemicellulose_area<-meta(pressed_spec_all)$hemicellulose*meta(pressed_spec_all)$LMA/1000
meta(pressed_spec_all)$cellulose_area<-meta(pressed_spec_all)$cellulose*meta(pressed_spec_all)$LMA/1000
meta(pressed_spec_all)$lignin_area<-meta(pressed_spec_all)$lignin*meta(pressed_spec_all)$LMA/1000
meta(pressed_spec_all)$chlA_area<-meta(pressed_spec_all)$chlA*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$chlB_area<-meta(pressed_spec_all)$chlB*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$car_area<-meta(pressed_spec_all)$car*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Al_area<-meta(pressed_spec_all)$Al*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Ca_area<-meta(pressed_spec_all)$Ca*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Cu_area<-meta(pressed_spec_all)$Cu*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Fe_area<-meta(pressed_spec_all)$Fe*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$K_area<-meta(pressed_spec_all)$K*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Mg_area<-meta(pressed_spec_all)$Mg*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Mn_area<-meta(pressed_spec_all)$Mn*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Na_area<-meta(pressed_spec_all)$Na*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$P_area<-meta(pressed_spec_all)$P*meta(pressed_spec_all)$LMA/10000
meta(pressed_spec_all)$Zn_area<-meta(pressed_spec_all)$Zn*meta(pressed_spec_all)$LMA/10000

######################################
## write files

saveRDS(fresh_spec_all,"ProcessedSpectralData/fresh_spec_all.rds")
saveRDS(pressed_spec_all,"ProcessedSpectralData/pressed_spec_all.rds")
saveRDS(ground_spec_all,"ProcessedSpectralData/ground_spec_all.rds")

######################################################
## read data again if needed

# fresh_spec_all<-readRDS("ProcessedSpectralData/fresh_spec_all.rds")
# pressed_spec_all<-readRDS("ProcessedSpectralData/pressed_spec_all.rds")
# ground_spec_all<-readRDS("ProcessedSpectralData/ground_spec_all.rds")

###############################################
## plot spectra quantiles

fresh_quantiles<-quantile(fresh_spec_all,probs=c(0.025,0.25,0.5,0.75,0.975))
pressed_quantiles<-quantile(pressed_spec_all,probs=c(0.025,0.25,0.5,0.75,0.975))
ground_quantiles<-quantile(ground_spec_all,probs=c(0.025,0.25,0.5,0.75,0.975))

fresh_CV<-apply(fresh_spec_all,2,function(x) sd(x)/mean(x))
pressed_CV<-apply(pressed_spec_all,2,function(x) sd(x)/mean(x))
ground_CV<-apply(ground_spec_all,2,function(x) sd(x)/mean(x))

fresh_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(fresh_quantiles)[1,],
                  ymax = as.matrix(fresh_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(fresh_quantiles)[2,],
                  ymax = as.matrix(fresh_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(fresh_quantiles)[3,]),size=1,color="black")+
  geom_line(aes(x=400:2400,y=fresh_CV),size=1,color="red")+
  geom_label(aes(x=400,y=0.93),hjust=0,label.size=NA,
             label="Fresh-leaf spectra",size=6,fill="white")+
  theme_bw()+
  theme(text = element_text(size=20),
#        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0,0.1,0.1,0),"in"))+
  labs(x="Wavelength (nm)",y="Reflectance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "A")

pressed_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(pressed_quantiles)[1,],
                  ymax = as.matrix(pressed_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(pressed_quantiles)[2,],
                  ymax = as.matrix(pressed_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(pressed_quantiles)[3,]),size=1,color="black")+
  geom_line(aes(x=400:2400,y=pressed_CV),size=1,color="red")+
  geom_label(aes(x=400,y=0.93),hjust=0,label.size=NA,
             label="Pressed-leaf spectra",size=6,fill="white")+
  theme_bw()+
  theme(text = element_text(size=20),
#        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0,0.1,0.1,0),"in"))+
  labs(x="Wavelength (nm)",y="Reflectance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "B")

ground_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(ground_quantiles)[1,],
                  ymax = as.matrix(ground_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(ground_quantiles)[2,],
                  ymax = as.matrix(ground_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(ground_quantiles)[3,]),size=1,color="black")+
  geom_line(aes(x=400:2400,y=ground_CV),size=1,color="red")+
  geom_label(aes(x=400,y=0.93),hjust=0,label.size=NA,
           label="Ground-leaf spectra",size=6,fill="white")+
  theme_bw()+
  theme(text = element_text(size=20),
#        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0,0.1,0.1,0),"in"))+
  labs(x="Wavelength (nm)",y="Reflectance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "C")

all_spec_plot<-ggplot()+
  geom_textline(aes(x=400:2400,y=as.matrix(fresh_quantiles)[3,]),
                linewidth=1,linetype="solid",label="Fresh",hjust=0.40,size=5)+
  geom_textline(aes(x=400:2400,y=as.matrix(pressed_quantiles)[3,]),
            linewidth=1,linetype="longdash",label="Pressed",hjust=0.45,size=5)+
  geom_textline(aes(x=400:2400,y=as.matrix(ground_quantiles)[3,]),
            linewidth=1,linetype="twodash",label="Ground",hjust=0.45,size=5)+
  theme_bw()+
  theme(text = element_text(size=20),
#        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.1,0,0),"in"))+
  labs(x="Wavelength (nm)",y="Reflectance")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "D")

pdf("Manuscript/Fig1.pdf",width=6,height=13)
fresh_spec_plot+pressed_spec_plot+ground_spec_plot+all_spec_plot+
  plot_layout(ncol = 1)
dev.off()

##################################
## other assorted plots

## spectra by discoloration stage
pressed_spec_disc<-aggregate(as.matrix(pressed_spec_all),
                             by=list(as.factor(meta(pressed_spec_all)$Discoloration)),
                             FUN=mean)
pressed_spec_disc_long<-melt(data = pressed_spec_disc,id.vars = "Group.1")
pressed_spec_disc_long$variable<-as.numeric(as.character(pressed_spec_disc_long$variable))

focal_palette=palette(brewer.pal(8,name="Set2")[c(3,4,5,6,8)])

pressed_disc_plot<-ggplot(pressed_spec_disc_long,
                          aes(x=variable,y=value,color=Group.1))+
  geom_line(size=1)+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0,0.1),"in"),
        legend.position = c(0.8,0.75))+
  labs(x="Wavelength (nm)",y="Reflectance",tag="A")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  guides(color=guide_legend(title="Discoloration"))+
  scale_color_manual(values=focal_palette)

## histogram of discoloration
disc_hist<-ggplot(meta(pressed_spec_all),
                  aes(x=as.factor(Discoloration)))+
  geom_bar()+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0,0.1),"in"))+
  labs(x="Discoloration",y="Count",tag="B")

pdf("Manuscript/disc_spec.pdf",width=7,height=10)
pressed_disc_plot+disc_hist+plot_layout(ncol=1)
dev.off()