# README: Pressed-leaf trait models

This repository is meant to support the manuscript [Reflectance spectroscopy allows rapid, accurate, and non-destructive estimates of functional traits from pressed leaves](https://www.biorxiv.org/content/10.1101/2021.04.21.440856v5) in press at _Methods in Ecology and Evolution_ and currently hosted at _bioRxiv_ (DOI: 10.1101/2021.04.21.440856). The repository is maintained by Shan Kothari (shan.kothari \[at\] umontreal \[dot\] ca). The current version of the repository is being archived at Zenodo (DOI forthcoming).

## Components

The repository contains R scripts numbered as belonging to 6 distinct 'stages' of analysis, as well as stage 00 (useful functions that are called in later scripts).

The stages are:

1. Reading reflectance spectra from fresh, pressed, and ground leaves; doing some processing on the pressed- and ground-leaf spectra; and attaching summary/trait data. Fresh-leaf spectra are already fairly processed when they are read in (see https://github.com/ShanKothari/CABO-trait-models)
2. Dividing up training and testing data for PLSR analyses.
3. Calibrating models using different kinds of data (\[A + E\] fresh-leaf spectra, \[B + D + F\] pressed-leaf spectra, \[C\] ground-leaf spectra) with different ranges (\[D\] restricted-range vs. \[all others\] full-range) and trait normalizations (\[E + F\] area-based vs \[all others\] mass-based chemical traits). Here I mostly follow the approach laid out by Burnett et al. (2021) _Journal of Experimental Botany_.
4. Plotting model predictions and VIP.
5. Externally validating the models on the Cedar Creek pressed-leaf dataset.
6. Running PLS-DA analyses to test the ability for fresh-, pressed-, and ground-leaf spectra to discriminate species.

The processed spectral and trait data products produced at the end of script 1, and read in at the beginning of script 2, are archived elsewhere (see **Associated data** below).

## How to use

This repository is *not* a software package or any sort of user-oriented product that people can use without further modification. It *is* meant to be a reasonably well-documented and faithful record of the analyses carried out in the two manuscripts listed above. Some analyses should be easily reproducible, with some modification, given the scripts and the archived data (see **Associated data** below). Users will have to (for example) change paths to the directories where they have saved the files locally.

## Associated data products

There are a few associated data products (DOIs forthcoming):

1. The main CABO dataset, including traits and spectra, are found [here](https://ecosis.org/package/cabo-2018-2019-leaf-level-spectra) at EcoSIS. Reflectance, transmittance, and absorptance spectra are available.
2. The complete Dessain project is not yet archived, but will be upon publication of manuscript (2).
3. [LOPEX](https://ecosis.org/package/leaf-optical-properties-experiment-database--lopex93-) and [ANGERS](https://ecosis.org/package/angers-leaf-optical-properties-database--2003-), used in the external validation, are available at those respective links
4. In script 10, I read in [fresh-](https://ecosis.org/package/fresh-leaf-cabo-spectra-from-herbarium-project), [pressed-](https://ecosis.org/package/pressed-leaf-cabo-spectra-from-herbarium-project), and [ground-](https://ecosis.org/package/pressed-leaf-cabo-spectra-from-herbarium-project)leaf spectra from another manuscript ([Kothari et al. 2022](https://www.biorxiv.org/content/10.1101/2021.04.21.440856), in press _Methods in Ecology and Evolution_) available at those respective links.

Even more raw data can be queried from the [CABO Data Portal](https://data.caboscience.org/leaf/).

In most of the scripts, I use the CRAN-hosted package `spectrolab v. 0.0.10` to handle the spectral data. In this package, the class `spectra` allows users to attach and retrieve metadata from spectral data using the function `meta()`. Below, you can find an example script that reads a .csv file, like our archived data, and turns it into an R `spectra` object like the ones I work with in the script.

```
library(spectrolab)

spec_df<-read.csv("mydata.csv")
name_var<-1 ## index for the column that contains sample names
meta_vars<-2:20 ## adjust as needed: indices for columns that contain metadata (including traits)
band_names<-400:2400 ## wavelengths of spectral bands corresponding to remaining columns

## you can also use the as_spectra command, but it's a bit more finicky 
## with data frames because the column names of bands must contain a letter
spec<-spectra(value = spec_df[,-c(name_vars,meta_vars)],
              band_names = 400:2400,
              names = spec_df[,name_var],
              meta = spec_df[,meta_vars])
```

## Questions?

Please feel free to reach me at _shan.kothari \[at\] umontreal \[dot\] ca_ with questions about the code or the data. I'm glad to help you adapt any of the approaches I use for your own purposes. If you draw heavily from my code, I would ask that (purely as a courtesy) you cite one of the two papers above⁠—whichever is more relevant.
