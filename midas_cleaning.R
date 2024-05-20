#***************** Up in the Air: The Effects of Air Pollution on Conflict in West Africa**********************

#        MASTER THESIS at the Georg-August University of Goettingen 

#                       Ankhbayar Delgerchuluun

#                               RSCRIPT
#                          ------------------
#                 Data cleaning for MiDAS DoD dust data

####_______________________README__________________________________________________####

#This R script carries out the data cleaning process of MiDAS DoD data. 
#Since the MiDAS datasets is daily and saved in different files for each day, I have to merge daily data into 
#monthly level and find the monthly mean. I also cropped the data set to West African region in this script 
# empirical results and abstract of the thesis 



####_________Required packages_____________________________________________________#####

###Load a required package
library(raster)
library(sp)
library(ncdf4)
library(terra)
### This R script is dedicated to import the dust dataset and create monthly average raster
##Source: https://zenodo.org/record/4244106
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation/AIC_estimation")
West<-shapefile('West.shp')
west_ext<-extent(West)
west_ext1<-extent(West)

####_________Data cleaning_____________________________________________________#####
###2003
setwd("D:/dust_midas/2003")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2003-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2003')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20030101.nc', 
                  varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                     varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  print(all_dust[i])
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2003<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2003, 'datasets/Dust/midas_wa_2003.tif', overwrite=TRUE)

###2004
setwd("D:/dust_midas/2004")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2004-01-01"), by = "day", length.out = 366)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2004')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20040101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:366){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2004<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2004, 'datasets/Dust/midas_wa_2004.tif', overwrite=TRUE)


###2005
setwd("D:/dust_midas/2005")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2005-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2005')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20050101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2005<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2005, 'datasets/Dust/midas_wa_2005.tif', overwrite=TRUE)


###2006
setwd("D:/dust_midas/2006")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2006-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2006')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20060101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2006<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2006, 'datasets/Dust/midas_wa_2006.tif', overwrite=TRUE)
rm(midas_daily)
###2007
setwd("D:/dust_midas/2007")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2007-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2007')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20070101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2007<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2007, 'datasets/Dust/midas_wa_2007.tif', overwrite=TRUE)
rm(midas_daily)

###2008
setwd("D:/dust_midas/2008")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2008-01-01"), by = "day", length.out = 366)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2008')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20080101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:366){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2008<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2008, 'datasets/Dust/midas_wa_2008.tif', overwrite=TRUE)
rm(midas_daily)

###2009
setwd("D:/dust_midas/2009")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2009-01-01"), by = "day", length.out = 364)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2009')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20090101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:364){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2009<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2009, 'datasets/Dust/midas_wa_2009.tif', overwrite=TRUE)
rm(midas_daily)

###2010
setwd("D:/dust_midas/2010")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2010-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2010')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20100101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2010<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2010, 'datasets/Dust/midas_wa_2010.tif', overwrite=TRUE)
rm(midas_daily)


###2011
setwd("D:/dust_midas/2011")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2011-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2011')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20110101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2011<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2011, 'datasets/Dust/midas_wa_2011.tif', overwrite=TRUE)
rm(midas_daily)


###2012
setwd("D:/dust_midas/2012")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2012-01-01"), by = "day", length.out = 366)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2012')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20120101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:366){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2012<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2012, 'datasets/Dust/midas_wa_2012.tif', overwrite=TRUE)
rm(midas_daily)

###2013
setwd("D:/dust_midas/2013")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2013-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2013')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20130101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2013<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2013, 'datasets/Dust/midas_wa_2013.tif', overwrite=TRUE)
rm(midas_daily)


###2014
setwd("D:/dust_midas/2014")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2014-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2014')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20140101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2014<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2014, 'datasets/Dust/midas_wa_2014.tif', overwrite=TRUE)
rm(midas_daily)

###2015
setwd("D:/dust_midas/2015")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2015-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2015')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20150101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2015<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2015, 'datasets/Dust/midas_wa_2015.tif', overwrite=TRUE)
rm(midas_daily)


###2016
setwd("D:/dust_midas/2016")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2016-01-01"), by = "day", length.out = 366)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2016')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20160101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:366){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2016<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2016, 'datasets/Dust/midas_wa_2016.tif', overwrite=TRUE)

rm(midas_daily)

###2017
setwd("D:/dust_midas/2017")
## defining the date and other auxaliary variables 
dates_day <- seq(as.Date("2017-01-01"), by = "day", length.out = 365)
date_month<-format(dates_day,"%Y-%m")
##Setting up the paths
dust_path<-paste0('D:/dust_midas/2017')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
length(all_dust)
## first version
midas_daily<-brick('MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-20170101.nc', 
                   varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
midas_daily<-crop(midas_daily, west_ext, snap="near")
midas_daily<-rast(midas_daily)
##loop
for (i in 2:365){
  midas_daily_<-brick(all_dust[i], 
                      varname="Modis-total-dust-optical-depth-at-550nm",stopIfNotEqualSpaced = FALSE)
  ##Cropping west africa
  midas_wa<-crop(midas_daily_, west_ext, snap="near")
  midas_wa<-rast(midas_wa)
  add(midas_daily)<- midas_wa
}

##Calculating monthly average
midas_wa_2017<-tapp(midas_daily, date_month, mean, na.rm=TRUE)
setwd("C:/Users/USER/OneDrive/Documents/3. THESIS/estimation")
writeRaster(midas_wa_2017, 'datasets/Dust/midas_wa_2017.tif', overwrite=TRUE)
