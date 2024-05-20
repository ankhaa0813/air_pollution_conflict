
#***************** Up in the Air: The Effects of Air Pollution on Conflict in West Africa**********************

#        MASTER THESIS at the Georg-August University of Goettingen 

#                       Ankhbayar Delgerchuluun

#                               R SCRIPT
#                          ------------------
#                            Data Cleaning

####_______________________README__________________________________________________####

#This R script carries out all pre-analysis process, such as data importing, cleaning, and merging. 
#Thus, it consists of 3 main sections and 35 sub-sections. In the beginning of the main sections, 
#I provided a very brief objective of the section.

####_________Required packages_____________________________________________________#####
#packages                                  

library(dplyr, warn.conflicts = FALSE)
library(sf)
library(spData)
library("sp")
library("ggplot2")
library(ncdf4)
library(raster)
library(terra)
library(stargazer)
library(lubridate)
library(geodata)
library("murdock")
library("osmdata")
library("OpenStreetMap")
library("psych")              
library("rbin")               

#Note for packages 
#Please be aware the conflicts between raster and terra packages. Terra package doesn't 
#recognize a object created by the raster package. It is also vise versa

####_______1.DATA IMPORTING _______________________________________________________####

# In this section of the code, I import the raster datas from different sources, crop them to West africa,
#crop time dimension between 2003 and 2009, and project them with with EGS4386 format. 
####_____________1.1 Defining the time dimension___________________________________####
start_date<-as.Date("2003-01-01")
end_date<-as.Date("2009-4-30")
month1 <- function(x) as.Date(cut(x, "month"))
months <- seq(month1(start_date), month1(end_date), "month")
months<-format(months,"%Y-%m")
head(months)
#Importing a shape file of West African map
West<-gadm(c("NGA","CPV", "GHA","BEN","TGO","GMB","BFA","GNB","LBR","MLI","MRT","NER","SEN","SLE","CIV","GIN"),0,'datasets/population')

a<-vect("datasets/ethnic/Ethnic_Groups_in_Africa.shp")
plot(a)
print(a)

####_____________1.2 Conflict data set from ACLED__________________________________####

#Calling raw data
data_raw_conflict<-read.csv("datasets/Conflict/ACLED/2003-01-01-2009-04-30-Western_Africa.csv")
unique_id<-unique(data_raw_conflict$event_id_cnty)

head(data_raw_conflict)
data_raw_conflict<-filter(data_raw_conflict, region=="Western Africa" | country=="Cabo Verde")
table(data_raw_conflict$country)
#There is no conflict for Cabo Verde
data_raw_conflict<-data.frame(data_raw_conflict$latitude, data_raw_conflict$longitude, 
                              data_raw_conflict$event_date, data_raw_conflict$fatalities, 
                              data_raw_conflict$year, data_raw_conflict$event_type, data_raw_conflict$country)
data_raw_conflict<-mutate(data_raw_conflict, count=1)
head(data_raw_conflict)
names(data_raw_conflict)[1]="latitude" #y axis
names(data_raw_conflict)[2]="longitude" #x axis
names(data_raw_conflict)[3]="event_date"
names(data_raw_conflict)[4]="fatalities"
names(data_raw_conflict)[5]="year"
names(data_raw_conflict)[6]="type"
names(data_raw_conflict)[7]="country"
write.csv(data_raw_conflict, "datasets/final_temp/conflict_check.csv", row.names=FALSE)
head(data_raw_conflict)
summary(data_raw_conflict)

#adding year-month variable 
data_raw_conflict<-mutate(data_raw_conflict, date_conflict=format(as.Date(data_raw_conflict$event_date), "%Y-%m"))
head(data_raw_conflict$date_conflict)
data_raw_conflict<-mutate(data_raw_conflict, date=as.Date(event_date))
data_raw_conflict<-filter(data_raw_conflict, date_conflict<="2009-04")
table(data_raw_conflict$date_conflict)
head(data_raw_conflict)
str(data_raw_conflict)
data_raw_conflict<-data_raw_conflict[order(data_raw_conflict$date),]
print(unique(data_raw_conflict$date_conflict))
write.csv(data_raw_conflict, "datasets/final_temp/conflict_check.csv", row.names=FALSE)
conflict_g<-data_raw_conflict %>%
  group_by(year, type)%>% summarise(conflict=sum(count), death=sum(fatalities))
write.csv(conflict_g, "datasets/Conflict/ACLED/conflict_g.csv", row.names=FALSE)

#creating a initial data points
data_2003 <- filter(data_raw_conflict, date_conflict=="2003-01")
table(data_2003$date)
sum(data_2003[, "fatalities"])
conflict_map<-cbind(data_2003$longitude, data_2003$latitude)
plot(conflict_map)
summary(conflict_map)
dim(conflict_map)
death<-data_2003$fatalities
table(death)


conflict<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")

conflict_death<-rasterize(conflict_map, conflict, death, fun="sum", na.rm=TRUE) 
conflict_count<-rasterize(conflict_map, conflict, death, fun="count", na.rm=TRUE) 

conflict_count
conflict_death

for (d in unique(data_raw_conflict$date_conflict) ){
  data_loop<-filter(data_raw_conflict, date_conflict==d)
  conflict_map_loop<-cbind(data_loop$longitude, data_loop$latitude)
  death_loop<-data_loop$fatalities
  conflict_death_loop<-rasterize(conflict_map_loop, conflict, death_loop, fun=sum, na.rm=TRUE) 
  names(conflict_death_loop)<-d
  add(conflict_death)<- conflict_death_loop
  conflict_count_loop<-rasterize(conflict_map_loop, conflict, death_loop, fun="count", na.rm=TRUE) 
  names(conflict_count_loop)<-d
  add(conflict_count)<- conflict_count_loop
}


conflict_death
conflict_count

names(conflict_death)
names(conflict_count)

conflict_count<-subset(conflict_count, 2:77)
conflict_death<-subset(conflict_death, 2:77)


conflict_death<-crop(conflict_death, West, snap="in", mask=T)
conflict_count<-crop(conflict_count, West, snap="in", mask=T)
conflict_death
conflict_count
temp_crs<-conflict_count

plot(conflict_count, y=17)
plot(conflict_death, y=17)
plot(conflict_death, y=14)

writeRaster(conflict_count, 'datasets/final_temp/conflict_count.tif', overwrite=TRUE)
writeRaster(conflict_death, 'datasets/final_temp/conflict_death.tif', overwrite=TRUE)


####_____________1.3 Air pollution PM2.5 from AGAP_________________________________#####

time_path<-paste0('datasets/PM25/Monthly/')
all_pm25<-list.files(time_path, full.names=TRUE, pattern=".nc$") #
length(all_pm25)  #76 months
print(all_pm25[5])

pm25_rast <- rast('datasets/PM25/Monthly/V5GL03.HybridPM25c_0p10.Global.200301-200301.nc')
for (i in 2:76){
  pm25_ <- rast(all_pm25[i])
  print(all_pm25[i]) #checking the order
  add(pm25_rast)<- pm25_
}
sources(pm25_rast) #checking the order
names(pm25_rast)<-months

plot(pm25_rast, y=2)
pm25_rast
summary(pm25_rast)
dim(pm25_rast)

pm25_rast<-aggregate(pm25_rast, fact=5, fun="mean", na.rm=TRUE)
pm25_rast
pm25_rast<-disagg(pm25_rast, fact=2)
pm25_rast
dim(pm25_rast)
#Cropping
pm25_wa<-crop(pm25_rast, West, snap="in", mask=T)
pm25_wa
dim(pm25_wa)
plot(pm25_wa, y=11)
pm25_wa<-project(pm25_wa, temp_crs, method = "bilinear")
writeRaster(pm25_wa, 'datasets/final_temp/pm25_rast_2.tif', overwrite=TRUE)
####_____________1.4 Dust data set MIDAS___________________________________________####

midas_final_path<-paste0('datasets/Dust/midas')
all_midas_final<-list.files(midas_final_path, full.names=TRUE, pattern=".tif$") #
length(all_midas_final)
print(all_midas_final[1])

midas<- rast('datasets/Dust/midas/midas_wa_2003.tif')
midas

plot(midas, y=11, main="Dust DoD 2003-11", colNA="black")
##loop
for (i in 2:7){
  midas_<-rast(all_midas_final[i])
  print(all_midas_final[i])
  add(midas)<- midas_
}

midas
#removing months after APR-2009
midas<-subset(midas, 1:76)
midas
dim(midas)
names(midas)
names(midas)<-months
plot(midas, y=2)

midas<-aggregate(midas, fact=5, fun="mean", na.rm=TRUE)
midas
midas<-disagg(midas, fact=2)
midas
names(midas)
midas= project(midas, temp_crs, method = "bilinear")
writeRaster(midas, 'datasets/final_temp/midas.tif', overwrite=TRUE)

####_____________1.5 Dust Dust fraction model______________________________________####

dust_path<-paste0('datasets/Dust/Columbia/')
all_dust<-list.files(dust_path, full.names=TRUE, pattern=".nc$") #
print(all_dust[2])
length(all_dust)

dust_model <- rast('datasets/Dust/Columbia/2003.nc')
dust_model_month <- tapp(dust_model, "months", fun=mean, na.rm=TRUE)

for (i in 2:7){
  dust_model_ <- rast(all_dust[i])
  print(all_dust[i])
  dust_model_ <- tapp(dust_model_, "months", fun=mean, na.rm=TRUE)
  add(dust_model_month)<- dust_model_
}

names(dust_model_month)
dust_model_month
plot(dust_model_month, 5)

dust_c_wa<-crop(dust_model_month, West, snap="in", mask=T)
plot(dust_c_wa, y=20)
names(dust_c_wa)<-months
names(dust_c_wa)
dust_c_wa = project(dust_c_wa, temp_crs, method = "bilinear")
writeRaster(dust_c_wa, 'datasets/final_temp/dust.tif', overwrite=TRUE)

####_____________1.6 Night time light data set_____________________________________####

ntl_path<-paste0('datasets/NTL/')
all_ntl<-list.files(ntl_path, full.names=TRUE, pattern=".tif$")
length(all_ntl) #
print(all_ntl[8])

ntl_model <- rast('datasets/NTL/Harmonized_DN_NTL_2003_calDMSP.tif')

for (i in 2:6){
  ntl_model_ <- rast(all_ntl[i])
  print(all_ntl[i])
  add(ntl_model)<- ntl_model_
}

ntl_model
ntl_model<-aggregate(ntl_model, fact=30, fun="max", na.rm=TRUE)
ntl_wa<-crop(ntl_model, West, snap="in", mask=T)
ntl_wa <- project(ntl_wa, temp_crs, method = "bilinear")
ntl_wa

ntl_2009 <- rast('datasets/NTL/Harmonized_DN_NTL_2009_calDMSP.tif')
ntl_2009<-aggregate(ntl_2009, fact=30, fun="max", na.rm=TRUE)
ntl_2009<-crop(ntl_2009, West, snap="in", mask=T)
ntl_2009 <- project(ntl_2009, temp_crs, method = "bilinear")

add(ntl_wa)<- ntl_2009
year1<-c(2003, 2004, 2005, 2006, 2007, 2008,2009)
print(names(ntl_wa))
names(ntl_wa)<-year1
writeRaster(ntl_wa, 'datasets/final_temp/ntl_rast.tif', overwrite=TRUE)

####_____________1.7 Population density____________________________________________####

year2<-c(2000,2005,2010)
pop <- rast('datasets/Population/gpw_v4_population_density_rev11_15_min.nc')
varnames(pop)

#ext(pop) <- ext(conflict_count)
pop_dens_world<-subset(pop, 1:3)
names(pop_dens_world)<-year2
max(pop_dens_world)
plot(pop_dens_world, 3)
pop_wa<-crop(pop, West, snap="in", mask=TRUE)
pop_dens<-subset(pop_wa, 1:3)
names(pop_dens)<-year2
pop_dens <- project(pop_dens, temp_crs, method = "bilinear")

summary(pop_dens_world)# Not reliable they used a sample instead of the whole dataset
summary(pop_dens)
plot(pop_dens, y=3)
writeRaster(pop_dens, 'datasets/final_temp/pop_dens.tif', overwrite=TRUE)



####_____________1.8 Weather data set from Copernicus center_______________________####

###Temperature
temp<-rast('datasets/Weather/copernicus_temp.nc')
temp

temp<-subset(temp, 1:76)
temp
names(temp)
names(temp)<-months
temp<-crop(temp, West, snap="in", mask=T)
temp<- project(temp, temp_crs, method = "bilinear")
writeRaster(temp, 'datasets/final_temp/temp.tif', overwrite=TRUE)

###Precipitation and wind_direction 
copernicus<-rast('datasets/wind_direction/copernicus.nc')
copernicus
#u-wind
u_wind<-subset(copernicus, 1:76)
u_wind
names(u_wind)
names(u_wind)<-months
u_wind<-crop(u_wind, West, snap="in", mask=T)
u_wind<- project(u_wind, temp_crs, method = "bilinear")
writeRaster(u_wind, 'datasets/final_temp/u_wind.tif', overwrite=TRUE)

#v-wind
v_wind<-subset(copernicus, 85:160)
v_wind
names(v_wind)
names(v_wind)<-months
v_wind<-crop(v_wind, West, snap="in", mask=T)
v_wind<- project(v_wind, temp_crs, method = "bilinear")
plot(v_wind, 2)
plot(v_wind, 6)
writeRaster(v_wind, 'datasets/final_temp/v_wind.tif', overwrite=TRUE)

#si_10
si_10<-subset(copernicus, 169:244)
si_10
names(si_10)
names(si_10)<-months
si_10<-crop(si_10, West, snap="in", mask=T)
si_10<- project(si_10, temp_crs, method = "bilinear")
writeRaster(si_10, 'datasets/final_temp/si_10.tif', overwrite=TRUE)

#tp
tp<-subset(copernicus, 253:328)
tp
names(tp)
names(tp)<-months
tp<-crop(tp, West, snap="in", mask=T)
tp<- project(tp, temp_crs, method = "bilinear")
plot(tp, 2)
plot(tp, 6)
writeRaster(tp, 'datasets/final_temp/tp.tif', overwrite=TRUE)

####_____________1.9 West Africa Map_______________________________________________####


admin1<-gadm(c("NGA","CPV", "GHA","BEN","TGO","GMB","BFA","GNB","LBR","MLI","MRT","NER","SEN","SLE","CIV","GIN"),1,'datasets/population')
admin2<-gadm(c("NGA", "GHA","BEN","TGO","GMB","BFA","GNB","LBR","MLI","MRT","NER","SEN","SLE","CIV","GIN"),2,'datasets/population')
ethnic<-vect("datasets/ethnic/Ethnic_Groups_in_Africa.shp")

b<-rast(ext(temp_crs), resolution=0.25, crs="EPSG:4326")
country<-rasterize(West, b, field="GID_0")
admin_1<-rasterize(admin1, b, field="NAME_1")
admin_2<-rasterize(admin2, b, field="NAME_2")
ethnic<-rasterize(ethnic, b, field="Ethnic_g")

plot(country)
plot(admin_2)
plot(admin_1)
plot(ethnic)

writeRaster(ethnic, 'datasets/final_temp/ethnic.tif', overwrite=TRUE)
writeRaster(country, 'datasets/final_temp/country.tif', overwrite=TRUE)
writeRaster(admin_1, 'datasets/final_temp/admin_1.tif', overwrite=TRUE)
writeRaster(admin_2, 'datasets/final_temp/admin_2.tif', overwrite=TRUE)


####_____________1.10 Crop land data set___________________________________________####

cropland<-rast('datasets/agg/cropland/af_cropland.tif')
plot(cropland, colNA="red")

cropland
cropland<-aggregate(cropland, fact=3, fun="mean", na.rm=TRUE)
cropland<-crop(cropland, West, snap="in", mask=T)
cropland <- project(cropland, temp_crs, method = "bilinear")

plot(cropland)
writeRaster(cropland, 'datasets/final_temp/cropland.tif', overwrite=TRUE)


####_____________1.11 Violence against civilians sub-data set from ACLED___________####

table(data_raw_conflict$type)
data_vac <- filter(data_raw_conflict, type=="Violence against civilians")
table(data_vac$type)
#creating a initial data points
data_2003 <- filter(data_vac, date_conflict=="2003-01")
table(data_2003$date)
sum(data_2003[, "fatalities"])
vac_map<-cbind(data_2003$longitude, data_2003$latitude)
summary(vac_map)
dim(vac_map)
death<-data_2003$fatalities
sum(death)


vac<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")

vac_death<-rasterize(vac_map, vac, death, fun=sum, na.rm=TRUE) 
vac_count<-rasterize(vac_map, vac, death, fun="count", na.rm=TRUE) 

vac_count
vac_death

for (d in unique(data_vac$date_conflict) ){
  data_loop<-filter(data_vac, date_conflict==d)
  vac_map_loop<-cbind(data_loop$longitude, data_loop$latitude)
  death_loop<-data_loop$fatalities
  vac_death_loop<-rasterize(vac_map_loop, vac, death_loop, fun=sum, na.rm=TRUE) 
  names(vac_death_loop)<-d
  add(vac_death)<- vac_death_loop
  vac_count_loop<-rasterize(vac_map_loop, vac, death_loop, fun="count", na.rm=TRUE) 
  names(vac_count_loop)<-d
  add(vac_count)<- vac_count_loop
}

vac_death
vac_count

names(vac_death)
names(vac_count )

vac_count<-subset(vac_count, 2:77)
vac_death<-subset(vac_death, 2:77)

vac_count<-crop(vac_count, West, snap="in", mask=T)
vac_death<-crop(vac_death, West, snap="in", mask=T)

plot(vac_count, y=3)
plot(vac_death, y=3)

writeRaster(vac_count, 'datasets/final_temp/vac_count.tif', overwrite=TRUE)
writeRaster(vac_death, 'datasets/final_temp/vac_death.tif', overwrite=TRUE)


####_____________1.12 Battles sub-data set from ACLED______________________________####

table(data_raw_conflict$type)
data_batt <- filter(data_raw_conflict, type=="Battles")
table(data_batt$type)
#creating a initial data points
data_2003 <- filter(data_batt, date_conflict=="2003-01")
table(data_2003$date)
sum(data_2003[, "fatalities"])
batt_map<-cbind(data_2003$longitude, data_2003$latitude)
summary(batt_map)
dim(batt_map)
death<-data_2003$fatalities
sum(death)


batt<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")

batt_death<-rasterize(batt_map, batt, death, fun=sum, na.rm=TRUE) 
batt_count<-rasterize(batt_map, batt, death, fun="count", na.rm=TRUE) 

batt_count
batt_death

for (d in unique(data_batt$date_conflict) ){
  data_loop<-filter(data_batt, date_conflict==d)
  batt_map_loop<-cbind(data_loop$longitude, data_loop$latitude)
  death_loop<-data_loop$fatalities
  batt_death_loop<-rasterize(batt_map_loop, batt, death_loop, fun=sum, na.rm=TRUE) 
  names(batt_death_loop)<-d
  add(batt_death)<- batt_death_loop
  batt_count_loop<-rasterize(batt_map_loop, batt, death_loop, fun="count", na.rm=TRUE) 
  names(batt_count_loop)<-d
  add(batt_count)<- batt_count_loop
}

batt_death
batt_count

names(batt_death)
names(batt_count )

batt_count<-subset(batt_count, 2:77)
batt_death<-subset(batt_death, 2:77)

batt_count<-crop(batt_count, West, snap="in", mask=T)
batt_death<-crop(batt_death, West, snap="in", mask=T)

plot(batt_count, y=3)
plot(batt_death, y=3)

writeRaster(batt_count, 'datasets/final_temp/batt_count.tif', overwrite=TRUE)
writeRaster(batt_death, 'datasets/final_temp/batt_death.tif', overwrite=TRUE)

####_____________1.13 Riots sub-data set from ACLED________________________________####



table(data_raw_conflict$type)
data_riot <- filter(data_raw_conflict, type=="Riots" | type=="Protests")
table(data_riot$type)
#creating a initial data points
data_2003 <- filter(data_riot, date_conflict=="2003-01")
table(data_2003$date)
sum(data_2003[, "fatalities"])
riot_map<-cbind(data_2003$longitude, data_2003$latitude)
summary(riot_map)
dim(riot_map)
death<-data_2003$fatalities
sum(death)


riot<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")

riot_death<-rasterize(riot_map, riot, death, fun=sum, na.rm=TRUE) 
riot_count<-rasterize(riot_map, riot, death, fun="count", na.rm=TRUE) 

riot_count
riot_death

for (d in unique(data_riot$date_conflict) ){
  data_loop<-filter(data_riot, date_conflict==d)
  riot_map_loop<-cbind(data_loop$longitude, data_loop$latitude)
  death_loop<-data_loop$fatalities
  riot_death_loop<-rasterize(riot_map_loop, riot, death_loop, fun=sum, na.rm=TRUE) 
  names(riot_death_loop)<-d
  add(riot_death)<- riot_death_loop
  riot_count_loop<-rasterize(riot_map_loop, riot, death_loop, fun="count", na.rm=TRUE) 
  names(riot_count_loop)<-d
  add(riot_count)<- riot_count_loop
}
print(unique(data_riot$date_conflict))
riot_death
riot_count

names(riot_death)
names(riot_count )

riot_count<-subset(riot_count, 2:77)
riot_death<-subset(riot_death, 2:77)

riot_count<-crop(riot_count, West, snap="in", mask=T)
riot_death<-crop(riot_death, West, snap="in", mask=T)

plot(riot_count, y=3)
plot(riot_death, y=3)

writeRaster(riot_count, 'datasets/final_temp/riot_count.tif', overwrite=TRUE)
writeRaster(riot_death, 'datasets/final_temp/riot_death.tif', overwrite=TRUE)


####_____________1.14 Alternative conflict data set from SCAD______________________####

#Calling raw data
scad_conflict<-read.csv("datasets/Conflict/SCAD/SCAD2018Africa_Final.csv")
head(scad_conflict)
scad_conflict<-filter(scad_conflict, styr<2010 & styr>2002)
table(scad_conflict$styr)
table(scad_conflict$country)
scad_conflict<-data.frame(scad_conflict$latitude, scad_conflict$longitude, 
                          scad_conflict$startdate, scad_conflict$duration, 
                          scad_conflict$styr, scad_conflict$escalation, scad_conflict$eventid, scad_conflict$countryname, scad_conflict$etype)
scad_conflict<-mutate(scad_conflict, count=1)
head(scad_conflict)
names(scad_conflict)[1]="latitude"
names(scad_conflict)[2]="longitude"
names(scad_conflict)[3]="start_date"
names(scad_conflict)[4]="duration"
names(scad_conflict)[5]="year"
names(scad_conflict)[6]="escalation"
names(scad_conflict)[7]="id"
names(scad_conflict)[8]="country"
names(scad_conflict)[9]="type"
#write.csv(scad_conflict, "datasets/final_temp/conflict_check.csv", row.names=FALSE)
head(scad_conflict)
summary(scad_conflict)
#managing some dates


scad_conflict<-mutate(scad_conflict, date_conflict=format(as.Date(scad_conflict$start_date), "%Y-%m"))
scad_conflict<-mutate(scad_conflict, date=as.Date(start_date))
scad_conflict<-filter(scad_conflict, date_conflict<="2009-04")
head(scad_conflict)
str(scad_conflict)
scad_conflict<-scad_conflict[order(scad_conflict$date),]
print(unique(scad_conflict$date_conflict))

#creating a initial data points
data_2003 <- filter(scad_conflict, date_conflict=="2003-01")
table(data_2003$date)
sum(data_2003[, "duration"]) #matched with xls
conflict_map<-cbind(data_2003$longitude, data_2003$latitude)
summary(conflict_map)
dim(conflict_map)
death<-data_2003$duration
sum(death)

conflict<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")

conflict_death<-rasterize(conflict_map, conflict, death, fun=sum, na.rm=TRUE) 
conflict_count<-rasterize(conflict_map, conflict, death, fun="count", na.rm=TRUE) 

conflict_count
conflict_death

for (d in unique(scad_conflict$date_conflict) ){
  data_loop<-filter(scad_conflict, date_conflict==d)
  conflict_map_loop<-cbind(data_loop$longitude, data_loop$latitude)
  death_loop<-data_loop$duration
  conflict_death_loop<-rasterize(conflict_map_loop, conflict, death_loop, fun=sum, na.rm=TRUE) 
  names(conflict_death_loop)<-d
  add(conflict_death)<- conflict_death_loop
  conflict_count_loop<-rasterize(conflict_map_loop, conflict, death_loop, fun="count", na.rm=TRUE) 
  names(conflict_count_loop)<-d
  add(conflict_count)<- conflict_count_loop
}

conflict_death
conflict_count

names(conflict_death)
names(conflict_count)

conflict_count<-subset(conflict_count, 2:77)
conflict_death<-subset(conflict_death, 2:77)

conflict_count<-crop(conflict_count, West, snap="in", mask=T)
conflict_death<-crop(conflict_death, West, snap="in", mask=T)

plot(conflict_count, y=3)
plot(conflict_death, y=3)

writeRaster(conflict_count, 'datasets/final_temp/scad_count.tif', overwrite=TRUE)
writeRaster(conflict_death, 'datasets/final_temp/scad_duration.tif', overwrite=TRUE)

####_____________1.15 Alternative conflict data set from UCDP______________________####

#Calling raw data
#ucdp_conflict<-read.csv("datasets/Conflict/UCDP/GEDEvent_v23_1.csv")
#head(ucdp_conflict)
#ucdp_conflict<-filter(ucdp_conflict, year<2010 & year>2002)
#write.csv(ucdp_conflict, "datasets/Conflict/UCDP/ucdp_period.csv", row.names=FALSE)
ucdp_conflict<-read.csv("datasets/Conflict/UCDP/ucdp_period.csv")
table(ucdp_conflict$year)
table(ucdp_conflict$country)
ucdp_conflict<-data.frame(ucdp_conflict$latitude, ucdp_conflict$longitude, 
                          ucdp_conflict$date_start, ucdp_conflict$deaths_a, 
                          ucdp_conflict$year, ucdp_conflict$type_of_violence, ucdp_conflict$date_end, ucdp_conflict$country)
ucdp_conflict<-mutate(ucdp_conflict, count=1)
head(ucdp_conflict)
names(ucdp_conflict)[1]="latitude"
names(ucdp_conflict)[2]="longitude"
names(ucdp_conflict)[3]="start_date"
names(ucdp_conflict)[4]="fatalities"
names(ucdp_conflict)[5]="year"
names(ucdp_conflict)[6]="type"
names(ucdp_conflict)[7]="date_end"
names(ucdp_conflict)[8]="country"

#write.csv(ucdp_conflict, "datasets/final_temp/conflict_check.csv", row.names=FALSE)
head(ucdp_conflict)
summary(ucdp_conflict)
#managing some dates

ucdp_conflict<-mutate(ucdp_conflict, date_conflict=format(as.Date(ucdp_conflict$start_date), "%Y-%m"))
ucdp_conflict<-mutate(ucdp_conflict, date=as.Date(start_date))
ucdp_conflict<-filter(ucdp_conflict, date_conflict<="2009-04")
head(ucdp_conflict)
str(ucdp_conflict)
ucdp_conflict<-ucdp_conflict[order(ucdp_conflict$date),]
print(unique(ucdp_conflict$date_conflict))

#creating a initial data points
data_2003 <- filter(ucdp_conflict, date_conflict=="2003-01")
table(data_2003$date)
sum(data_2003[, "fatalities"])
conflict_map<-cbind(data_2003$longitude, data_2003$latitude)
summary(conflict_map)
dim(conflict_map)
death<-data_2003$fatalities
sum(death)


conflict<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")

conflict_death<-rasterize(conflict_map, conflict, death, fun=sum, na.rm=TRUE) 
conflict_count<-rasterize(conflict_map, conflict, death, fun="count", na.rm=TRUE) 

conflict_count
conflict_death


for (d in unique(ucdp_conflict$date_conflict) ){
  data_loop<-filter(ucdp_conflict, date_conflict==d)
  conflict_map_loop<-cbind(data_loop$longitude, data_loop$latitude)
  death_loop<-data_loop$fatalities
  conflict_death_loop<-rasterize(conflict_map_loop, conflict, death_loop, fun=sum, na.rm=TRUE) 
  names(conflict_death_loop)<-d
  add(conflict_death)<- conflict_death_loop
  conflict_count_loop<-rasterize(conflict_map_loop, conflict, death_loop, fun="count", na.rm=TRUE) 
  names(conflict_count_loop)<-d
  add(conflict_count)<- conflict_count_loop
}

conflict_death
conflict_count

names(conflict_death)
names(conflict_count)

conflict_count<-subset(conflict_count, 2:77)
conflict_death<-subset(conflict_death, 2:77)

conflict_count<-crop(conflict_count, West, snap="in", mask=T)
conflict_death<-crop(conflict_death, West, snap="in", mask=T)

plot(conflict_count, y=3)
plot(conflict_death, y=3)

writeRaster(conflict_count, 'datasets/final_temp/ucdp_count.tif', overwrite=TRUE)
writeRaster(conflict_death, 'datasets/final_temp/ucdp_death.tif', overwrite=TRUE)

####_______2.CONVERTING RASTER DATA INTO A DATA FRAME______________________________####

#I convert 3 dimensional raster data into 2 dimensional data frame and merge all 
# different data sets into a single data frame that can be directly used for empirical analysis and figures.

####_____________2.1 PM2.5_________________________________________________________####

pm25 <- rast('datasets/final_temp/pm25_rast_2.tif')
pm25
plot(pm25, 3)
pm25_df<-as.data.frame(pm25,na.rm=TRUE, xy=TRUE)
dim(pm25_df)
pm25_df <- reshape(pm25_df,  varying=c(names(pm25_df[3:78])),
                   v.names="pm25",
                   timevar= "months",
                   times = c(names(pm25_df[3:78])),
                   direction = "long")
dim(pm25_df)
summary(pm25_df)
table(pm25_df$months)
rownames(pm25_df)<-1:nrow(pm25_df)

####_____________2.2 Intensity of conflict_________________________________________####
conflict_d <- rast('datasets/final_temp/conflict_death.tif')
conflict_d
plot(conflict_d,3)
conflict_d_df<-as.data.frame(conflict_d, na.rm=FALSE, xy=TRUE)
dim(conflict_d_df)
conflict_d_df<- reshape(conflict_d_df,  varying=c(names(conflict_d_df[3:78])),
                        v.names="cdeath",
                        timevar= "months",
                        times = c(names(conflict_d_df[3:78])),
                        direction = "long")
dim(conflict_d_df)
table(conflict_d_df$months)
summary(conflict_d_df)
conflict_d_df["cdeath"][is.na(conflict_d_df["cdeath"])]<-0
str(conflict_d_df)

pm25_df$x<-round(pm25_df$x,3)
pm25_df$y<-round(pm25_df$y,3)
conflict_d_df$x<-round(conflict_d_df$x,3)
conflict_d_df$y<-round(conflict_d_df$y,3)

data_merge<-inner_join(pm25_df, conflict_d_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id.x,id.y))
summary(data_merge)
dim(data_merge)

####_____________2.3 Outbreak of conflict__________________________________________####

conflict_c <- rast('datasets/final_temp/conflict_count.tif')
conflict_c
plot(conflict_c, y=3)

conflict_c_df<-as.data.frame(conflict_c, na.rm=FALSE, xy=TRUE)
dim(conflict_c_df)
conflict_c_df<- reshape(conflict_c_df,  varying=c(names(conflict_c_df[3:78])),
                        v.names="ccount",
                        timevar= "months",
                        times = c(names(conflict_c_df[3:78])),
                        direction = "long")
dim(conflict_c_df)
summary(conflict_c_df)
conflict_c_df<-mutate(conflict_c_df, conflict=ifelse(ccount>=0,1,0))
conflict_c_df["ccount"][is.na(conflict_c_df["ccount"])]<-0
conflict_c_df["conflict"][is.na(conflict_c_df["conflict"])]<-0
summary(conflict_c_df)
conflict_c_df$x<-round(conflict_c_df$x,3)
conflict_c_df$y<-round(conflict_c_df$y,3)

data_merge<-inner_join(data_merge, conflict_c_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)
dim(data_merge)
####_____________2.4 Dust fraction (model)_________________________________________####

dust<- rast('datasets/final_temp/dust.tif')
dust
plot(dust,2)
dust_df<-as.data.frame(dust, na.rm=FALSE, xy=TRUE)
dim(dust_df)
dust_df<- reshape(dust_df,  varying=c(names(dust_df[3:78])),
                  v.names="dust",
                  timevar= "months",
                  times = c(names(dust_df[3:78])),
                  direction = "long")

dim(dust_df)
dust_df["dust"][is.na(dust_df["dust"])]<-0
summary(dust_df)
dust_df$x<-round(dust_df$x,3)
dust_df$y<-round(dust_df$y,3)

data_merge<-inner_join(data_merge, dust_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)
dim(data_merge)

####_____________2.5 DoD MIDAS_____________________________________________________####

midas<- rast('datasets/final_temp/midas.tif')
midas
plot(midas,2)

midas_df<-as.data.frame(midas, na.rm=FALSE, xy=TRUE)
dim(midas_df)
midas_df<- reshape(midas_df,  varying=c(names(midas_df[3:78])),
                   v.names="midas",
                   timevar= "months",
                   times = c(names(midas_df[3:78])),
                   direction = "long")
dim(midas_df)
summary(midas_df)
midas_df$x<-round(midas_df$x,3)
midas_df$y<-round(midas_df$y,3)

data_merge<-inner_join(data_merge, midas_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)
dim(data_merge)

####_____________2.6 Weather variables_____________________________________________#### 

#temp
temp<- rast('datasets/final_temp/temp.tif')
temp
plot(temp, y=22)

temp_df<-as.data.frame(temp, na.rm=FALSE, xy=TRUE)
dim(temp_df)
temp_df<- reshape(temp_df,  varying=c(names(temp_df[3:78])),
                  v.names="temp",
                  timevar= "months",
                  times = c(names(temp_df[3:78])),
                  direction = "long")
summary(temp_df)
dim(temp_df)
temp_df$x<-round(temp_df$x,3)
temp_df$y<-round(temp_df$y,3)
data_merge<-inner_join(data_merge, temp_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)
dim(data_merge)

#_uwind

u_wind<- rast('datasets/final_temp/u_wind.tif')
u_wind
plot(u_wind,2)

u_wind_df<-as.data.frame(u_wind, na.rm=FALSE, xy=TRUE)
dim(u_wind_df)
u_wind_df<- reshape(u_wind_df,  varying=c(names(u_wind_df[3:78])),
                    v.names="u10",
                    timevar= "months",
                    times = c(names(u_wind_df[3:78])),
                    direction = "long")
dim(u_wind_df)

u_wind_df$x<-round(u_wind_df$x,3)
u_wind_df$y<-round(u_wind_df$y,3)
data_merge<-inner_join(data_merge, u_wind_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))

head(data_merge)
summary(data_merge)
dim(data_merge)


#_v_wind
v_wind<- rast('datasets/final_temp/v_wind.tif')
plot(v_wind, y=3)
dim(v_wind)

v_wind_df<-as.data.frame(v_wind, na.rm=FALSE, xy=TRUE)
dim(v_wind_df)
v_wind_df<- reshape(v_wind_df,  varying=c(names(v_wind_df[3:78])),
                    v.names="v10",
                    timevar= "months",
                    times = c(names(v_wind_df[3:78])),
                    direction = "long")
dim(v_wind_df)
v_wind_df$x<-round(v_wind_df$x,3)
v_wind_df$y<-round(v_wind_df$y,3)

data_merge<-inner_join(data_merge, v_wind_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)

#_Wind speed si_10
si_10<- rast('datasets/final_temp/si_10.tif')
plot(si_10, y=8)
dim(si_10)

si_10_df<-as.data.frame(si_10, na.rm=FALSE, xy=TRUE)
dim(si_10_df)
si_10_df<- reshape(si_10_df,  varying=c(names(si_10_df[3:78])),
                   v.names="si_10",
                   timevar= "months",
                   times = c(names(si_10_df[3:78])),
                   direction = "long")
dim(si_10_df)
si_10_df$x<-round(si_10_df$x,3)
si_10_df$y<-round(si_10_df$y,3)

data_merge<-inner_join(data_merge, si_10_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)

#_Precipitation
tp<- rast('datasets/final_temp/tp.tif')
tp 
plot(tp,2)
plot(tp,7)

tp_df<-as.data.frame(tp, na.rm=FALSE, xy=TRUE)
dim(tp_df)
tp_df<- reshape(tp_df,  varying=c(names(tp_df[3:78])),
                v.names="tp",
                timevar= "months",
                times = c(names(tp_df[3:78])),
                direction = "long")
dim(tp_df)
tp_df$x<-round(tp_df$x,3)
tp_df$y<-round(tp_df$y,3)

data_merge<-inner_join(data_merge, tp_df,  by=c("months", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)
dim(data_merge)

####_____________2.7 Population density____________________________________________####

pop<- rast('datasets/final_temp/pop_dens.tif')
data_2003<-read.csv("datasets/Conflict/ACLED/data_2003.csv")

#creation of filler raster
pop_temp<-rast(resolution=c(0.25, 0.25), crs="EPSG:4326")
conflict_map<-cbind(data_2003$latitude, data_2003$longitude)
death<-NA # converting it to the raster frame I had from FLEXdataframe
pop_gap<-rasterize(conflict_map, pop_temp, death, fun=sum, na.rm=FALSE)
pop_gap<-crop(pop_gap, West, snap="in", mask=T)
summary(pop_gap)


pop_2000<-subset(pop, 1:1)
add(pop_2000)<- pop_gap
add(pop_2000)<- pop_gap
add(pop_2000)<- pop_gap
add(pop_2000)<- pop_gap
summary(pop_2000)
pop_2005<-subset(pop, 2:2) 
add(pop_2000)<- pop_2005
add(pop_2000)<- pop_gap
add(pop_2000)<- pop_gap
add(pop_2000)<- pop_gap
add(pop_2000)<- pop_gap
pop_2010<-subset(pop, 3:3) 
add(pop_2000)<- pop_2010

names(pop_2000)
summary(pop_2000)
year2<-c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010)
names(pop_2000)<-year2

pop_df<-as.data.frame(pop_2000, na.rm=FALSE, xy=TRUE)
summary(pop_df)
pop_df<- reshape(pop_df,  varying=c(names(pop_df[3:13])),
                 v.names="pop",
                 timevar= "year",
                 times = c(names(pop_df[3:13])),
                 direction = "long")
summary(pop_df)

pop_df["pop"][is.na(pop_df["pop"])]<-999
write.csv(pop_df, "datasets/final_temp/pop_df.csv", row.names=FALSE)
pop_df<-read.csv("datasets/final_temp/pop_df_2.csv")
describeBy(pop_df$pop, group=pop_df$year)

pop_df$x<-round(pop_df$x,3)
pop_df$y<-round(pop_df$y,3)
summary(pop_df)
pop_df$year<-as.numeric(pop_df$year)

data_merge<-mutate(data_merge, date=ym(months))
data_merge<-mutate(data_merge, year=as.numeric(format(date, "%Y")))
head(pop_df)

data_merge<-left_join(data_merge, pop_df,  by=c("year", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)
dim(data_merge)

####_____________2.8 Nighttime light_______________________________________________####
ntl<- rast('datasets/final_temp/ntl_rast.tif')
plot(ntl,5)
summary(ntl)
dim(ntl)

ntl_df<-as.data.frame(ntl, na.rm=FALSE, xy=TRUE)
dim(ntl_df)
ntl_df<- reshape(ntl_df,  varying=c(names(ntl_df[3:9])),
                 v.names="ntl",
                 timevar= "year",
                 times = c(names(ntl_df[3:9])),
                 direction = "long")
dim(ntl_df)
ntl_df$x<-round(ntl_df$x,3)
ntl_df$y<-round(ntl_df$y,3)
summary(ntl_df)
ntl_df$year<-as.numeric(ntl_df$year)

data_merge<-inner_join(data_merge, ntl_df,  by=c("year", "x", "y"))
data_merge<-subset(data_merge, select=-c(id))
summary(data_merge)


####_____________2.9 Ethnicity map_________________________________________________####
ethnic<- rast('datasets/final_temp/ethnic.tif')
dim(ethnic)
ethnic
names(ethnic)<-"ethnic"
plot(ethnic, colNA="red")

ethnic<-as.data.frame(ethnic, na.rm=FALSE, xy=TRUE)
summary(ethnic)
ethnic$x<-round(ethnic$x,3)
ethnic$y<-round(ethnic$y,3)
table(ethnic$ethnic)
data_merge<-inner_join(data_merge, ethnic,  by=c("x", "y"))
summary(data_merge)

####_____________2.10 Map of administration level___________________________________####
#country
country<- rast('datasets/final_temp/country.tif')
names(country)<-"ccode"
plot(country, colNA="red")
country

country<-as.data.frame(country, na.rm=FALSE, xy=TRUE)
summary(country)
dim(country)
country$x<-round(country$x,3)
country$y<-round(country$y,3)
table(country$ccode)
data_merge<-inner_join(data_merge, country,  by=c("x", "y"))
summary(data_merge)

#admin1
admin1<- rast('datasets/final_temp/admin_1.tif')
names(admin1)<-"admin1"
plot(admin1, colNA="red")
dim(admin1)

admin1<-as.data.frame(admin1, na.rm=FALSE, xy=TRUE)
summary(admin1)
admin1$x<-round(admin1$x,3)
admin1$y<-round(admin1$y,3)

data_merge<-inner_join(data_merge, admin1,  by=c("x", "y"))
summary(data_merge)
nrow(data_merge)

#admin2
admin2<- rast('datasets/final_temp/admin_2.tif')
names(admin2)<-"admin2"
plot(admin2, colNA="red")

admin2<-as.data.frame(admin2, na.rm=FALSE, xy=TRUE)
summary(admin2)
admin2$x<-round(admin2$x,3)
admin2$y<-round(admin2$y,3)

data_merge<-inner_join(data_merge, admin2,  by=c("x", "y"))
summary(data_merge)
#Ethnic, country, admin1, and admin2 are the factor variable
nrow(data_merge)
str(data_merge)

####_____________2.11 Alternative data set SCAD____________________________________####
#scad_death


scad_d <- rast('datasets/final_temp/scad_duration.tif')
scad_d
plot(scad_d, y=3)

scad_d_df<-as.data.frame(scad_d, na.rm=FALSE, xy=TRUE)
dim(scad_d_df)
scad_d_df<- reshape(scad_d_df,  varying=c(names(scad_d_df[3:78])),
                    v.names="scad_duration",
                    timevar= "months",
                    times = c(names(scad_d_df[3:78])),
                    direction = "long")
dim(scad_d_df)
table(scad_d_df$months)
summary(scad_d_df)
scad_d_df["scad_duration"][is.na(scad_d_df["scad_duration"])]<-0
str(scad_d_df)
rownames(scad_d_df)<-1:nrow(scad_d_df)


scad_d_df$x<-round(scad_d_df$x,3)
scad_d_df$y<-round(scad_d_df$y,3)
data_merge2<-inner_join(data_merge, scad_d_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)


###___scad_c

scad_c <- rast('datasets/final_temp/scad_count.tif')
plot(scad_c, y=2)
scad_c


scad_c_df<-as.data.frame(scad_c, na.rm=FALSE, xy=TRUE)
dim(scad_c_df)
scad_c_df<- reshape(scad_c_df,  varying=c(names(scad_c_df[3:78])),
                    v.names="scad_count",
                    timevar= "months",
                    times = c(names(scad_c_df[3:78])),
                    direction = "long")
dim(scad_c_df)
scad_c_df<-mutate(scad_c_df, scad=ifelse(scad_count>=0,1,0))
scad_c_df["scad_count"][is.na(scad_c_df["scad_count"])]<-0
scad_c_df["scad"][is.na(scad_c_df["scad"])]<-0
summary(scad_c_df)

scad_c_df$x<-round(scad_c_df$x,3)
scad_c_df$y<-round(scad_c_df$y,3)
data_merge2<-inner_join(data_merge2, scad_c_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
nrow(data_merge2)
####_____________2.12 Alternative data set UCPD____________________________________####
#ucpd_death

ucdp_d <- rast('datasets/final_temp/ucdp_death.tif')
ucdp_d
plot(ucdp_d, y=3)

ucdp_d_df<-as.data.frame(ucdp_d, na.rm=FALSE, xy=TRUE)
dim(ucdp_d_df)
ucdp_d_df<- reshape(ucdp_d_df,  varying=c(names(ucdp_d_df[3:78])),
                    v.names="ucdp_death",
                    timevar= "months",
                    times = c(names(ucdp_d_df[3:78])),
                    direction = "long")
dim(ucdp_d_df)
table(ucdp_d_df$months)
summary(ucdp_d_df)
ucdp_d_df["ucdp_death"][is.na(ucdp_d_df["ucdp_death"])]<-0
str(ucdp_d_df)
rownames(ucdp_d_df)<-1:nrow(ucdp_d_df)

ucdp_d_df$x<-round(ucdp_d_df$x,3)
ucdp_d_df$y<-round(ucdp_d_df$y,3)
data_merge2<-inner_join(data_merge2, ucdp_d_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
dim(data_merge2)

###___ucdp_c

ucdp_c <- rast('datasets/final_temp/ucdp_count.tif')
plot(ucdp_c, y=2)
ucdp_c

ucdp_c_df<-as.data.frame(ucdp_c, na.rm=FALSE, xy=TRUE)
dim(ucdp_c_df)
ucdp_c_df<- reshape(ucdp_c_df,  varying=c(names(ucdp_c_df[3:78])),
                    v.names="ucdp_count",
                    timevar= "months",
                    times = c(names(ucdp_c_df[3:78])),
                    direction = "long")
dim(ucdp_c_df)
ucdp_c_df<-mutate(ucdp_c_df, ucdp=ifelse(ucdp_count>=0,1,0))
ucdp_c_df["ucdp_count"][is.na(ucdp_c_df["ucdp_count"])]<-0
ucdp_c_df["ucdp"][is.na(ucdp_c_df["ucdp"])]<-0
summary(ucdp_c_df)

ucdp_c_df$x<-round(ucdp_c_df$x,3)
ucdp_c_df$y<-round(ucdp_c_df$y,3)
data_merge2<-inner_join(data_merge2, ucdp_c_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
dim(data_merge2)

####_____________2.13 Cropland data set ___________________________________________####
cland<- rast('datasets/final_temp/cropland.tif')
names(cland)<-"cland"
cland
plot(cland, colNA="red")

cland<-as.data.frame(cland, na.rm=FALSE, xy=TRUE)
dim(cland)
summary(cland)
cland$x<-round(cland$x,3)
cland$y<-round(cland$y,3)

data_merge2<-inner_join(data_merge2, cland,  by=c("x", "y"))
summary(data_merge2)
nrow(data_merge2)

####_____________2.14 Conflict sub-data sets VAC___________________________________####
#vac_death

vac_d <- rast('datasets/final_temp/vac_death.tif')
vac_d

vac_d_df<-as.data.frame(vac_d, na.rm=FALSE, xy=TRUE)
dim(vac_d_df)
vac_d_df<- reshape(vac_d_df,  varying=c(names(vac_d_df[3:78])),
                   v.names="vac_death",
                   timevar= "months",
                   times = c(names(vac_d_df[3:78])),
                   direction = "long")
dim(vac_d_df)
table(vac_d_df$months)
summary(vac_d_df)
vac_d_df["vac_death"][is.na(vac_d_df["vac_death"])]<-0
str(vac_d_df)
rownames(vac_d_df)<-1:nrow(vac_d_df)

vac_d_df$x<-round(vac_d_df$x,3)
vac_d_df$y<-round(vac_d_df$y,3)

data_merge2<-inner_join(data_merge2, vac_d_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
dim(data_merge2)

###___vac_c

vac_c <- rast('datasets/final_temp/vac_count.tif')
plot(vac_c, y=2)
vac_c

vac_c_df<-as.data.frame(vac_c, na.rm=FALSE, xy=TRUE)
vac_c_df<- reshape(vac_c_df,  varying=c(names(vac_c_df[3:78])),
                   v.names="vac_count",
                   timevar= "months",
                   times = c(names(vac_c_df[3:78])),
                   direction = "long")
vac_c_df<-mutate(vac_c_df, vac=ifelse(vac_count>=0,1,0))
vac_c_df["vac_count"][is.na(vac_c_df["vac_count"])]<-0
vac_c_df["vac"][is.na(vac_c_df["vac"])]<-0
summary(vac_c_df)

vac_c_df$x<-round(vac_c_df$x,3)
vac_c_df$y<-round(vac_c_df$y,3)
data_merge2<-inner_join(data_merge2, vac_c_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
dim(data_merge2)

#batt_death
batt_d <- rast('datasets/final_temp/batt_death.tif')

batt_d_df<-as.data.frame(batt_d, na.rm=FALSE, xy=TRUE)
batt_d_df<- reshape(batt_d_df,  varying=c(names(batt_d_df[3:78])),
                    v.names="batt_death",
                    timevar= "months",
                    times = c(names(batt_d_df[3:78])),
                    direction = "long")

table(batt_d_df$months)
summary(batt_d_df)
batt_d_df["batt_death"][is.na(batt_d_df["batt_death"])]<-0
str(batt_d_df)
rownames(batt_d_df)<-1:nrow(batt_d_df)

batt_d_df$x<-round(batt_d_df$x,3)
batt_d_df$y<-round(batt_d_df$y,3)
data_merge2<-inner_join(data_merge2, batt_d_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
dim(data_merge2)

###___batt_c

batt_c <- rast('datasets/final_temp/batt_count.tif')
plot(batt_c, y=2)

batt_c_df<-as.data.frame(batt_c, na.rm=FALSE, xy=TRUE)
dim(batt_c_df)
batt_c_df<- reshape(batt_c_df,  varying=c(names(batt_c_df[3:78])),
                    v.names="batt_count",
                    timevar= "months",
                    times = c(names(batt_c_df[3:78])),
                    direction = "long")

batt_c_df<-mutate(batt_c_df, batt=ifelse(batt_count>=0,1,0))
batt_c_df["batt_count"][is.na(batt_c_df["batt_count"])]<-0
batt_c_df["batt"][is.na(batt_c_df["batt"])]<-0
summary(batt_c_df)

batt_c_df$x<-round(batt_c_df$x,3)
batt_c_df$y<-round(batt_c_df$y,3)

data_merge2<-inner_join(data_merge2, batt_c_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
summary(data_merge2)
dim(data_merge2)

#riot_death
riot_d <- rast('datasets/final_temp/riot_death.tif')

riot_d_df<-as.data.frame(riot_d, na.rm=FALSE, xy=TRUE)
dim(riot_d_df)
riot_d_df<- reshape(riot_d_df,  varying=c(names(riot_d_df[3:77])),
                    v.names="riot_death",
                    timevar= "months",
                    times = c(names(riot_d_df[3:77])),
                    direction = "long")

table(riot_d_df$months)
summary(riot_d_df)
riot_d_df["riot_death"][is.na(riot_d_df["riot_death"])]<-0
str(riot_d_df)
rownames(riot_d_df)<-1:nrow(riot_d_df)

riot_d_df$x<-round(riot_d_df$x,3)
riot_d_df$y<-round(riot_d_df$y,3)
data_merge2<-left_join(data_merge2, riot_d_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
data_merge2["riot_death"][is.na(data_merge2["riot_death"])]<-0
summary(data_merge2)
dim(data_merge2)

###___riot_c

riot_c <- rast('datasets/final_temp/riot_count.tif')
plot(riot_c, y=2)

riot_c_df<-as.data.frame(riot_c, na.rm=FALSE, xy=TRUE)
dim(riot_c_df)
riot_c_df<- reshape(riot_c_df,  varying=c(names(riot_c_df[3:77])),
                    v.names="riot_count",
                    timevar= "months",
                    times = c(names(riot_c_df[3:77])),
                    direction = "long")

riot_c_df<-mutate(riot_c_df, riot=ifelse(riot_count>=0,1,0))
riot_c_df["riot_count"][is.na(riot_c_df["riot_count"])]<-0
riot_c_df["riot"][is.na(riot_c_df["riot"])]<-0
summary(riot_c_df)

riot_c_df$x<-round(riot_c_df$x,3)
riot_c_df$y<-round(riot_c_df$y,3)
data_merge2<-left_join(data_merge2, riot_c_df,  by=c("months", "x", "y"))
data_merge2<-subset(data_merge2, select=-c(id))
data_merge2["riot_count"][is.na(data_merge2["riot_count"])]<-0
data_merge2["riot"][is.na(data_merge2["riot"])]<-0
summary(data_merge2)
dim(data_merge2)

####_____________2.15 Health expenditure data set__________________________________####
#share of the out of the pocket expenditure

pocket_exp<-read.csv("datasets/health/pocket_exp.csv")
pocket_exp<-as.data.frame(pocket_exp)
head(pocket_exp)
str(pocket_exp)
summary(pocket_exp)

pocket_exp<-mutate(pocket_exp, year=as.numeric(year))
data_merge2<-inner_join(data_merge2, pocket_exp,  by=c("ccode", "year"))
summary(data_merge2)
table(data_merge2$ccode)
dim(data_merge2)

####_____________2.16 Defense expenditure__________________________________________####
#share of the out of the pocket expenditure

defense_exp<-read.csv("datasets/defense/defense.csv")
defense_exp<-as.data.frame(defense_exp)
head(defense_exp)
str(defense_exp)
summary(defense_exp)

defense_exp<-mutate(defense_exp, year=as.numeric(year))
data_merge2<-inner_join(data_merge2, defense_exp,  by=c("ccode", "year"))
summary(data_merge2)
table(data_merge2$ccode)
dim(data_merge2)

####_____________2.17 Transmission:Good governance index___________________________####

cgi<-read.csv("datasets/Institution/cgi/cgi_data.csv")
cgi<-as.data.frame(cgi)
head(cgi, 10)
str(cgi)
summary(cgi)

cgi<-mutate(cgi, year=as.numeric(year))
data_merge2<-inner_join(data_merge2, cgi,  by=c("ccode", "year"))
summary(data_merge2)
table(data_merge2$ccode)
dim(data_merge2)

#write.csv(data_merge2, "datasets/data_check.csv", row.names=FALSE)
####_______3.DATA CLEANING ________________________________________________________####

# In this section, I prepare the merged data set into data visualization and empirical analysis

####_____________3.1 Removing NA's_________________________________________________####
data_merge2<-data_merge2[!is.na(data_merge2$pop),]
dim(data_merge2)
data_merge2<-mutate(data_merge2, midas=ifelse(midas<0,0,midas))
dim(data_merge2)
data_merge2<-data_merge2[!is.na(data_merge2$pm25),]
dim(data_merge2)
data_merge2<-data_merge2[!is.na(data_merge2$midas),]
dim(data_merge2)
####_____________3.2 Creating new variables________________________________________####
data_merge2<-data_merge2 %>% group_by(x, y) %>%
  mutate(grid = cur_group_id())%>% ungroup() %>%
  arrange(grid)
data_merge2<-mutate(data_merge2, month=month(date))

#https://stats.stackexchange.com/questions/24131/using-r-and-plm-to-estimate-fixed-effects-models-that-include-interactions-with
#https://stackoverflow.com/questions/40367855/iv-estimation-with-cluster-robust-standard-errors-using-the-plm-package-in-r
#Exploring the dataset 
str(data_merge2)
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}
data_merge2<-mutate(data_merge2, ihs_vac=ihs(vac_death))
data_merge2<-mutate(data_merge2, ihs_riot=ihs(riot_death))
data_merge2<-mutate(data_merge2, ihs_batt=ihs(batt_death))
data_merge2<-mutate(data_merge2, ihs_ntl=ihs(ntl))
data_merge2<-mutate(data_merge2, ihs_pop=ihs(pop))
data_merge2<-mutate(data_merge2, ihs_cdeath=ihs(cdeath))
data_merge2<-mutate(data_merge2, pm25_mg=pm25/1000)
data_merge2<-mutate(data_merge2, pm25_mg_sq=pm25_mg*pm25_mg)
data_merge2<-mutate(data_merge2, ihs_ucdp=ihs(ucdp_death))
data_merge2<-mutate(data_merge2, ihs_scad=ihs(scad_duration))
data_merge2<-mutate(data_merge2, agg=ifelse(cland>0,1,0))

bins<-rbin_quantiles(data_merge2, y, temp, 5)
data_merge2<-rbin_create(data_merge2,temp, bins)
data_merge2<-data_merge2 %>% rename("temp1"="temp_<_297.476379394531",
                                    "temp2"="temp_<_299.822143554688",
                                    "temp3"="temp_<_302.038024902344",
                                    "temp4"="temp_<_305.165313720703",
                                    "temp5"="temp_>=_305.165313720703")
head(data_merge2)
str(data_merge2)
summary(data_merge2)
####_____________3.3 Calculating averages and creating dummies_____________________####
#urban
data_merge2<-data_merge2 %>% group_by(year)%>% mutate(pop_mean=mean(pop, na.rm=TRUE)) #%>% select( months, crop_share) %>% head()
data_merge2<-mutate(data_merge2, urban=ifelse(pop>pop_mean,1,0))
table(data_merge2$urban, group=data_merge2$year)
describeBy(data_merge2$ccount, group=data_merge2$urban)
rowsum(data_merge2$ccount, data_merge2$urban)

#Health expenditure
data_merge2<-data_merge2 %>% group_by(year)%>% mutate(pocket_mean=mean(pocket)) #%>% select( months, crop_share) %>% head()
data_merge2<-mutate(data_merge2, pocket_high=ifelse(pocket>pocket_mean,1,0))
describeBy(data_merge2$pocket_mean, group=data_merge2$year)

#Defense expenditure
data_merge2<-data_merge2 %>% group_by(year)%>% mutate(defense_mean=mean(defense, na.rm=TRUE)) #%>% select( months, crop_share) %>% head()
data_merge2<-mutate(data_merge2, defense_high=ifelse(defense>defense_mean,1,0))
describeBy(data_merge2$defense_mean, group=data_merge2$year)

#Harmattan
data_merge2<-mutate(data_merge2, harmattan=ifelse(month<4 | month>10, 1,0))
table(data_merge2$harmattan)
describeBy(data_merge2$harmattan, group=data_merge2$month)
describeBy(data_merge2$pop_mean, group=data_merge2$year)

#rule of law
data_merge2<-data_merge2 %>% group_by(year)%>% mutate(rule_mean=mean(rule)) #%>% select( months, crop_share) %>% head()
data_merge2<-mutate(data_merge2, rol_high=ifelse(rule>rule_mean,1,0))
describeBy(data_merge2$rule_mean, group=data_merge2$year)

#ntl
data_merge2<-data_merge2 %>% group_by(year)%>% mutate(ntl_mean=mean(ntl)) #%>% select( months, crop_share) %>% head()
data_merge2<-mutate(data_merge2, ntl_bright=ifelse(ntl>ntl_mean,1,0))
describeBy(data_merge2$ntl_mean, group=data_merge2$year)


write.csv(data_merge2, "datasets/data_01.csv", row.names=FALSE)
#data<-read.csv("datasets/data_01.csv")
describeBy(data$pm25_mg, group=data$year)
describeBy(data$pocket, group=data$ccode)
describeBy(data$defense, group=data$code)
describeBy(data$defense, group=data$year)


