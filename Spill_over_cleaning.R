
#***************** Up in the Air: The Effects of Air Pollution on Conflict in West Africa**********************

#        MASTER THESIS at the Georg-August University of Goettingen 

#                       Ankhbayar Delgerchuluun

#                               RSCRIPT
#                          ------------------
#                            Data Cleaning for 0.5 vs. 0.5

####_______________________README__________________________________________________####

#This R script carries out all pre-analysis process to clean the data in 0.5 vs. 0.5 grid level to 
# check the spillover effect. 
#Thus, it consists of 3 main sections and 22 sub-sections. In the beginning of the main sections, 
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
#crop time dimension between 2003 and 2009, and prject all them with with EGS4386 format. 
####_____________1.1 Defining the time dimension___________________________________####
start_date<-as.Date("2003-01-01")
end_date<-as.Date("2009-4-30")
month1 <- function(x) as.Date(cut(x, "month"))
months <- seq(month1(start_date), month1(end_date), "month")
months<-format(months,"%Y-%m")
head(months)
#Importing a shape file of West African map
West<-gadm(c("NGA","CPV", "GHA","BEN","TGO","GMB","BFA","GNB","LBR","MLI","MRT","NER","SEN","SLE","CIV","GIN"),0,'datasets/population')


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
write.csv(data_raw_conflict, "datasets/spill_temp/conflict_check.csv", row.names=FALSE)
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
write.csv(data_raw_conflict, "datasets/spill_temp/conflict_check.csv", row.names=FALSE)
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


conflict<-rast(resolution=c(0.5, 0.5), crs="EPSG:4326")

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
#year<-as.numeric(format(months, "%Y"))
#conflict_death_mean <- tapp(conflict_death, year, fun=sum, na.rm=TRUE)
#conflict_count_mean <- tapp(conflict_count, year, fun=sum, na.rm=TRUE)
#plot(conflict_death_mean, y=2)
#plot(conflict_count_mean, y=2)

plot(conflict_count, y=17)
plot(conflict_death, y=17)
plot(conflict_death, y=14)

writeRaster(conflict_count, 'datasets/spill_temp/conflict_count.tif', overwrite=TRUE)
writeRaster(conflict_death, 'datasets/spill_temp/conflict_death.tif', overwrite=TRUE)


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
dim(pm25_rast)
#Cropping
pm25_wa<-crop(pm25_rast, West, snap="in", mask=T)
pm25_wa
dim(pm25_wa)
plot(pm25_wa, y=19)
pm25_wa<-project(pm25_wa, temp_crs, method = "bilinear")
writeRaster(pm25_wa, 'datasets/spill_temp/pm25_rast_2.tif', overwrite=TRUE)
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
midas<-crop(midas, West, snap="in", mask=T)
midas<-aggregate(midas, fact=5, fun="mean", na.rm=TRUE)
midas
names(midas)
midas= project(midas, temp_crs, method = "bilinear")
writeRaster(midas, 'datasets/spill_temp/midas.tif', overwrite=TRUE)

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
plot(ntl_model, y=2)
ntl_model<-aggregate(ntl_model, fact=60, fun="max", na.rm=TRUE)
ntl_wa<-crop(ntl_model, West, snap="in", mask=T)
ntl_wa <- project(ntl_wa, temp_crs, method = "bilinear")
ntl_wa

ntl_2009 <- rast('datasets/NTL/Harmonized_DN_NTL_2009_calDMSP.tif')
ntl_2009<-aggregate(ntl_2009, fact=60, fun="max", na.rm=TRUE)
ntl_2009<-crop(ntl_2009, West, snap="in", mask=T)
ntl_2009 <- project(ntl_2009, temp_crs, method = "bilinear")

add(ntl_wa)<- ntl_2009
names(ntl_wa)
year1<-c(2003, 2004, 2005, 2006, 2007, 2008,2009)
print(names(ntl_wa))
names(ntl_wa)<-year1
plot(ntl_wa, y=7)
writeRaster(ntl_wa, 'datasets/spill_temp/ntl_rast.tif', overwrite=TRUE)

####_____________1.7 Population density____________________________________________####

year2<-c(2000,2005,2010)
pop <- rast('datasets/Population/gpw_v4_population_density_rev11_15_min.nc')
varnames(pop)

#ext(pop) <- ext(conflict_count)
pop_dens_world<-subset(pop, 1:3)
names(pop_dens_world)<-year2
max(pop_dens_world)
plot(pop_dens_world, 3)
pop<-aggregate(pop, fact=2, fun="mean", na.rm=TRUE)
pop
pop_wa<-crop(pop, West, snap="in", mask=TRUE)
pop_dens<-subset(pop_wa, 1:3)
names(pop_dens)<-year2
pop_dens <- project(pop_dens, temp_crs, method = "bilinear")

summary(pop_dens_world)# Not reliable they used a sample instead of the whole dataset
summary(pop_dens)
plot(pop_dens, y=3)
writeRaster(pop_dens, 'datasets/spill_temp/pop_dens.tif', overwrite=TRUE)



####_____________1.8 Weather data set from Copernicus center_______________________####

###Temperature
temp<-rast('datasets/Weather/copernicus_temp.nc')
temp
#u-wind
temp<-subset(temp, 1:76)
temp
names(temp)
names(temp)<-months
temp<-aggregate(temp, fact=2, fun="mean", na.rm=TRUE)
temp
temp<-crop(temp, West, snap="in", mask=T)
temp<- project(temp, temp_crs, method = "bilinear")
writeRaster(temp, 'datasets/spill_temp/temp.tif', overwrite=TRUE)

###Precipitation and wind_direction 
copernicus<-rast('datasets/wind_direction/copernicus.nc')
copernicus<-aggregate(copernicus, fact=2, fun="mean", na.rm=TRUE)
copernicus
#u-wind
u_wind<-subset(copernicus, 1:76)
u_wind
names(u_wind)
names(u_wind)<-months
u_wind<-crop(u_wind, West, snap="in", mask=T)
u_wind<- project(u_wind, temp_crs, method = "bilinear")
plot(u_wind,y=2)
writeRaster(u_wind, 'datasets/spill_temp/u_wind.tif', overwrite=TRUE)

#v-wind
v_wind<-subset(copernicus, 85:160)
v_wind
names(v_wind)
names(v_wind)<-months
v_wind<-crop(v_wind, West, snap="in", mask=T)
v_wind<- project(v_wind, temp_crs, method = "bilinear")
plot(v_wind, 2)
plot(v_wind, 6)
writeRaster(v_wind, 'datasets/spill_temp/v_wind.tif', overwrite=TRUE)

#si_10
si_10<-subset(copernicus, 169:244)
si_10
names(si_10)
names(si_10)<-months
si_10<-crop(si_10, West, snap="in", mask=T)
si_10<- project(si_10, temp_crs, method = "bilinear")
writeRaster(si_10, 'datasets/spill_temp/si_10.tif', overwrite=TRUE)

#tp
tp<-subset(copernicus, 253:328)
tp
names(tp)
names(tp)<-months
tp<-crop(tp, West, snap="in", mask=T)
tp<- project(tp, temp_crs, method = "bilinear")
plot(tp, 2)
plot(tp, 6)
writeRaster(tp, 'datasets/spill_temp/tp.tif', overwrite=TRUE)

####_____________1.9 West Africa Map_______________________________________________####


admin1<-gadm(c("NGA","CPV", "GHA","BEN","TGO","GMB","BFA","GNB","LBR","MLI","MRT","NER","SEN","SLE","CIV","GIN"),1,'datasets/population')
admin2<-gadm(c("NGA", "GHA","BEN","TGO","GMB","BFA","GNB","LBR","MLI","MRT","NER","SEN","SLE","CIV","GIN"),2,'datasets/population')
ethnic<-vect("datasets/ethnic/Ethnic_Groups_in_Africa.shp")

b<-rast(ext(temp_crs), resolution=0.5, crs="EPSG:4326")
country<-rasterize(West, b, field="GID_0")
admin_1<-rasterize(admin1, b, field="NAME_1")
admin_2<-rasterize(admin2, b, field="NAME_2")
ethnic<-rasterize(ethnic, b, field="Ethnic_g")

plot(country)
plot(admin_2)
plot(admin_1)
plot(ethnic)

writeRaster(ethnic, 'datasets/spill_temp/ethnic.tif', overwrite=TRUE)
writeRaster(country, 'datasets/spill_temp/country.tif', overwrite=TRUE)
writeRaster(admin_1, 'datasets/spill_temp/admin_1.tif', overwrite=TRUE)
writeRaster(admin_2, 'datasets/spill_temp/admin_2.tif', overwrite=TRUE)


####_______2.CONVERTING RASTER DATA INTO A DATA FRAME______________________________####

#I convert 3 dimensional raster data into 2 dimensional data frame and merge all 
# different data sets into a single data frame that can be directly used for empirical analysis and figures.

####_____________2.1 PM2.5_________________________________________________________####

pm25 <- rast('datasets/spill_temp/pm25_rast_2.tif')
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
conflict_d <- rast('datasets/spill_temp/conflict_death.tif')
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

conflict_c <- rast('datasets/spill_temp/conflict_count.tif')
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
####_____________2.5 DoD MIDAS_____________________________________________________####

midas<- rast('datasets/spill_temp/midas.tif')
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
temp<- rast('datasets/spill_temp/temp.tif')
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

u_wind<- rast('datasets/spill_temp/u_wind.tif')
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
v_wind<- rast('datasets/spill_temp/v_wind.tif')
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
si_10<- rast('datasets/spill_temp/si_10.tif')
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
tp<- rast('datasets/spill_temp/tp.tif')
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

pop<- rast('datasets/spill_temp/pop_dens.tif')
data_2003<-read.csv("datasets/Conflict/ACLED/data_2003.csv")

#creation of filler raster
pop_temp<-rast(resolution=c(0.5, 0.5), crs="EPSG:4326")
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
write.csv(pop_df, "datasets/spill_temp/pop_df_05.csv", row.names=FALSE)
pop_df<-read.csv("datasets/spill_temp/pop_df_2_05.csv")
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
ntl<- rast('datasets/spill_temp/ntl_rast.tif')
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

####_____________2.10 Map of administration level___________________________________####
#country
country<- rast('datasets/spill_temp/country.tif')
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


####_______3.DATA CLEANING ________________________________________________________####

# In this section, I prepare the merged data set into data visualization and empirical analysis

####_____________3.1 Removing NA's_________________________________________________####
data_merge<-data_merge[!is.na(data_merge$pop),]
dim(data_merge)
data_merge<-mutate(data_merge, midas=ifelse(midas<0,0,midas))
dim(data_merge)
data_merge<-data_merge[!is.na(data_merge$pm25),]
dim(data_merge)
data_merge<-data_merge[!is.na(data_merge$midas),]
dim(data_merge)
####_____________3.2 Creating new variables________________________________________####
data_merge<-data_merge %>% group_by(x, y) %>%
  mutate(grid = cur_group_id())%>% ungroup() %>%
  arrange(grid)
data_merge<-mutate(data_merge, month=month(date))

#https://stats.stackexchange.com/questions/24131/using-r-and-plm-to-estimate-fixed-effects-models-that-include-interactions-with
#https://stackoverflow.com/questions/40367855/iv-estimation-with-cluster-robust-standard-errors-using-the-plm-package-in-r
#Exploring the dataset 
str(data_merge)
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}
data_merge<-mutate(data_merge, ihs_ntl=ihs(ntl))
data_merge<-mutate(data_merge, ihs_pop=ihs(pop))
data_merge<-mutate(data_merge, ihs_cdeath=ihs(cdeath))
data_merge<-mutate(data_merge, pm25_mg=pm25/1000)
data_merge<-mutate(data_merge, pm25_mg_sq=pm25_mg*pm25_mg)

####_____________3.3 Calculating averages and creating dummies_____________________####
#urban
data_merge<-data_merge %>% group_by(year)%>% mutate(pop_mean=mean(pop, na.rm=TRUE)) #%>% select( months, crop_share) %>% head()
data_merge<-mutate(data_merge, urban=ifelse(pop>pop_mean,1,0))
table(data_merge$urban, group=data_merge$year)
describeBy(data_merge$ccount, group=data_merge$urban)
rowsum(data_merge$ccount, data_merge$urban)
#Harmattan
data_merge<-mutate(data_merge, harmattan=ifelse(month<4 | month>10, 1,0))
table(data_merge$harmattan)
describeBy(data_merge$harmattan, group=data_merge$month)
describeBy(data_merge$pop_mean, group=data_merge$year)

#ntl
data_merge<-data_merge %>% group_by(year)%>% mutate(ntl_mean=mean(ntl)) #%>% select( months, crop_share) %>% head()
data_merge<-mutate(data_merge, ntl_bright=ifelse(ntl>ntl_mean,1,0))
describeBy(data_merge$ntl_mean, group=data_merge$year)


write.csv(data_merge, "datasets/data_05.csv", row.names=FALSE)



