#############################
### Extract forest cover  ###
#############################

##### Example for 2014, repeat for all years/tiles

## Load final XSS database and extract centroid IDs and points 
df <- read.csv("forest_tiles.csv")         # read in csv data to match forest raster with GPS data 
df <- df[c("cent_id", "fr14", "X", "Y")]   # substitute in column for required tile raster
df <- df %>% filter(!is.na(fr14))          # filter and extract from dataset data to match tile
View(df)

## Load raster data
fr14 <- raster('Forest10N_110E2014.tif')

## Set coordinates
coordinates(df) <- c("X", "Y")

## Extract land cover around houses
for(i in 1:length(buffer_distances)){
  
  ## Extract data within buffer
  fr <- extract(fr14, df, buffer=buffer_distances[i])      # load appropriate raster tile
  
  ## Summarise by proportions
  sm <- lapply(fr, function(x){prop.table(table(x))})
  
  # convert to data frame
  ev <- data.frame(id = rep(df$tile_id, lapply(sm, length)),
                   cover = names(unlist(sm)),
                   percent = unlist(sm), 
                   X = rep(df$X, lapply(sm, length)), 
                   Y = rep(df$Y, lapply(sm, length))
  )
  
  ev$dist <- buffer_distances[i]
  
  # write output
  write.csv(ev, paste0("C:/Users/Kimberly Fornace/Documents/EJ_Pknowlesi_2021/outputs/fr14_", buffer_distances[i], ".csv"))    # substitute tile code
  
}

######################################
### Extract fragmentation metrics  ###
######################################

install.packages('pkgdir')

library(tidyr)
library(gstat)
library(rgeos)
library(scales)
library(landscapemetrics)
library(raster)
library(sp)
library(raster)
library(dplyr)
library(sf)
library(rgdal)
library(raster)
library(lsm)
library(landscapemetrics)
library(ggplot2)
library(SDMTools)

## RTools (installing)
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
Sys.which("make")

install.packages("remotes")
remotes::install_github("matthewkling/ecoclim")

## SDMTools (installing) 
install.packages("remotes")
library(remotes)
remotes::install_version("SDMTools", "1.1-221")

install.packages('rgeos')
library(rgeos)

###### Example for 2014 forest cover and 2014 data points (repeat for all years/tiles)

# read in raster, save by tile code 
fr14 <- raster("forest_R_2014.tif")

# input vector of distances (in meters)
buffer_distances <- c(5000, 10000, 20000)

## Load database, extract IDs and points to match raster tile
df <- read.csv("forest_tiles_2014extra.csv")
df <- df[c("tile_id", "fr14", "X", "Y")]
df <- df %>% filter(!is.na(fr14))

## Set coordinates
coordinates(df) <- c("X", "Y")
proj4string(df) <- CRS("+proj=longlat +datum=WGS84")
df <- spTransform(df, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"))

## If WG84, re-project to Pseudo-Mercator and reclassify raster into 0/1
fr14 <- projectRaster(fr14, crs="+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs")
fr14 <- reclassify(fr14, c(-Inf, 0.5, 0, 0.5, Inf, 1))

## Loop through all buffer sizes
for(i in 1:length(buffer_distances)){
  
  ## Loop through each point
  for(j in 1:nrow(df)){
    
    ## Clip points by buffer radius
    pbuf <- gBuffer(df[j,], width=buffer_distances[i])
    buf <- mask(fr14, pbuf)
    buf <- crop(buf, pbuf)
    
    ## Extract fragstats
    frg <- PatchStat(buf, cellsize = 30, latlon = FALSE) 
    frg$tile_id <- factor(df[j,]$tile_id)
    frg$dist <- buffer_distances[i]
    
    # Merge to data
    if(exists("frg.df")){
      frg.df <- rbind(frg.df, frg)
    }
    
    if(!exists("frg.df")){
      frg.df <- frg
    }
    
    print(j)
    gc()
    
  }
  
  ## Write csv, remove frg.df
  write.csv(frg.df, paste0("C:/R/EJ_Pknowlesi_2021/outputs/frag14_", buffer_distances[i], ".csv"))    # substitute tile code
  rm(frg.df)
  
  
}
