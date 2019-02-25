## This script generates EVI2 GeoTIFF for each image in the ARD stack
## Still working on topographic correction component (need to modify for each LC)

require(raster)
require(rgdal)
require(gdalUtils)
require(rgeos)
require(jsonlite)
require(spatialEco)
library(foreach)
library(iterators)
library(doParallel)
library(landsat)

#Register the parallel backend
registerDoParallel(16)

args = commandArgs(trailingOnly=T) 
tile_name = args[1]  

#tile_name <- 'h17v02'

H <- as.numeric(substring(tile_name,2,3))
V <- as.numeric(substring(tile_name,5,6))

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))

## GENERATE ARD TILE AND SAVE TO DIRECTORY

ard_tiles <- readOGR('/projectnb/modislc/projects/landsat_sentinel/ARD/CONUS_ARD_grid/',
  'conus_ard_grid')
tile <- ard_tiles[which(ard_tiles$h==H & ard_tiles$v==V),]
# Determine lat/lon extent of tile for NED download
tile_latlon <- spTransform(tile,CRS("+proj=longlat +datum=WGS84"))
extent(tile_latlon)
tile_proj <- spTransform(tile,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
writeOGR(tile_proj,paste(getwd(),'/','SHP',sep=''),
  tile_name,driver="ESRI Shapefile",overwrite=TRUE)


## GENERATE DEM, SLOPE, ASPECT and LC MAPS

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/SHP',sep=''))
src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_dem_aea_proj.tif'
DEM <- gdalwarp(src_data,
  dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM/DEM.tif',sep=''),
  t_srs=projection(tile_proj),ts=c(5000,5000),
  cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
  crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)

src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_slope_aea_proj.tif'
slope <- gdalwarp(src_data,
  dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM/slope.tif',sep=''),
  t_srs=projection(tile_proj),ts=c(5000,5000),
  cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
  crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)
slope_vals <- getValues(slope)

src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_aspect_aea_proj.tif'
aspect <- gdalwarp(src_data,
  dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM/aspect.tif',sep=''),
  t_srs=projection(tile_proj),ts=c(5000,5000),
  cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
  crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)
aspect_vals <- getValues(aspect)

src_data2 <- '/projectnb/modislc/data/lc_database/regional/united_states/NLCD2006_landcover_4-20-11_se5/nlcd2006_landcover_4-20-11_se5.img'
LC <- gdalwarp(src_data2,dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/NLCD/nlcd.tif',sep=''),
  t_srs=projection(tile_proj),ts=c(5000,5000),
  cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
  crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)
LC_vals <- getValues(LC)
LC_types <- as.numeric(rownames(table(LC_vals)))


## LOAD IN SURFACE REFLECTANCE BANDS, QA LAYER, SOLAR ZENITH & AZIMUTH

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/IMG',sep=''))
in_dirs <- list.files(path=getwd(),pattern=glob2rx("L*"),
  full.names=T,include.dirs=T)

# Check for images with missing evi outputs
setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx("evi2_topocorr.tif"),full.names=T,include.dirs=T,recursive=TRUE)
ydoy <- substring(in_dirs_tile,76,83)
ydoy2 <- substring(in_dirs,76,83)
w <- which(!(ydoy2 %in% ydoy))

all_count <- foreach(i = w, .combine = rbind) %dopar% {

  # Landsat 8
  if (substr(in_dirs[i],64,64)==8){
    nir <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB5.tif',sep='')) # near infrared reflectance
    red <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB4.tif',sep='')) # red reflectance
    QA <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_PIXELQA.tif',sep='')) # pixel QA

    s <- stack(red,nir,QA)
    band_vals <- getValues(s)
    band_vals[band_vals<=0] <- NA

    qa <- band_vals[,3]
    w <- which(qa!=322 & qa!=386 & qa!=834 & qa!=898 & qa!=1346)
    band_vals[w,1:2] <- NA

    # Landsat 4-7
  } else {
    nir <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB4.tif',sep='')) # near infrared reflectance
    red <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB3.tif',sep='')) # red reflectance
    QA <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_PIXELQA.tif',sep='')) # pixel QA

    s <- stack(red,nir,QA)
    band_vals <- getValues(s)
    band_vals[band_vals<=0] <- NA

    qa <- band_vals[,3]
    w <- which(qa!=66 & qa!=130)
    band_vals[w,1:2] <- NA
  }
  
  evi2 <- 2.5*(band_vals[,2]/10000 - band_vals[,1]/10000)/(band_vals[,2]/10000 + 2.4*band_vals[,1]/10000 + 1)
  evi2_map <- setValues(nir,round(evi2*10000))
  writeRaster(evi2_map,filename=paste(in_dirs[i],'/evi2',sep=""),
    format='GTiff',overwrite=TRUE)

  ndvi <- (band_vals[,2]/10000 - band_vals[,1]/10000)/(band_vals[,2]/10000 + band_vals[,1]/10000)
  ndvi_map <- setValues(nir,round(ndvi*10000))
  writeRaster(ndvi_map,filename=paste(in_dirs[i],'/ndvi',sep=""),
    format='GTiff',overwrite=TRUE)
  
  # TOPOGRAPHIC CORRECTION 
  SAA <- 0.01*raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
    '_C01_V01_SOA4.tif',sep='')) # solar azimuth angle
  SZA <- 0.01*raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
    '_C01_V01_SOZ4.tif',sep='')) # solar zenith angle
  
  SAA_vals <- getValues(SAA)
  SZA_vals <- getValues(SZA)

  b3_corr <- as.numeric(matrix(NA,ncell(SZA),1))
  b4_corr <- as.numeric(matrix(NA,ncell(SZA),1))

  for (n in 1:length(LC_types)){
    print(c(i,n))
    wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,1])==0)
    if (length(wLC) > 4){
      b3_corr_LC <- topocorr_v3(x=band_vals[wLC,1],slope=slope_vals[wLC], 
        aspect=aspect_vals[wLC], sunelev=90-SZA_vals[wLC], sunazimuth=SAA_vals[wLC],
        method='rotational')
      b3_corr[wLC] <- b3_corr_LC
    }
    
    wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,2])==0)
    if (length(wLC) > 4){
      b4_corr_LC <- topocorr_v3(x=band_vals[wLC,2],slope=slope_vals[wLC], 
        aspect=aspect_vals[wLC], sunelev=90-SZA_vals[wLC], sunazimuth=SAA_vals[wLC],
        method='rotational')
      b4_corr[wLC] <- b4_corr_LC
    }
  }
  
  evi2 <- 2.5*(b4_corr/10000 - b3_corr/10000)/(b4_corr/10000 + 2.4*b3_corr/10000 + 1)
  evi2_map <- setValues(nir,round(evi2*10000))
  writeRaster(evi2_map,filename=paste(in_dirs[i],'/evi2_topocorr',sep=""),
    format='GTiff',overwrite=TRUE)
  
  ndvi <- (b4_corr/10000 - b3_corr/10000)/(b4_corr/10000 + b3_corr/10000)
  ndvi_map <- setValues(nir,round(ndvi*10000))
  writeRaster(ndvi_map,filename=paste(in_dirs[i],'/ndvi_topocorr',sep=""),
    format='GTiff',overwrite=TRUE)

  count <- i
}


