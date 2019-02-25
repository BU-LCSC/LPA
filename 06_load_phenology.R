### This script generates maps of annual/long-term average spring and autumn phenology 
### for each Landsat ARD tile

library(rgdal)
library(raster)

args = commandArgs(trailingOnly=T)
tile_name = args[1]

#tile_name <- 'h23v07'

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/IMG',sep=''))
in_dirs <- list.files(path=getwd(),pattern=glob2rx("L*"),
  full.names=T,include.dirs=T)

img <- raster(paste(in_dirs[1],'/evi2.tif',sep=''))

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/NLCD',sep=''))
LC <- raster('nlcd.tif')
LC_vals <- getValues(LC)

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO',sep=''))
tmp_files <- list.files(pattern = "evi2_phenology_topo*", recursive = TRUE,full.names = TRUE)
chunk <- unlist(lapply(tmp_files,
  function(x) na.omit(as.numeric(unlist(strsplit(unlist(x), "[^0-9]+"))))[2]))

which(!(seq(1,200) %in% chunk))

phen <- matrix(NA,ncell(img),4)
spr_annual <- matrix(NA,ncell(img),36)
aut_annual <- matrix(NA,ncell(img),36)

for (j in 1:length(tmp_files)){
  print(j)

  load(tmp_files[j])

  phen[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,1:4]

  spr_annual[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,5:40]
  aut_annual[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,41:76]
}

# phen_sd <- apply(phen_annual,1,sd,na.rm=TRUE)
# phen_sd[LC_vals==0] <- NA
# phen_sd[which(phen[,2]<0.9)] <- NA
# phen_sd <- round(phen_sd*100)
# s <- setValues(img,phen_sd)
# writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
#   tile_name,'/MAPS/spr_sd',sep=""),format='GTiff',overwrite=TRUE)

phen[LC_vals==0] <- NA
phen[phen<=0] <- NA
phen[which(phen[,2]<0.90),3:4] <- NA

phen[,2] <- round(phen[,2]*100)

names <- c('nobs','rsquare','spr','aut')
for (k in 1:4){
  print(k)
  s <- setValues(img,phen[,k])
  name <- paste(names[k],'.tif',sep="")
  writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/MAPS/',name,sep=""),format='GTiff',overwrite=TRUE)
}


for (k in 1:36){
  print(k)
  yr <- 1981+k
  s <- setValues(img,spr_annual[,k])
  name <- paste('spr',yr,'.tif',sep="")
  writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/MAPS/',name,sep=""),format='GTiff',overwrite=TRUE)
}

for (k in 1:36){
  print(k)
  yr <- 1981+k
  s <- setValues(img,aut_annual[,k])
  name <- paste('aut',yr,'.tif',sep="")
  writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/MAPS/',name,sep=""),format='GTiff',overwrite=TRUE)
}

