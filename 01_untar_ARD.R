##########################################################################################
### Untar Landsat ARD files and generate directories for further processing 

library(foreach)
library(iterators)
library(doParallel)

#Register the parallel backend
registerDoParallel(16)

args = commandArgs(trailingOnly=T)
tile_name <- args[1]

#tile_name <- 'h19v06'

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,sep=''))

dir.create(paste(getwd(),'/DEM/',sep=''))
dir.create(paste(getwd(),'/IMG/',sep=''))
dir.create(paste(getwd(),'/MAPS/',sep=''))
dir.create(paste(getwd(),'/NLCD/',sep=''))
dir.create(paste(getwd(),'/PHENO/',sep=''))
dir.create(paste(getwd(),'/SHP/',sep=''))

in_dirs_SR <- list.files(path=getwd(),pattern=glob2rx("*SR*"),
  full.names=T,include.dirs=T)
in_dirs_TA <- list.files(path=getwd(),pattern=glob2rx("*TA*"),
  full.names=T,include.dirs=T)

dir_names <- substr(in_dirs_TA,57,88)
for (i in 1:length(dir_names)){
  dir.create(paste(getwd(),'/IMG/',dir_names[i],sep=''))
}

TA_all <- foreach(i = 1:length(in_dirs_TA), .combine = rbind) %dopar% {
  print(i)

  dir_name_TA <- substr(in_dirs_TA[i],57,88)  
  untar(in_dirs_TA[i],exdir = paste(getwd(),'/IMG/',dir_name_TA,sep=''))
  unlink(in_dirs_TA[i])
  in_dirs_TAB <- list.files(path=paste(getwd(),'/IMG/',dir_name_TA,sep=''),
    pattern=glob2rx("*TAB*"),full.names=T,include.dirs=T)
  unlink(in_dirs_TAB)

  TA <- i
}

SR_all <- foreach(i = 1:length(in_dirs_SR), .combine = rbind) %dopar% {
  print(i)
  
  dir_name_SR <- substr(in_dirs_SR[i],57,88)  
  untar(in_dirs_SR[i],exdir = paste(getwd(),'/IMG/',dir_name_SR,sep=''))
  unlink(in_dirs_SR[i])
  
  SR <- i
}


