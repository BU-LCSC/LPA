## This script calculates phenometrics using the Landsat Phenology algorithm 

library(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)
library(zoo)
library(RColorBrewer)
require(data.table)
library(lubridate)

#Register the parallel backend
registerDoParallel(8)

# Calls Landsat Phenology Algorithm
source(file='/usr3/graduate/emelaas/Code/GitHub/LPA/landsat_pheno_pixel.R')

args = commandArgs(trailingOnly=T)
tile_name <- args[1]
chunk = as.numeric(args[2])

#chunk <- 1
#tile_name <- 'h25v15'

#---------------------------------------------------------------------

system.time({
  #Find all images for each site
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx("evi2_topocorr.tif"),full.names=T,include.dirs=T,recursive=TRUE)
  
  r <- stack(in_dirs_tile)
  c <- crop(r,extent(r, 25*(chunk-1)+1, 25*(chunk), 1, 5000))
  vi_vals <- getValues(c)
  vi_vals[vi_vals<=0] <- NA
  
  #Generate index of ancillary info for all images
  doy <- yday(ymd(substring(in_dirs_tile,76,83)))
  yr <- as.numeric(substring(in_dirs_tile,76,79))
  sat <- substring(in_dirs_tile,64,64)
  
  index <- data.frame(sat,yr,doy,t(vi_vals))
  colnames(index) <- c('sat','yr','doy','vi')
  index <- index[order(yr,doy),]
  
  vi_vals <- index[,-c(1,2,3)]
  
  dt.evi <- data.table(vi_vals,keep.rownames=FALSE)
  info <- as.matrix(index[,1:3])
  sat <- as.numeric(info[,1])
  yr <- as.numeric(info[,2])
  doy <- as.numeric(info[,3])
  
  prd <- doy
  prd[yr>=1982 & yr<=1986] <- 1
  prd[yr>=1987 & yr<=1989] <- 2
  prd[yr>=1990 & yr<=1992] <- 3
  prd[yr>=1993 & yr<=1995] <- 4
  prd[yr>=1996 & yr<=1998] <- 5
  prd[yr>=1999 & yr<=2001] <- 6
  prd[yr>=2002 & yr<=2004] <- 7
  prd[yr>=2005 & yr<=2007] <- 8
  prd[yr>=2008 & yr<=2010] <- 9
  prd[yr>=2011 & yr<=2014] <- 10
  prd[yr>=2015 & yr<=2017] <- 11
  
  
  all_pheno <- foreach(i = 1:ncol(vi_vals), .combine = rbind) %dopar% {
    if (i%%10000==0) print(i)
    
    dYR <- 2050
    evi_ts <- as.matrix(dt.evi[,..i]/10000)
    pheno_metrics <- try(Landsat_Phenology(evi_ts,doy,yr,prd,dYR))
  }
  
  all_pheno[all_pheno==0] <- NA
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO',sep=''))
  save(all_pheno,file=paste('evi2_phenology_topocorr_',chunk,sep=""))
})
