##########################################################################################
#LANDSAT PHENOLOGY ALGORITHM RUN ON EACH PIXEL WITH SUFFICIENT NO. OF OBSERVATIONS

#Phenology algorithm parameters:
#evi_ts - EVI time series
#doy - Julian day of year
#year - Year
#prd - three-year period 
#dYR - disturbance year 
Landsat_Phenology <- function(evi_ts,doy,yr,prd,dYR){
  
  #normalized EVI threshold for long-term mean phenology date used to calculate annual phenology
  thresh <- 0.5 
  
  #Initialize Phenology matrices
  phenoSPR <- matrix(0,1,36) #annual spring phenology (1982-2017)
  phenoAUT <- matrix(0,1,36) #annual autumn phenology (1982-2017)
  ltmSPR <- matrix(0,1,1) #long-term mean spring phenology date
  ltmAUT <- matrix(0,1,1) #long-term mean autumn phenology date
  bline <- matrix(0,1,1) #long-term mean winter background EVI
  hline <- matrix(0,1,1) #long-term mean summer maximum EVI
  rsmooth <- matrix(0,1,1) #correlation between observed and smoothed EVI
  nobs <- matrix(0,1,1) #number of available EVI observations used
  
  SPR_thresh = matrix(0,1,1)
  AUT_thresh = matrix(0,1,1)
  SPRsmooth_max = matrix(0,1,1)
  SPRsmooth_min = matrix(0,1,1)
  AUTsmooth_max = matrix(0,1,1)
  AUTsmooth_min = matrix(0,1,1)
  
  info_spr <- NA
  info_aut <- NA
    
  #Set erroneous EVI values to NA
  EVI <- evi_ts
  EVI[EVI==-Inf | EVI==Inf]<-NA
  EVI[EVI<0 | EVI>1]<-NA
  
  #Remove NAs and, if necessary, observations after disturbance year (dYR)
  if (dYR!=0) pos <- which(is.na(EVI)==0 & yr<dYR)
  if (dYR==0) pos <- which(is.na(EVI)==0)
  EVI <- EVI[pos]
  DOY <- doy[pos]
  YR <- yr[pos]
  PRD <- prd[pos]
  nobs <- length(which(is.na(EVI)==0)) #count number of available observations
  
  #Normalize EVI time series using 10% and 90% quantiles for EVI during three-year windows
  DF <- data.frame(PRD,EVI)
  colnames(DF) <- c('yr','evi')
  q1 <- with(DF,tapply(EVI,PRD,quantile,probs=0.10))
  q2 <- with(DF,tapply(EVI,PRD,quantile,probs=0.90))
  
  quant1 <- matrix(NA,11,1)
  quant2 <- matrix(NA,11,1)
  quant1[as.numeric(names(q1))] <- q1
  quant2[as.numeric(names(q2))] <- q2
  
  EVImax <- quant2[PRD]
  EVImin <- quant1[PRD]
  EVInorm <- (EVI-EVImin)/(EVImax-EVImin)
  
  EVInorm[EVInorm==Inf | EVInorm==-Inf]<-NA
  
  pos <- which(is.na(EVInorm)==0 & DOY<=365)
  EVInorm <- EVInorm[pos]
  DOY <- DOY[pos]
  YR <- YR[pos]
  
  if (nobs > 100){
    #Sort DOY/EVI by DOY and compute winter background EVI for improved smoothing spline fit
    x <- cbind(DOY,EVInorm) #Spring background
    x <- x[order(x[,1]),]
    start <- seq(1,x[1,1],1) 
    end <- seq(x[nrow(x),1],365,1)
    YRextra <- c(YR,rep(2012,each=length(start)),rep(2012,each=length(end)))
    DOYextra <- c(DOY,start,end)
    EVIextra <- c(EVInorm,matrix(0,length(start)+length(end),1))
    
    #Fit smoothing spline through EVI and DOY data
    fit <- smooth.spline(DOYextra,EVIextra,spar=0.55)
    EVIsmooth <- data.frame(predict(fit,x=1:365))
    rsmooth <- cor(EVIsmooth[DOY,2],EVInorm)
    
    if (rsmooth>0.85){
      #Separate spline into spring and autumn segments using annual maximum      
      pkval <- which.max(EVIsmooth[,2])
      SPRsmooth <- EVIsmooth[1:pkval,2]
      AUTsmooth <- EVIsmooth[(pkval+1):365,2];
      
      #Compute half-maximum of spring logistic for "ruling in" image dates (points) for
      #anamoly calculation
      SPR_thresh <- which.min(abs(SPRsmooth-thresh))
      SPR_halfmax <- which.min(abs(SPRsmooth-0.5))
      
      #Find anomalies inside of designated box
      SPRsmooth_max <- 1
      SPRsmooth_min <- 0
      box_max <- SPRsmooth_max-0.2*(SPRsmooth_max-SPRsmooth_min)
      box_min <- SPRsmooth_min+0.2*(SPRsmooth_max-SPRsmooth_min)
      
      #Generate a matrix with candidate spring phenology observations
      if (is.na(box_max)==0 && is.na(box_min)==0){
        info <- cbind(YR,DOY,EVInorm,matrix(0,length(DOY),1))
        info <- info[order(info[,1],info[,2]),]
        info <- rbind(matrix(0,1,4),info)
        pos <- which(is.na(info[,3])==0 & info[,3]>box_min & info[,3]<box_max &
            info[,2]<SPR_halfmax+20 & info[,2]>SPR_halfmax-20)
        if (length(pos)>4){
          k <- 1
          for (p in 1:length(pos)){
            smooth_ratio <- abs(SPRsmooth-info[pos[p],3])
            info[pos[p],4] <- which.min(smooth_ratio)
            
            if (k==1)
              info_spr <- info[pos[p],]
            else
              info_spr <- rbind(info_spr,info[pos[p],])
            
            k <- k+1
          }
          
          #Remove repeat observations
          w <- which((info[pos,4]<info[(pos-1),4] & info[pos,1]==info[(pos-1),1])==1)
          if (length(w)>0) info_spr <- info_spr[-w,]
        }
      }
      
      #Compute half-maximum of spring logistic for "ruling in" image dates (points) for
      #anamoly calculation
      AUT_thresh <- which.min(abs(AUTsmooth-thresh))+pkval
      AUT_halfmax <- which.min(abs(AUTsmooth-0.5))+pkval
      
      #Find anomalies inside of designated box
      AUTsmooth_max <- 1
      AUTsmooth_min <- 0
      box_max <- AUTsmooth_max-0.3*(AUTsmooth_max-AUTsmooth_min)
      box_min <- AUTsmooth_min+0.2*(AUTsmooth_max-AUTsmooth_min)
      
      #Generate matrix with candidate autumn phenology observations
      if (is.na(box_max)==0 && is.na(box_min)==0){
        info <- cbind(YR,DOY,EVInorm,matrix(0,length(DOY),1))
        info <- info[order(info[,1],info[,2]),]
        info <- rbind(info,matrix(0,1,4))
        pos <- which(is.na(info[,3])==0 & info[,3]>box_min & info[,3]<box_max &
            info[,2]<AUT_halfmax+20 & info[,2]>AUT_halfmax-20)
        if (length(pos)>4){
          k <- 1
          for (p in 1:length(pos)){
            smooth_ratio <- abs(AUTsmooth-info[pos[p],3])
            info[pos[p],4] <- which.min(smooth_ratio)+pkval
            
            if (k==1)
              info_aut <- info[pos[p],]
            else
              info_aut <- rbind(info_aut,info[pos[p],])
            
            k <- k+1
          }
          w <- which((info[pos,4]<info[(pos+1),4] & info[pos,1]==info[(pos+1),1])==1)
          if (length(w)>0) info_aut <- info_aut[-w,]
        }
      }
      
      #Calculate interannual phenology dates by taking the distance
      #between each candidate observation and where the same magnitude 
      #of EVI occurs on the spline
      if (exists('info_spr') == 1 && exists('info_aut') == 1){
        if (ncell(info_spr) > 20 & ncell(info_aut) > 20) {
          info_spr <- cbind(info_spr,SPR_thresh+(info_spr[,2]-info_spr[,4]))
          info_aut <- cbind(info_aut,AUT_thresh+(info_aut[,2]-info_aut[,4]))
          
          for (y in 1982:2017){
            pos1 <- which(info_spr[,1] == y)
            pos2 <- which(info_aut[,1] == y)
            
            if (length(pos1) > 0 && ncol(info_spr)==5){
              w <- which.min(abs(info_spr[pos1,4]-0.5))
              phenoSPR[1,y-1981] <- ceiling(mean(info_spr[pos1[w],5]))
            }
            if (length(pos2) > 0 && ncol(info_aut)==5){
              w <- which.min(abs(info_aut[pos2,4]-0.5))
              phenoAUT[1,y-1981] <- ceiling(mean(info_aut[pos2[w],5]))
            }
          }
        }
        
        ltmSPR <- SPR_thresh
        ltmAUT <- AUT_thresh
        
        remove(info_spr,info_aut) 
      }
    }
  }
  
  pheno_matrix <- cbind(nobs,rsmooth,ltmSPR,ltmAUT,phenoSPR,phenoAUT)
  
  if (is.character(pheno_matrix) == 1){
    pheno_matrix <- matrix(NA,1,76)
  } 
  
  return(pheno_matrix)
}

