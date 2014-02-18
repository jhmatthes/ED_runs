#This script aggregates the output from the initial ED2 spin-up 
#to calculate the semi-analytical solution

#Load libraries
ok = require(chron,lib.loc="/usr4/spclpgm/jmatthes/"); if (! ok) stop("Package chron is not available...")
ok = require(ncdf,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package ncdf is not available...")
ok = require(maps,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package maps is not available...")
ok = require(sp,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")
ok = require(raster,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")
ok = require(colorspace,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")

#Set sites
sites <- c("PBL","PHA","PMB","PDL","PHO","PUN")

pft <- c(5,6,8,9,10,11)
suffix <- "g01.h5"
dpy <- c(31,28,31,30,31,30,31,31,30,31,30,31)

for(s in 1:length(sites)){
  #Set directories
  dat.dir    <- paste("/projectnb/cheas/paleon/ED_runs/phase1a_runs/",sites[s],"/analy4/",sep="")
  match.files <- grep("-E-",list.files(dat.dir))
  files <- list.files(dat.dir)
  mon.files  <- files[match.files] #monthly files only
  
  #Get time window
  yeara <- as.numeric(strsplit(mon.files,"-")[[1]][3]) #first year
  yearz <- as.numeric(strsplit(mon.files,"-")[[length(mon.files)]][3]) #last year
  montha <- as.numeric(strsplit(mon.files,"-")[[1]][4]) #first month
  monthz <- as.numeric(strsplit(mon.files,"-")[[length(mon.files)]][4]) #first month
  
  #Set up storage
  atm_pres <- atm_vpd <- atm_pre <- atm_tmp <- atm_rnt <- atm_trsp <- atm_swat <- npp_site <- vector()
  cbal.pft <- cbalbr.pft <- PAR.pft <- lai.pft <- matrix(nrow=((yearz-yeara-1)*12+(12-montha+1)),
                                                         ncol=length(pft))
  fnpp <- matrix(ncol=8,nrow=(yearz-yeara+1))
  #monthly storage
  npp <- leaf.drop <- npp.leaf <- npp.croot <- npp.froot <- npp.sapwood <- npp.wood <- npp.seed <- matrix(ncol=12,nrow=(yearz-yeara+1))
  
  #loop over months and aggregate monthly mean data
  for(y in yeara:yearz){
    
    #calculate month start/end based on year 
    if (y == yeara){
      month.begin = montha
    }else{
      month.begin = 1
    }
    if (y == yearz){
      month.end = monthz
    }else{
      month.end = 12
    }
    
    for(m in month.begin:month.end){   
  
      #Make the file name. 
      year.now  = sprintf("%4.4i",y)
      month.now = sprintf("%2.2i",m)
      day.now   = sprintf("%2.2i",0)
      hour.now  = sprintf("%6.6i",0)
      file.now  = paste(sites[s],"spin","-E-",year.now,"-",month.now,"-",day.now,"-"
                        ,hour.now,"-",suffix,sep="")
      
      cat(" - Reading file :",file.now,"...","\n")
      now <- open.ncdf(paste(dat.dir,file.now,sep=""))
      
      #Calculate monthly NPP pools
      yind                <- y-yeara+1
      npp[yind,m]         <- get.var.ncdf(now,"MMEAN_NPPDAILY")/365*dpy[m]
      npp.leaf[yind,m]    <- get.var.ncdf(now,"MMEAN_NPPLEAF")/365*dpy[m]
      npp.croot[yind,m]   <- get.var.ncdf(now,"MMEAN_NPPCROOT")/365*dpy[m]
      npp.froot[yind,m]   <- get.var.ncdf(now,"MMEAN_NPPFROOT")/365*dpy[m]
      npp.sapwood[yind,m] <- get.var.ncdf(now,"MMEAN_NPPSAPWOOD")/365*dpy[m]
      npp.wood[yind,m]    <- get.var.ncdf(now,"MMEAN_NPPWOOD")/365*dpy[m]
      npp.seed[yind,m]    <- get.var.ncdf(now,"MMEAN_NPPSEEDS")/365*dpy[m]
      leaf.drop[yind,m]   <- mean(get.var.ncdf(now,"MMEAN_LEAF_DROP_CO"))

      close.ncdf(now)
    }
    
    #Find fraction of total annual NPP transferred to pools
    yind <- y-yeara+1
    fnpp[yind,1]        <- sum(npp.leaf[yind,],npp.croot[yind,],npp.froot[yind,],npp.sapwood[yind,],
                               npp.wood[yind,],npp.seed[yind,])/sum(npp[yind,]) 
    fnpp[yind,2]        <- sum(npp.leaf[yind,])/sum(npp[yind,]) 
    fnpp[yind,3]        <- sum(npp.croot[yind,])/sum(npp[yind,]) 
    fnpp[yind,4]        <- sum(npp.froot[yind,])/sum(npp[yind,]) 
    fnpp[yind,5]        <- sum(npp.sapwood[yind,])/sum(npp[yind,]) 
    fnpp[yind,6]        <- sum(npp.wood[yind,])/sum(npp[yind,]) 
    fnpp[yind,7]        <- sum(npp.seed[yind,])/sum(npp[yind,]) 
    
    #Get annual ecosystem C pools [kgC/m2]
    file.now  = paste(sites[s],"spin","-Y-",year.now,"-00-",day.now,"-"
                      ,hour.now,"-",suffix,sep="")
    
    cat(" - Reading file :",file.now,"...","\n")
    now <- open.ncdf(paste(dat.dir,file.now,sep=""))
    
    balive[y-yeara+1]    <- sum(get.var.ncdf(now,'BALIVE')) 
    bdead[y-yeara+1]     <- sum(get.var.ncdf(now,'BDEAD'))
    bleaf[y-yeara+1]     <- sum(get.var.ncdf(now,'BLEAF')) 
    broot[y-yeara+1]     <- sum(get.var.ncdf(now,'BROOT')) 
    bsapa[y-yeara+1]     <- sum(get.var.ncdf(now,'BSAPWOODA')) 
    bsapb[y-yeara+1]     <- sum(get.var.ncdf(now,'BSAPWOODB')) 
    sfast[y-yeara+1]     <- get.var.ncdf(now,'FAST_SOIL_C') 
    sslow[y-yeara+1]     <- get.var.ncdf(now,'SLOW_SOIL_C') 
    sstruc[y-yeara+1]    <- get.var.ncdf(now,'STRUCTURAL_SOIL_C')
    
  }
}
