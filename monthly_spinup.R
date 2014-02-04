#This code processes and plots annual output from the ED2 model spinup period
#Jaclyn Hatala Matthes, 1/30/14
#jaclyn.hatala.matthes@gmail.com

#Load libraries
ok = require(chron,lib.loc="/usr4/spclpgm/jmatthes/"); if (! ok) stop("Package chron is not available...")
ok = require(ncdf,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package ncdf is not available...")
ok = require(maps,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package maps is not available...")
ok = require(sp,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")
ok = require(raster,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")
ok = require(colorspace,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")

#Set sites
sites <- c("PBL","PHA","PMB","PDL","PHO","PUN")

#Some plotting stuff
plot.path  <- paste("/projectnb/cheas/paleon/ED_runs/phase1a_runs/",sites[s],"/plots/",sep="")
plot.place <- paste(sites[s],": Spin-up",sep="")
plot.ptsz   = 14     # Default font size
plot.width  = 9.7    # Window width
plot.height = 6.0    # Window.heights

#Names and colors for all PFTs, so it works for tropical and temperate.
pft.names = c("C4 grass"          ,"Early tropical"    ,"Mid tropical"      
              ,"Late tropical"     ,"Temperate C3 Grass","North Pine"        
              ,"South Pine"        ,"Late conifer"      ,"Early hardwood"    
              ,"Mid hardwood"      ,"Late hardwood"     ,"C3 crop"           
              ,"C3 pasture"        ,"C4 crop"           ,"C4 pasture"        
              ,"C3 grass"          ,"Araucaria"         ,"Total"             )
pft.cols  = c("gold"              ,"chartreuse"        ,"chartreuse4"       
              ,"#004E00"           ,"mediumpurple1"     ,"deepskyblue"       
              ,"mediumturquoise"   ,"royalblue4"        , "darkorange"       
              ,"orangered"         ,"firebrick4"         , "purple4"          
              ,"darkorchid1"       ,"darkgoldenrod"     ,   "khaki"          
              ,"lightgoldenrod3"   ,"steelblue3"        ,   "grey22"         )
n.pft     = length(pft.names) - 1

pft <- c(5,6,8,9,10,11)
suffix <- "g01.h5"

for(s in 1:length(sites)){
  #Set directories
  dat.dir    <- paste("/projectnb/cheas/paleon/ED_runs/phase1a_runs/",sites[s],"/analy2/",sep="")
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
      
      #calculate monthly index for storage
      if(y==yeara){
        mon.ind <- m-montha+1
      }else{
        mon.ind <- ((y-1)-yeara)*12+m+montha
      }
      print(mon.ind)
      
      #Grab global & patch level met variables.
      npoly     <- get.var.ncdf(now,'NPOLYGONS_GLOBAL')
      nsites    <- get.var.ncdf(now,'NSITES_GLOBAL')
      npatches  <- get.var.ncdf(now,'NPATCHES_GLOBAL')
      ncohorts  <- get.var.ncdf(now,'NCOHORTS_GLOBAL')
      atm_pres[mon.ind]      <- get.var.ncdf(now,'MMEAN_ATM_PRSS')
      atm_vpd[mon.ind]       <- get.var.ncdf(now,'MMEAN_ATM_VPDEF')
      atm_pre[mon.ind]       <- get.var.ncdf(now,'MMEAN_PCPG')
      atm_tmp[mon.ind]       <- get.var.ncdf(now,'MMEAN_ATM_TEMP')-273.15
      atm_rnt[mon.ind]       <- get.var.ncdf(now,'MMEAN_RNET')
      atm_trsp[mon.ind]      <- get.var.ncdf(now,'MMEAN_TRANSP')
      atm_swat[mon.ind]      <- get.var.ncdf(now,'MMEAN_SOIL_WATER')
      npp_site[mon.ind]      <- get.var.ncdf(now,"MMEAN_NPPDAILY")
      
      #Grab cohort level variables.21600/
      ipft      <- get.var.ncdf(now,'PFT')
      dbh       <- get.var.ncdf(now,'DBH')
      nplant    <- get.var.ncdf(now,'NPLANT')
      lai       <- get.var.ncdf(now,'LAI_CO')
      npp       <- get.var.ncdf(now,'MMEAN_NPPDAILY_CO')
      npp.croot <- get.var.ncdf(now,'MMEAN_NPPCROOT_CO')
      npp.froot <- get.var.ncdf(now,'MMEAN_NPPFROOT_CO')
      npp.leaf  <- get.var.ncdf(now,'MMEAN_NPPLEAF_CO')
      npp.leaf  <- get.var.ncdf(now,'MMEAN_NPPWOOD_CO')
      gpp       <- get.var.ncdf(now,'MMEAN_GPP_CO')
      
      rleaf <- get.var.ncdf(now,'MMEAN_LEAF_RESP_CO')
      rroot <- get.var.ncdf(now,'MMEAN_ROOT_RESP_CO')
      rgrow <- get.var.ncdf(now,'MMEAN_GROWTH_RESP_CO')
      rstor <- get.var.ncdf(now,'MMEAN_STORAGE_RESP_CO')
      rvleaf <- get.var.ncdf(now,'MMEAN_VLEAF_RESP_CO')
      mainl <- get.var.ncdf(now,'MMEAN_LEAF_MAINTENANCE')
      mainr <- get.var.ncdf(now,'MMEAN_ROOT_MAINTENANCE')
      leafdrop <- get.var.ncdf(now,'MMEAN_LEAF_DROP_CO')
      agb <- get.var.ncdf(now,'AGB_CO')
      cbal <- ((gpp - rleaf-rroot-rgrow-rstor-rvleaf-mainl-mainr-leafdrop)/agb)
      cbr_bar <- get.var.ncdf(now,'MMEAN_CB')
      PAR       <- get.var.ncdf(now,'MMEAN_PAR_L')
      
      #if any PFTs go extinct, make placeholders
      if(length(unique(ipft))<length(pft)){
        tmp  <- (length(pft)-length(unique(ipft)))
        ipft <- c(ipft,pft[!(pft %in% ipft)])
        cbal <- c(cbal,rep(0,tmp))
        lai  <- c(lai,rep(0,tmp))
        cbr_bar  <- c(cbr_bar,rep(0,tmp))
        PAR <- c(PAR,rep(0,tmp))
      }
      
      cbal.pft[mon.ind,] <- tapply(cbal,ipft,mean)
      cbalbr.pft[mon.ind,] <- tapply(cbr_bar,ipft,mean)
      lai.pft[mon.ind,] <- tapply(lai,ipft,sum)
      PAR.pft[mon.ind,] <- tapply(PAR,ipft,mean)
      
      close.ncdf(now)
    }
  }  
  
  dates <- seq(as.Date(paste(yeara,montha,"1",sep="/")), as.Date(paste(yearz-1,12,"1",sep="/")), "months")

#  png(paste(plot.path,sites[s],'_AGBbyPFT','.png',sep=''),width=900,height=600)
#  pdf(paste(plot.path,sites[s],"_spinup",sep=''))
  
  #Monthly met drivers
  par(mfrow=c(2,2))
  plot(dates,atm_pres[1:length(atm_pres)-1],xlab="spin-up date",main=paste(sites[s],": spin-up drivers"),
       ylab="Pressure [Pa]")
  plot(dates,atm_vpd[1:length(atm_pres)-1],xlab="spin-up date",main=paste(sites[s],": spin-up drivers"),
       ylab="Pressure [Pa]")
  plot(dates,atm_tmp[1:length(atm_pres)-1],xlab="spin-up date",main=paste(sites[s],": spin-up drivers"),
       ylab="Temperature [C]")
  plot(dates,atm_swat[1:length(atm_pres)-1],xlab="spin-up date",main=paste(sites[s],": spin-up drivers"),
       ylab="Soil water [m/m]")
  
  #Monthly LAI by PFT
  par(mfrow=c(1,1))
  plot(dates,lai.pft[,1],col=pft.cols[5],pch=16,ylim=range(lai.pft),
       xlab="spin-up date",ylab="Monthly LAI",
       main=paste(sites[s],": Spin-up",sep=""))
  for(p in 2:ncol(lai.pft)){
    points(dates,lai.pft[,p],col=pft.cols[4+p],pch=16)
  }
  legend(dates[2],max(lai.pft)-mean(lai.pft),pft.names[sort(unique(ipft))],col=pft.cols[5:10],pch=16)

  #Monthly NPP
  plot(dates,npp_site[1:length(npp_site)-1],col=pft.cols[5],pch=16,ylim=range(lai.pft),
       xlab="spin-up date",ylab="Mean daily NPP [kg m-2]",
       main=paste(sites[s],": Spin-up",sep=""))
  
  #Monthly C balance by PFT
  plot(dates,cbal.pft[,1],col=pft.cols[5],pch=16,ylim=c(-5,5),
       xlab="spin-up date",ylab="Monthly C Balance",
       main=paste(sites[s],": Spin-up",sep=""))
  for(p in 2:ncol(cbal.pft)){
    points(dates,cbal.pft[,p],col=pft.cols[4+p],pch=16)
  }
  legend(dates[2],max(lai.pft)-mean(lai.pft),pft.names[sort(unique(ipft))],col=pft.cols[5:10],pch=16)
  
  #Monthly C balance ratio by PFT
  plot(dates,cbalbr.pft[,1],col=pft.cols[5],pch=16,ylim=c(0,0.02),
       xlab="spin-up date",ylab="Monthly C Balance",
       main=paste(sites[s],": Spin-up",sep=""))
  for(p in 2:ncol(cbalbr.pft)){
    points(dates,cbalbr.pft[,p],col=pft.cols[4+p],pch=16)
  }
  legend(dates[2],0.8,pft.names[sort(unique(ipft))],col=pft.cols[5:10],pch=16)
  
#  dev.off()
  rm(dat.dir,yeara,yearz,files,ann.files,match.files)
}

