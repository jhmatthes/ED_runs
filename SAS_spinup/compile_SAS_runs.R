#This file compiles the ensemble output from the SAS runs 
#Jaclyn Hatala Matthes, 2/18/14
#jaclyn.hatala.matthes@gmail.com
#sites <- c("PBL","PHA","PMB","PHO","PUN","PDL")
#site.coord <- c("lat46.5lon-94.5","lat42.5lon-72.5","lat43.5lon-82.5",
#                "lat45.5lon-68.5","lat46.5lon-89.5","lat47.5lon-95.5")
sites <- "PDL"
site.coord <- "lat47.5lon-95.5"

#Load libraries
ok = require(chron,lib.loc="/usr4/spclpgm/jmatthes/"); if (! ok) stop("Package chron is not available...")
ok = require(ncdf,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package ncdf is not available...")
ok = require(colorspace,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")

#Setup analysis file structure
base  <- "/projectnb/cheas/paleon/ED_runs/p1a_spin_042214/"
out   <- "/projectnb/cheas/paleon/ED_runs/SAS_spinup/phase1a_spinup/"
#sites <- c("PBL","PHA","PMB","PDL","PHO","PUN")
blckyr<- 1 #number of years to chunk data by
niter <- length(list.dirs(paste(base,sites[2],"/",sep=""),recursive=FALSE)) #iterations/site 
pft   <- c(5,6,8,9,10,11) #set of PFTs used in analysis
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)
sufx  <- "g01.h5"

#constants from ED2 for SAS solution for soil pools
fsc_loss <- 1.0/11.0
ssc_loss <- 1.0/100.2
ssl_loss <- 1.0/5.4
resp_temperature_increase <- 0.23503
rel_soil_moist <- 0.5
resp_opt_water            <- 0.8938
resp_water_below_opt      <- 5.0786
Lc                        <- 0.049787
c2n_slow                  <- 10.0
c2n_structural            <- 150.0
r_stsc                    <- 0.3
soil_tempk <- c(278.5,279.3, 280.0, 277.3, 276.8, 277.2)+10 #mean temp per site
temperature_limitation = resp_temperature_increase * exp(308.56 * (1./56.02 - 1./(soil_tempk-227.15)))
#water_limitation <- exp((rel_soil_mois - resp_opt_water) * resp_water_below_opt)
water_limitation <- rel_soil_moist*4.0893 + rel_soil_moist^2*-3.1681 - 0.3195897
A_decomp <- temperature_limitation * water_limitation 

#First loop over analy files (faster than histo) to aggregate initial 
#.css and .pss files for each site
for(s in sites){
  
  #Set directories
  dat.dir    <- paste(base,s,"/analy/",sep="")
  match.files <- grep("-Y-",list.files(dat.dir))
  files <- list.files(dat.dir)
  ann.files  <- files[match.files] #yearly files only
  
  #Get time window
  yeara <- as.numeric(strsplit(ann.files,"-")[[1]][3]) #first year
#  yearz <- as.numeric(strsplit(ann.files,"-")[[length(ann.files)]][3]) #last year
  yearz <- 1999

  #storage
  pss.big <- matrix(nrow=floor((yearz-yeara+1)/blckyr)+1,ncol=14) #save every 10 yrs
  colnames(pss.big) <- c("site","year","patch","dst","age","area","water","fsc","stsc","stsl",
                         "ssc","psc","msn","fsn")
  
  for (y in yeara:yearz){
    if(y%%blckyr == 0){
      cat(" - Reading file :",ann.files[y-yeara+1],"...","\n")
      now <- open.ncdf(paste(dat.dir,ann.files[y-yeara+1],sep=""))
      
      #Grab variable to see how many cohorts there are
      ipft      <- get.var.ncdf(now,'PFT')
      
      #organize into .css variables
      css.tmp <- matrix(nrow=length(ipft),ncol=10)
      css.tmp[,1] <- rep(1850,length(ipft))
      css.tmp[,2] <- rep(floor((y-yeara)/blckyr)+1,length(ipft))
      css.tmp[,3] <- 1:length(ipft)
      css.tmp[,4] <- get.var.ncdf(now,'DBH')
      css.tmp[,5] <- get.var.ncdf(now,'HITE')
      css.tmp[,6] <- ipft
      css.tmp[,7] <- get.var.ncdf(now,'NPLANT')
      css.tmp[,8] <- get.var.ncdf(now,'BDEAD')
      css.tmp[,9] <- get.var.ncdf(now,'BALIVE')
      css.tmp[,10] <- rep(-999,length(ipft))
      colnames(css.tmp) <- c("year","patch","cohort","dbh","ht","pft","n","bdead","balive","Avgrg")
      
      #save big .css matrix
      if(y==yeara){
        css.big <- css.tmp
      } else{
        css.big <- rbind(css.big,css.tmp)
      }
      
      #save .pss variables
      ind <- (y-yeara)/blckyr + 1
      pss.big[ind,1]  <- 1
      pss.big[ind,2]  <- 1850
      pss.big[ind,3]  <- floor((y-yeara)/blckyr)+1
      pss.big[ind,4]  <- 1
      pss.big[ind,5]  <- y-yeara+1
      pss.big[ind,6]  <- get.var.ncdf(now,"AREA")
      pss.big[ind,7]  <- 0.1
      pss.big[ind,8]  <- get.var.ncdf(now,"FAST_SOIL_C")
      pss.big[ind,9]  <- get.var.ncdf(now,"STRUCTURAL_SOIL_C")
      pss.big[ind,10] <- get.var.ncdf(now,"STRUCTURAL_SOIL_L")
      pss.big[ind,11] <- get.var.ncdf(now,"SLOW_SOIL_C")
      pss.big[ind,12] <- 0
      pss.big[ind,13] <- get.var.ncdf(now,"AVG_MSN")
      pss.big[ind,14] <- get.var.ncdf(now,"AVG_FSN")
      
      close.ncdf(now)
    }
  }
#}


#Second loop over histo files (much slower than analy) to aggregate soil inputs
#for steady-state solution

#storage
fsc_in_y <- ssc_in_y <- ssl_in_y <- fsn_in_y <- pln_up_y <- vector()
fsc_in_m <- ssc_in_m <- ssl_in_m <- fsn_in_m <- pln_up_m <-  vector()

#for(s in sites){
  dat.dir    <- paste(base,s,"/histo/",sep="")
  match.files <- grep("-S-",list.files(dat.dir))
  files <- list.files(dat.dir)
  mon.files  <- files[match.files] #monthly files only
  
  #Get time window
  yeara <- as.numeric(strsplit(mon.files,"-")[[1]][3]) #first year
#  yearz <- as.numeric(strsplit(mon.files,"-")[[length(mon.files)]][3]) #last year
  yearz  <- 1999
  monthz <- 12
  montha <- as.numeric(strsplit(mon.files,"-")[[1]][4]) #first month
#  monthz <- as.numeric(strsplit(mon.files,"-")[[length(mon.files)]][4]) #first month
  
  for (y in yeara:yearz){
    if(y%%blckyr == 0){
      
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
        year.now  <-sprintf("%4.4i",y)
        month.now <- sprintf("%2.2i",m)
        day.now   <- sprintf("%2.2i",1)
        hour.now  <- sprintf("%6.6i",0)
        
        dat.dir     <- paste(base,s,"/histo/",sep="")
        file.now    <- paste(s,"spin","-S-",year.now,"-",month.now,"-",day.now,"-"
                             ,hour.now,"-",sufx,sep="")
        
        cat(" - Reading file :",file.now,"...","\n")
        now <- open.ncdf(paste(dat.dir,file.now,sep=""))
        
        fsc_in_m[m-month.begin+1] <- get.var.ncdf(now,"FSC_IN")*dpm[m] #kg/(m2*day) --> kg/(m2*month)
        ssc_in_m[m-month.begin+1] <- get.var.ncdf(now,"SSC_IN")*dpm[m]
        ssl_in_m[m-month.begin+1] <- get.var.ncdf(now,"SSL_IN")*dpm[m]
        fsn_in_m[m-month.begin+1] <- get.var.ncdf(now,"FSN_IN")*dpm[m]
        pln_up_m[m-month.begin+1] <- get.var.ncdf(now,"TOTAL_PLANT_NITROGEN_UPTAKE")*dpm[m]
        close.ncdf(now)
      }
      ind <- (y-yeara)/blckyr + 1
      fsc_in_y[ind] <- sum(fsc_in_m,na.rm=TRUE)
      ssc_in_y[ind] <- sum(ssc_in_m,na.rm=TRUE)
      ssl_in_y[ind] <- sum(ssl_in_m,na.rm=TRUE)
      fsn_in_y[ind] <- sum(fsn_in_m,na.rm=TRUE)
      pln_up_y[ind] <- sum(pln_up_m,na.rm=TRUE)
    }
  }
  
  #Calculate steady-state soil pools
#  fsc_ss <- median(fsc_in_y)/(fsc_loss * A_decomp)
#  ssc_ss <- median(ssc_in_y)/(ssc_loss * A_decomp)
#  ssl_ss <- median(ssl_in_y)/(ssl_loss * A_decomp * Lc)
#  fsn_ss <- median(fsn_in_y)/(fsc_loss * A_decomp)

  fsc_ss <- fsc_in_y[11]/(fsc_loss * A_decomp)
  ssc_ss <- ssc_in_y[11]/(ssc_loss * A_decomp)
  ssl_ss <- ssl_in_y[11]/(ssl_loss * A_decomp * Lc)
  fsn_ss <- fsn_in_y[11]/(fsc_loss * A_decomp)
  
  #fast_N_loss + slow_C_loss
  msn_med  <- fsc_loss*A_decomp*fsn_in_y[11]+ (ssc_loss * A_decomp)/c2n_slow 
  
  #ED2: csite%mineralized_N_loss  = csite%total_plant_nitrogen_uptake(ipa)             
  # + csite%today_Af_decomp(ipa) * Lc * K1 * csite%structural_soil_C(ipa)                     
  # * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
  msn_loss <- pln_up_y[11] + A_decomp*Lc*ssl_loss*ssl_in_y[11]*
    ((1.0-r_stsc)/c2n_slow - 1.0/c2n_structural)
  
  msn_ss   <- msn_med/msn_loss
#}

pss.big[,3] <- 1:nrow(pss.big)
pss.big[,6] <- dgeom(seq(1,nrow(pss.big)*blckyr,by=blckyr),0.05)
pss.big[,8] <- rep(fsc_ss[1],nrow(pss.big))
pss.big[,9] <- rep(ssl_ss[1],nrow(pss.big))
pss.big[,10] <- rep(ssl_ss[1],nrow(pss.big))
pss.big[,11] <- rep(ssc_ss[1],nrow(pss.big))
pss.big[,13] <- rep(msn_ss[1],nrow(pss.big))
pss.big[,14] <- rep(fsn_ss[1],nrow(pss.big))
write.table(css.big,file=paste(out,s,"spin",site.coord[which(s==sites)],".css",sep=""),row.names=FALSE,append=FALSE,
            col.names=TRUE,quote=FALSE)
write.table(pss.big,file=paste(out,s,"spin",site.coord[which(s==sites)],".pss",sep=""),row.names=FALSE,append=FALSE,
            col.names=TRUE,quote=FALSE)
}