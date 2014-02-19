#This file compiles the ensemble output from the SAS runs 
#Jaclyn Hatala Matthes, 2/18/14
#jaclyn.hatala.matthes@gmail.com

#Load libraries
ok = require(chron,lib.loc="/usr4/spclpgm/jmatthes/"); if (! ok) stop("Package chron is not available...")
ok = require(ncdf,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package ncdf is not available...")
ok = require(colorspace,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")

#Setup analysis file structure
base  <- "/projectnb/cheas/paleon/ED_runs/phase1a_spinup/"
out   <- "/projectnb/cheas/paleon/ED_runs/SAS_spinup/phase1a_spinup/"
sites <- c("PBL","PHA","PMB","PDL","PHO","PUN")
niter <- length(list.dirs(paste(base,sites[1],"/",sep=""),recursive=FALSE)) #iterations/site 
pft   <- c(5,6,8,9,10,11) #set of PFTs used in analysis
dpm <- c(31,28,31,30,31,30,31,31,30,31,30,31)
sufx  <- "g01.h5"

#constants from ED2 for SAS solution for soil pools
fsc_loss <- 1.0/11.0
ssc_loss <- 1.0/100.2
ssl_loss <- 1.0/5.4
resp_temperature_increase <- 0.0757
resp_opt_water            <- 0.8938
resp_water_below_opt      <- 5.0786
Lc                        <- 0.049787
c2n_slow                  <- 10.0
c2n_structural            <- 150.0
r_stsc                    <- 0.3
soil_tempk <- c(278.5,279.3, 280.0, 277.3, 276.8, 277.2) #mean temp per site
temperature_limitation = resp_temperature_increase * exp(308.56 * (1./56.02 - 1./(soil_tempk-227.15)))
water_limitation <- exp((0.85 - resp_opt_water) * resp_water_below_opt)
A_decomp <- temperature_limitation * water_limitation

#First loop over analy files (faster than histo) to aggregate initial 
#.css and .pss files for each site
for(s in sites){
  
  #Set directories
  dat.dir    <- paste(base,s,"/spin01/analy/",sep="")
  match.files <- grep("-Y-",list.files(dat.dir))
  files <- list.files(dat.dir)
  ann.files  <- files[match.files] #yearly files only
  
  #Get time window
  yeara <- as.numeric(strsplit(ann.files,"-")[[1]][3]) #first year
  yearz <- as.numeric(strsplit(ann.files,"-")[[length(ann.files)]][3]) #last year

  #storage
  pss.big <- matrix(nrow=(yearz-yeara+1),ncol=14)
  colnames(pss.big) <- c("site","year","patch","dst","age","area","water","fsc","stsc","stsl",
                         "ssc","psc","msn","fsn")
  
  for (y in yeara:yearz){
    cat(" - Reading file :",ann.files[y-yeara+1],"...","\n")
    now <- open.ncdf(paste(dat.dir,ann.files[y-yeara+1],sep=""))
    
    #Grab variable to see how many cohorts there are
    ipft      <- get.var.ncdf(now,'PFT')
    
    #organize into .css variables
    css.tmp <- matrix(nrow=length(ipft),ncol=10)
    css.tmp[,1] <- rep(y,length(ipft))
    css.tmp[,2] <- rep(y-yeara+1,length(ipft))
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
    
    #if any PFTs go extinct, make placeholders for averaging
    if(length(unique(ipft))<length(pft)){
      tmp  <- (length(pft)-length(unique(ipft)))
      ipft <- c(ipft,pft[!(pft %in% ipft)])
      agb  <- c(agb,rep(0,tmp))
      lai  <- c(lai,rep(0,tmp))
      bsa  <- c(bsa,rep(0,tmp))
    }
    
    #save .pss variables
    pss.big[(y-yeara+1),1]  <- 1
    pss.big[(y-yeara+1),2]  <- y
    pss.big[(y-yeara+1),3]  <- y-yeara+1
    pss.big[(y-yeara+1),4]  <- 1
    pss.big[(y-yeara+1),5]  <- y-yeara+1
    pss.big[(y-yeara+1),6]  <- get.var.ncdf(now,"AREA")
    pss.big[(y-yeara+1),7]  <- 0.1
    pss.big[(y-yeara+1),8]  <- get.var.ncdf(now,"FAST_SOIL_C")
    pss.big[(y-yeara+1),9]  <- get.var.ncdf(now,"STRUCTURAL_SOIL_C")
    pss.big[(y-yeara+1),10] <- get.var.ncdf(now,"STRUCTURAL_SOIL_L")
    pss.big[(y-yeara+1),11] <- get.var.ncdf(now,"SLOW_SOIL_C")
    pss.big[(y-yeara+1),12] <- 0
    pss.big[(y-yeara+1),13] <- get.var.ncdf(now,"AVG_MSN")
    pss.big[(y-yeara+1),14] <- get.var.ncdf(now,"AVG_FSN")
    
    close.ncdf(now)
  }
  write.table(css.big,file=paste(out,s,"spin.css",sep=""),row.names=FALSE,append=FALSE,
              col.names=TRUE,quote=FALSE)
  write.table(pss.big,file=paste(out,s,"spin.pss",sep=""),row.names=FALSE,append=FALSE,
              col.names=TRUE,quote=FALSE)
}


#Second loop over histo files (much slower than analy) to aggregate soil inputs
#for steady-state solution

#storage
fsc_in_y <- ssc_in_y <- ssl_in_y <- fsn_in_y <- pln_up_y <- vector()
fsc_in_m <- ssc_in_m <- ssl_in_m <- fsn_in_m <- pln_up_m <-  vector()

for(s in sites){
  dat.dir    <- paste(base,s,"/spin01/histo/",sep="")
  match.files <- grep("-S-",list.files(dat.dir))
  files <- list.files(dat.dir)
  mon.files  <- files[match.files] #monthly files only
  
  #Get time window
  yeara <- as.numeric(strsplit(mon.files,"-")[[1]][3]) #first year
  yearz <- as.numeric(strsplit(mon.files,"-")[[length(mon.files)]][3]) #last year
  montha <- as.numeric(strsplit(mon.files,"-")[[1]][4]) #first month
  monthz <- as.numeric(strsplit(mon.files,"-")[[length(mon.files)]][4]) #first month
  
  for (y in yeara:yearz){
    
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
      
      dat.dir     <- paste(base,s,"/spin01/histo/",sep="")
      file.now    <- paste(s,"01spin","-S-",year.now,"-",month.now,"-",day.now,"-"
                           ,hour.now,"-",sufx,sep="")
      
      cat(" - Reading file :",file.now,"...","\n")
      now <- open.ncdf(paste(dat.dir,file.now,sep=""))
      
      fsc_in_m[m-month.begin+1] <- get.var.ncdf(now,"FSC_IN")*dpm[m] #kg/(m2*day) --> kg/(m2*month)
      ssc_in_m[m-month.begin+1] <- get.var.ncdf(now,"SSC_IN")*dpm[m]
      ssl_in_m[m-month.begin+1] <- get.var.ncdf(now,"SSL_IN")*dpm[m]
      fsn_in_m[m-month.begin+1] <- get.var.ncdf(now,"FSN_IN")*dpm[m]
      pln_up_m[m-month.begin+1] <- get.var.ncdf(now,"TOTAL_PLANT_NITROGEN_UPTAKE")*dpm[m]
      
    }
    fsc_in_y[y-yeara+1] <- sum(fsc_in_m,na.rm=TRUE)
    ssc_in_y[y-yeara+1] <- sum(ssc_in_m,na.rm=TRUE)
    ssl_in_y[y-yeara+1] <- sum(ssl_in_m,na.rm=TRUE)
    fsn_in_y[y-yeara+1] <- sum(fsn_in_m,na.rm=TRUE)
    pln_up_y[y-yeara+1] <- sum(pln_up_m,na.rm=TRUE)
    
    close.ncdf(now)
  }
  
  #Calculate steady-state soil pools
  fsc_ss <- median(fsc_in_y)/(fsc_loss * A_decomp)
  ssc_ss <- median(ssc_in_y)/(ssc_loss * A_decomp)
  ssl_ss <- median(ssl_in_y)/(ssl_loss * A_decomp * Lc)
  fsn_ss <- median(fsn_in_y)/(fsc_loss * A_decomp)
  
  #fast_N_loss + slow_C_loss
  msn_med  <- fsc_loss*A_decomp*median(fsn_in_y)+ (ssc_loss * A_decomp)/c2n_slow 
  
  #ED2: csite%mineralized_N_loss  = csite%total_plant_nitrogen_uptake(ipa)             
  # + csite%today_Af_decomp(ipa) * Lc * K1 * csite%structural_soil_C(ipa)                     
  # * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
  msn_loss <- median(pln_up_y) + A_decomp*Lc*ssl_loss*median(ssl_in_y)*
    ((1.0-r_stsc)/c2n_slow - 1.0/c2n_structural)
  
  msn_ss   <- msn_med/msn_loss
}
