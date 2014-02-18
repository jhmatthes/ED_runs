#This file compiles the ensemble output from the SAS runs 
#Jaclyn Hatala Matthes, 2/18/14
#jaclyn.hatala.matthes@gmail.com

#Load libraries
ok = require(chron,lib.loc="/usr4/spclpgm/jmatthes/"); if (! ok) stop("Package chron is not available...")
ok = require(ncdf,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package ncdf is not available...")
ok = require(colorspace,lib.loc="/usr4/spclpgm/jmatthes/") ; if (! ok) stop("Package raster is not available...")

#Setup analysis file structure
base  <- "/projectnb/cheas/paleon/ED_runs/phase1a_spinup/"
sites <- c("PBL","PHA","PMB","PDL","PHO","PUN")
niter <- length(list.dirs(paste(base,sites[1],"/",sep=""),recursive=FALSE)) #iterations/site
yeara <- 1850 #actually is 850, but had to trick ED
yearz <- 1855 
pft   <- c(5,6,8,9,10,11) #set of PFTs used in analysis
sufx  <- "g01.h5"

#Pretty PFT plotting params
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

#Set directories
s <- 1
  #Set up storage
  atm_pres <- atm_vpd <- atm_pre <- atm_tmp <- atm_rnt <- atm_trsp <- atm_swat <- npp_site <- vector()
  cbal.pft <- cbalbr.pft <- PAR.pft <- lai.pft <- matrix(nrow=((yearz-yeara-1)*12+(12-montha+1)),
                                                         ncol=length(pft))
  fnpp <- matrix(ncol=8,nrow=(yearz-yeara+1))
  #monthly storage
  npp <- leaf.drop <- npp.leaf <- npp.croot <- npp.froot <- npp.sapwood <- npp.wood <- npp.seed <- matrix(ncol=12,nrow=(yearz-yeara+1))
  
  gb.pft <- lai.pft <- bsa.pft <- matrix(nrow=(yearz-yeara+1),ncol=length(pft))
  balive <- broot <- bleaf <- bsapa <- bsapb <- sfast <- sslow <- sstruc <- vector(length=(yearz-yeara+1))
  #loop over years and aggregate annual data
  
fsc_in_y <- ssc_in_y <- ssl_in_y <- vector()

for (y in yeara:yearz){
  
  #Make the file name. 
  year.now  <-sprintf("%4.4i",y)
  month.now <- sprintf("%2.2i",6)
  day.now   <- sprintf("%2.2i",1)
  hour.now  <- sprintf("%6.6i",0)
  
  for(n in 1:niter){
    iter        <- sprintf("%02i",n)
    dat.dir     <- paste(base,sites[s],"/spin",iter,"/histo/",sep="")
    file.now    <- paste(sites[s],iter,"spin","-S-",year.now,"-",month.now,"-",day.now,"-"
                         ,hour.now,"-",sufx,sep="")
    
    cat(" - Reading file :",file.now,"...","\n")
    now <- open.ncdf(paste(dat.dir,file.now,sep=""))
    
    #get soil C params for SAS
    fsc_in[n] <- get.var.ncdf(now,"FSC_IN")
    ssc_in[n] <- get.var.ncdf(now,"SSC_IN")
    ssl_in[n] <- get.var.ncdf(now,"SSL_IN")
    
#       #Grab cohort level variables.
#       ipft      <- get.var.ncdf(now,'PFT')
#       dbh       <- get.var.ncdf(now,'DBH')
#       nplant    <- get.var.ncdf(now,'NPLANT')
#       lai       <- get.var.ncdf(now,'LAI_CO')
#       hgt       <- get.var.ncdf(now,'HITE')
#       agb       <- get.var.ncdf(now,'AGB_CO')
#       bsa       <- get.var.ncdf(now,'BA_CO')
#       
#       ncohorts  <- get.var.ncdf(now,'NCOHORTS_GLOBAL')
#       ncohort.per.patch    <- get.var.ncdf(now,'PACO_N')
#       
#       #if any PFTs go extinct, make placeholders
#       if(length(unique(ipft))<length(pft)){
#         tmp  <- (length(pft)-length(unique(ipft)))
#         ipft <- c(ipft,pft[!(pft %in% ipft)])
#         agb  <- c(agb,rep(0,tmp))
#         lai  <- c(lai,rep(0,tmp))
#         bsa  <- c(ba,rep(0,tmp))
#       }
#       
#       agb.pft[(y-yeara+1),] <- tapply(agb,ipft,sum)
#       lai.pft[(y-yeara+1),] <- tapply(lai,ipft,sum)
#       bsa.pft[(y-yeara+1),] <- tapply(bsa,ipft,sum)
#     
#       #Average annual carbon pools [kgC/m2]
#       balive[y-yeara+1]    <- sum(get.var.ncdf(now,'BALIVE')) #avg living biomass
#       bleaf[y-yeara+1]     <- sum(get.var.ncdf(now,'BLEAF')) 
#       broot[y-yeara+1]     <- sum(get.var.ncdf(now,'BROOT')) 
#       bsapa[y-yeara+1]     <- sum(get.var.ncdf(now,'BSAPWOODA')) 
#       bsapb[y-yeara+1]     <- sum(get.var.ncdf(now,'BSAPWOODB')) 
#       sfast[y-yeara+1]     <- get.var.ncdf(now,'FAST_SOIL_C') 
#       sslow[y-yeara+1]     <- get.var.ncdf(now,'SLOW_SOIL_C') 
#       sstruc[y-yeara+1]    <- get.var.ncdf(now,'STRUCTURAL_SOIL_C')
#       
    close.ncdf(now)
  }
  fsc_in_y[y-yeara+1] <- mean(fsc_in,na.rm=TRUE)
  ssc_in_y[y-yeara+1] <- mean(ssc_in,na.rm=TRUE)
  ssl_in_y[y-yeara+1] <- mean(ssl_in,na.rm=TRUE)
}

  


