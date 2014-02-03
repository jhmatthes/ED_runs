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

for(s in 1:length(sites)){
  #Set directories
  dat.dir    <- paste("/projectnb/cheas/paleon/ED_runs/phase1a_runs/",sites[s],"/analy/",sep="")
  match.files <- grep("-Y-",list.files(dat.dir))
  files <- list.files(dat.dir)
  ann.files  <- files[match.files] #yearly files only
  
  #Get time window
  yeara <- as.numeric(strsplit(ann.files,"-")[[1]][3]) #first year
  yearz <- as.numeric(strsplit(ann.files,"-")[[length(ann.files)]][3]) #last year
  
  agb.pft <- lai.pft <- bsa.pft <- matrix(nrow=(yearz-yeara+1),ncol=length(pft))
  balive <- broot <- bleaf <- bsapa <- bsapb <- sfast <- sslow <- sstruc <- vector(length=(yearz-yeara+1))
  #loop over years and aggregate annual data
  
  for (y in yeara:yearz){
    cat(" - Reading file :",ann.files[y-yeara+1],"...","\n")
    now <- open.ncdf(paste(dat.dir,ann.files[y-yeara+1],sep=""))
    
    #Grab global & patch level variables.
    npoly     <- get.var.ncdf(now,'NPOLYGONS_GLOBAL')
    nsites    <- get.var.ncdf(now,'NSITES_GLOBAL')
    npatches  <- get.var.ncdf(now,'NPATCHES_GLOBAL')
    ncohorts  <- get.var.ncdf(now,'NCOHORTS_GLOBAL')
    ncohort.per.patch    <- get.var.ncdf(now,'PACO_N')
    patch.area.per.site  <- get.var.ncdf(now,'AREA')  
    site.area.per.poly   <- get.var.ncdf(now,'AREA_SI')
    ind.patch.per.site   <- get.var.ncdf(now,'SIPA_ID')
    ind.cohort.per.patch <- get.var.ncdf(now,'PACO_ID')
    lat       <- get.var.ncdf(now,'LATITUDE')
    lon       <- get.var.ncdf(now,'LONGITUDE')
    
    #Grab cohort level variables.
    ipft      <- get.var.ncdf(now,'PFT')
    dbh       <- get.var.ncdf(now,'DBH')
    nplant    <- get.var.ncdf(now,'NPLANT')
    lai       <- get.var.ncdf(now,'LAI_CO')
    hgt       <- get.var.ncdf(now,'HITE')
    agb       <- get.var.ncdf(now,'AGB_CO')
    bsa       <- get.var.ncdf(now,'BA_CO')
    
    #if any PFTs go extinct, make placeholders
    if(length(unique(ipft))<length(pft)){
      tmp  <- (length(pft)-length(unique(ipft)))
      ipft <- c(ipft,pft[!(pft %in% ipft)])
      agb  <- c(agb,rep(0,tmp))
      lai  <- c(lai,rep(0,tmp))
      bsa  <- c(bsa,rep(0,tmp))
    }
    
    agb.pft[(y-yeara+1),] <- tapply(agb,ipft,sum)
    lai.pft[(y-yeara+1),] <- tapply(lai,ipft,sum)
    bsa.pft[(y-yeara+1),] <- tapply(bsa,ipft,sum)
    
    #Average annual carbon pools [kgC/m2]
    balive[y-yeara+1]    <- sum(get.var.ncdf(now,'BALIVE')) #avg living biomass
    bleaf[y-yeara+1]     <- sum(get.var.ncdf(now,'BLEAF')) 
    broot[y-yeara+1]     <- sum(get.var.ncdf(now,'BROOT')) 
    bsapa[y-yeara+1]     <- sum(get.var.ncdf(now,'BSAPWOODA')) 
    bsapb[y-yeara+1]     <- sum(get.var.ncdf(now,'BSAPWOODB')) 
    sfast[y-yeara+1]     <- get.var.ncdf(now,'FAST_SOIL_C') 
    sslow[y-yeara+1]     <- get.var.ncdf(now,'SLOW_SOIL_C') 
    sstruc[y-yeara+1]    <- get.var.ncdf(now,'STRUCTURAL_SOIL_C')
       
    close.ncdf(now)
  }
  years <- as.character((yeara:yearz)-1000)
  year.date <- as.Date(years,"%Y")
#  png(paste(plot.path,sites[s],'_AGBbyPFT','.png',sep=''),width=900,height=600)
#  pdf(paste(plot.path,sites[s],"_spinup",sep=''))
  plot(year.date,agb.pft[,1],col=pft.cols[5],pch=16,ylim=range(agb.pft),
       xlab="spin-up date",ylab="Annual aboveground biomass [kg/m2]",
       main=paste(sites[s],": Spin-up",sep=""))
  for(p in 2:ncol(agb.pft)){
    points(year.date,agb.pft[,p],col=pft.cols[4+p],pch=16)
  }
  legend(year.date[2],max(agb.pft)-mean(agb.pft),pft.names[sort(unique(ipft))],col=pft.cols[5:10],pch=16)

  cpools1 <- rbind(balive,bleaf,broot,bsapa,bsapb)
  plot(year.date,cpools1[1,],col=pft.cols[1],pch=16,ylim=range(c(balive,bleaf,broot)),
       xlab="spin-up date",ylab="Annual mean C pool [kg/m2]",
       main=paste(sites[s],": Spin-up",sep=""))
  for(p in 2:5){
    points(year.date,cpools1[p,],col=pft.cols[p],pch=16)
  }
  names <- c("balive","bleaf","broot","bsapa","bsapb")
  legend(year.date[2],max(cpools1)-mean(cpools1),names,col=pft.cols[1:nrow(cpools1)],pch=16)
  
  cpools2 <- rbind(sfast,sslow,sstruc)
  plot(year.date,cpools2[1,],col=pft.cols[1],pch=16,ylim=range(c(sfast,sslow,sstruc)),
       xlab="spin-up date",ylab="Annual mean C pool [kg/m2]",
       main=paste(sites[s],": Spin-up",sep=""))
  for(p in 2:nrow(cpools2)){
    points(year.date,cpools2[p,],col=pft.cols[p],pch=16)
  }
  names <- c("sfast","sslow","sstruc")
  legend(year.date[2],max(cpools2)-mean(cpools2),names,col=pft.cols[1:nrow(cpools2)],pch=16)
  
#  dev.off()
  
  rm(dat.dir,yeara,yearz,files,ann.files,match.files)
}

