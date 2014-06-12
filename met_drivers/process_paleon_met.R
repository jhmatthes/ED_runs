#Process netCDF PalEON met files into HDF5 files for ED2:
#
#This works with files in the format site_metvar_year_month.nc (i.e. Ha1_lwdown_1850_01.nc)
#It loads the netCDF file, formats the data into HDF5 format, and renames variables and the date 
#to be in the ED2 HDF5 format with the correct dimensions.  
#
#It requires the rhdf5 library, which is not available on CRAN, but by can be installed locally:
#>source("http://bioconductor.org/biocLite.R")
#>biocLite("rhdf5")
#
#on GEO CLUSTER (local install of rhdf5): 
#install.packages("/usr4/spclpgm/jmatthes/zlibbioc_1.6.0.tar.gz",repos=NULL,type="source",lib="/usr4/spclpgm/jmatthes/")
#install.packages("/usr4/spclpgm/jmatthes/rhdf5_2.4.0.tar.gz",repos=NULL,type="source",lib="/usr4/spclpgm/jmatthes/")
#
#Jaclyn Hatala Matthes, 1/7/14
#jaclyn.hatala.matthes@gmail.com

library(ncdf,lib.loc="/usr4/spclpgm/jmatthes/")
library(rhdf5,lib.loc="/usr4/spclpgm/jmatthes/")
library(abind,lib.loc="/usr4/spclpgm/jmatthes/")

in.path  <- "/projectnb/cheas/paleon/met_regional/phase1a_met_drivers_v4/"
out.path <- "/projectnb/cheas/paleon/ED_runs/met_drivers/phase1a_met/bias_corr_v3/"
dir.create(file.path(out.path), showWarnings = FALSE)

#sites <- c("PBL","PDL","PHA","PHO","PMB","PUN")
sites <- c("PHA","PHO","PMB","PUN")
orig.vars <- c("lwdown","precipf","psurf","qair","swdown","tair","wind")
ed2.vars  <- c("dlwrf","prate","pres","sh","vbdsf","tmp","ugrd")
month.txt <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

for(s in 1:length(sites)){
  in.dir  <- paste(in.path,sites[s],"/",sep="")
  out.dir <- paste(out.path,sites[s],"/",sep="")
  dir.create(file.path(paste(out.path,sites[s],"/",sep="")), showWarnings = FALSE)
  print(sites[s])

  for(v in 1:length(orig.vars)){
    var.path <- paste(in.dir,orig.vars[v],"/",sep="")
    in.files <- list.files(var.path)
    dir.create(file.path(paste(out.path,sites[s],"/",ed2.vars[v],"/",sep="")), showWarnings = FALSE)
    print(orig.vars[v])    

    for(f in 1:length(in.files)){
      
      #open and read netcdf file
      nc.file <- open.ncdf(paste(var.path,in.files[f],sep=""))
      var     <- get.var.ncdf(nc.file,orig.vars[v])
      time    <- get.var.ncdf(nc.file,"time")
      lat     <- get.var.ncdf(nc.file,"lat")
      lon     <- get.var.ncdf(nc.file,"lon")
      close.ncdf(nc.file)
      
      var <- array(var,dim=c(length(var),1,1))                   
      
      #process year and month for naming
      filesplit <- strsplit(in.files[f],"_")
      year  <- as.numeric(filesplit[[1]][3])+1000 ###CAREFUL - I HAD TO ADD 1000 FOR PALEON
      monthsplit <- strsplit(filesplit[[1]][4],".nc")
      month <- monthsplit[[1]]
      month.num <- as.numeric(month)
      
      #write HDF5 file
      out.file <- paste(out.path,sites[s],"/",ed2.vars[v],"/",sites[s],"_",ed2.vars[v],"_",year,month.txt[month.num],".h5",sep="")
      h5createFile(out.file)
      h5write(var,out.file,ed2.vars[v])
      h5write(time,out.file,"time")
      h5write(lon,out.file,"lon")
      h5write(lat,out.file,"lat")
      #  H5Fclose(paste(out.path,site,"_dlwrf_",year,"_",month.txt[month.num],".h5",sep=""))
    }
  }
}