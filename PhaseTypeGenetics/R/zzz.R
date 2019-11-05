.onLoad <- function(libname=NULL, pkgname){
  
  expmpackage <-ifelse(nzchar(system.file(package = "expm")), TRUE, FALSE)
  partitionspackage <-ifelse(nzchar(system.file(package = "partitions")), TRUE, FALSE)
  
  if(expmpackage & partitionspackage){
    
    cat("This package requires the packages *expm* and *partitions*.\n Both packages are installed on your device and will now be loaded...\n")
    
    library(expm)
    
  }else if(expmpackage & partitionspackage==FALSE){
    
    cat("This package requires the packages *expm* and *partitions*. \n
It seems that the package *partitions* is not installed on your device.\n
You need to install this package in order to be able to use all functions.\n
The package *exmp* is installed on your device and will now be loaded...\n")
    
    library(expm)
    
  }else if(expmpackage==FALSE & partitionspackage){
    
    cat("This package requires the packages *expm* and *partitions*.\n
It seems that the package *expm* is not installed on your device.\n
You need to install this package in order to be able to use all functions.\n")
    
  }else{
    
    cat("This package requires the packages *expm* and *partitions*.\n
It seems that none of the packages are installed on your device.\n
Please install both packages in order to be able to use all
functions.\n")
  }
}
