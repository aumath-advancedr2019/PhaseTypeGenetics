.onLoad <- function(){
  
  expmpackage <-ifelse(nzchar(system.file(package = "expm")), TRUE, FALSE)
  partitionpackage <-ifelse(nzchar(system.file(package = "partitions")), TRUE, FALSE)
  
  if(expmpackage & partitionpackage){
    
    cat("This package requires the packages *expm* and *partition*.\n Both packages are installed on your device and will now be loaded...\n")
    
    library(partitions)
    library(expm)
    
  }else if(expmpackage & partitionpackage==FALSE){
    
    cat("This package requires the packages *expm* and *partition*. \n
It seems that the package *partition* is not installed on your device.\n
You need to install this package in order to be able to use all functions.\n
The package *exmp* is installed on your device and will now be loaded...\n")
    
    library(expm)
    
  }else if(expmpackage==FALSE & partitionpackage){
    
    cat("This package requires the packages *expm* and *partition*.\n
It seems that the package *expm* is not installed on your device.\n
You need to install this package in order to be able to use all functions.\n
The package *partitions* is installed on your device and will now be loaded...\n")
    
    library(partitions)
    
  }else{
    
    cat("This package requires the packages *expm* and *partition*.\n
It seems that none of the packages are installed on your device.\n
Please install both packages in order to be able to use all
functions.\n")
  }
}
