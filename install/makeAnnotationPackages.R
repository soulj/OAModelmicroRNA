library(ensembldb)
library(AnnotationHub)
library(AnnotationForge)



# Load the annotation resource.
ah <- AnnotationHub()

makeEnsemblDBPackage <- function(speciesName,version=103){
  
ahDb<- query(ah, pattern = c(speciesName, "EnsDb",version))[[1]]

package<-makeEnsembldbPackage(ensdb = dbfile(dbconn(ahDb)), version = "0.0.1",
                     maintainer = "Jamie Soul <jamie.soul@liverpool.ac.uk>",
                     author = "Jamie Soul")
return(package)
}

species<-c("Homo sapiens","Mus musculus")

#create the ensembldb packages
lapply(species,makeEnsemblDBPackage)

#install all the packages
install.packages(list.files(pattern = "EnsDb*"),repos = NULL)


