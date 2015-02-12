



##Define the packages needed
cran.packages <- c("ecodist", "devtools", "MASS", "pROC", "pheatmap", "pensim", "caret", "logging", "optparse", "knitr", "caret", "multtest", "gplots", "ggplot2", "mvtnorm")

bioc.packages <- c("affy", "genefilter", "biomaRt", "limma")

needed.packages <- c(cran.packages, bioc.packages)

if (!require(BiocInstaller))
    stop("You need to install Bioconductor, which includes BiocInstaller.")


for (pkg in needed.packages){
    if(!require(package=pkg, character.only=TRUE)){
        print(paste("Need to install", pkg)) 
        biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
    }
}

if (!require( ClassDiscovery )){
    source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
    oompaLite()
}

if( !require(princyr) || 
   package_version( sessionInfo(package="princyr")$otherPkgs$princyr$Version ) < 
            package_version("0.3.0")){
    library(devtools)
    install_github("princyr", username="lwaldron")
    library(princyr)
}


if( !require(LeviRmisc) || 
   package_version( sessionInfo(package="LeviRmisc")$otherPkgs$LeviRmisc$Version ) < 
            package_version("0.16.1")){
    library(devtools)
    install_url("https://bitbucket.org/lwaldron/curatedovariandata/downloads/LeviRmisc_0.16.1.tar.gz")
}




