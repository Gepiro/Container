setwd("/home")
source("functions.R")
library("SIMLR")
library("argparser")
library(dplyr)
library("vioplot")
library("Publish")

p <- arg_parser("clusteringZero") 
p <- add_argument(p, "matrixName", help="matrix count name")
p <- add_argument(p, "nCluster", help="number of cluster desired")
p <- add_argument(p, "format", help="matrix format like csv, txt...")
p <- add_argument(p, "separator", help="matrix separator ")
p <- add_argument(p, "log10", help="1 or 0 if is matrix is already in log10 or if is not")
p <- add_argument(p, "clustering", help="Clulstering method: SIMLR tSne griph")
p <- add_argument(p, "perplexity", help="Number of close neighbors each point has")
p <- add_argument(p, "seed", help="Seed necessary for the reproducibility")

argv <- parse_args(p)

options(bitmapType='cairo')
Sys.setenv("DISPLAY"=":0.0")

matrixName=argv$matrixName
format=argv$format
separator=argv$separator
logTen=as.numeric(argv$log10)
clusteringMethod=argv$clustering
perplexity=as.numeric(argv$perplexity)
seed=as.numeric(argv$seed)

set.seed(seed)
dir.create(paste("./../scratch/",matrixName,sep=""))

if(clusteringMethod == "SIMLR" || clusteringMethod == "tSne"){
  nCluster=as.numeric(argv$nCluster)

  system(paste("cp ","/scratch","/",matrixName,".",format," ","/scratch","/",matrixName,sep=""))
  setwd(paste("./../scratch/",matrixName,sep=""))
  clustering(matrixName,nCluster,logTen,format,separator,clusteringMethod,perplexity)
  setwd("./../../../home")
  
  setwd(paste("./../scratch/",matrixName,"/",sep=""))
  
  silhouettePlot(matrixName,nCluster,format,separator)
  setwd("./../..")
  
}else{
  if(!is.null(argv$nCluster)){cat(paste("\nWARNING: nCluster with ",clusteringMethod," is suppose to be null\n"))}
  setwd(paste("./../scratch/",matrixName,"/",sep=""))

  nCluster=clustering(matrixName,nCluster=0,logTen,format,separator,clusteringMethod,perplexity)
  system("rm Lvis_*")
  system("rm ./../Lvis*")

  silhouettePlot(matrixName,nCluster,format,separator)
  setwd("./../..")
  
}

system(paste("echo ", nCluster, " > nCluster.txt", sep=""))
system(paste("cp nCluster.txt ./data", sep=""))
system("chmod -R 777 ./scratch")
system("chmod -R 777 ./data")
