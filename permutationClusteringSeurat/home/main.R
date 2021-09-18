 
 setwd("/home")
 
source("functions.R")
 library(Seurat)
 library("argparser")
 library(dplyr)
 library("vioplot")
 library("Publish")
 
p <- arg_parser("permutation")
p <- add_argument(p, "matrixName", help="matrix count name")
p <- add_argument(p, "percent", help="Percentage of cell removed for bootstrap algorithm ")
p <- add_argument(p, "format", help="matrix format like csv, txt...")
p <- add_argument(p, "separator", help="matrix separator ")
p <- add_argument(p, "logTen", help="1 or 0 if is matrix is already in log10 or if is not")
p <- add_argument(p, "pcaDimensions", help="PCA dimension for seurat first number")
p <- add_argument(p, "seed", help="Seed necessary for the reproducibility")
p <- add_argument(p, "sparse", help="Seed necessary for the reproducibility")
 p <- add_argument(p, "resolution", help="Seed necessary for the reproducibility")
 p <- add_argument(p, "index", help="Seed necessary for the reproducibility")


argv <- parse_args(p)


#argv=list()
#argv$matrixName="annotated_setPace_10000"
#argv$nPerm=2
#argv$permAtTime=2
#argv$percent=10
#argv$format="txt"
#argv$separator="tab"
#argv$log10=0
#argv$seed=111
#argv$pcaDimensions=5




options(bitmapType='cairo')
Sys.setenv("DISPLAY"=":0.0")
matrixName=argv$matrixName
percent=as.numeric(argv$percent)
format=argv$format
separator=argv$separator
logTen=as.numeric(argv$logTen)
 pcaDimensions=as.numeric(argv$pcaDimensions)
seed=as.numeric(argv$seed)
set.seed(seed)
sparse=argv$sparse
 resolution=argv$resolution
 index=argv$index
dir.create(paste("./../scratch/",matrixName,sep=""))
 


  setwd(paste("./../scratch/",matrixName,"/",sep=""))
nCluster=clustering(matrixName,percent,nCluster=0,logTen,format,separator,pcaDimensions, index)

setwd("./../../../home")
  setwd(paste("./../scratch/",matrixName,"/",sep=""))
  silhouettePlot(matrixName,nCluster,format,separator, index)

#dir.create("./../../data/Results")

#system("cp -r ./../* ./../../data/Results")
setwd("./../..")
#system("rm -r ./scratch/*")





system("chmod -R 777 ./scratch") 
system("chmod -R 777 ./data")
