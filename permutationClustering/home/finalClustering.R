library("argparser")

p <- arg_parser("elements")
p <- add_argument(p, "matrixName", help="matrix count name")
p <- add_argument(p, "nCluster", help="number of clusters wanted")
p <- add_argument(p, "format", help="matrix format like csv, txt...")
p <- add_argument(p, "separator", help="matrix separator")
p <- add_argument(p, "clustering", help="Clulstering method: SIMLR tSne griph")
argv <- parse_args(p)

options(bitmapType='cairo')
Sys.setenv("DISPLAY"=":0.0")

matrixName=argv$matrixName
nCluster=as.numeric(argv$nCluster)
format=argv$format
separator=argv$separator
clusteringMethod=argv$clustering

if(separator=="tab"){separator2="\t"}else{separator2=separator}

cluster_p=sapply(list.files("./Permutation/",pattern="clusterB_*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]})
killedC=sapply(list.files("./Permutation/",pattern="killC*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]})
write.table(as.matrix(cluster_p,col.names=1),paste(matrixName,"_final_clusterP.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)
write.table(as.matrix(killedC,col.names=1),paste(matrixName,"_final_killedCell.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)
system(paste("cp ", matrixName,"_final_clusterP.",format, " ./data/",sep=""))
system(paste("cp ", matrixName,"_final_killedCell.",format, " ./data/",sep=""))

if(clusteringMethod == "griph"){
  pdf("hist.pdf")
  clusters=apply(cluster_p,2,FUN=function(x){max(x)})
  hist(clusters,xlab="nCluster",breaks=length(unique(cluster_p)))
  dev.off()
}
