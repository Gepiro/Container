source("skeleton_finalClustering.R") 
library(rCASC)
system("tar -zxvf testSCumi_mm10.tar.gz")
file=paste(getwd(),"testSCumi_mm10.csv",sep="/")
permutationFolder=paste(getwd(),"/Permutation", sep="")
finalClustering(group="docker", permutation.folder=permutationFolder, file=file, nCluster=3, separator=",", clustering="tSne")

