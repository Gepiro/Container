source("permutationClusteringSWP.R") 
library(rCASC)
system("tar -zxvf testSCumi_mm10.tar.gz")
file=paste(getwd(),"testSCumi_mm10.csv",sep="/")
permutationClusteringSWP(group="docker", scratch.folder=getwd(), file=file, percent=10, nCluster=8, separator=",", logTen=0, clustering="SIMLR", perplexity=10 , seed=1111, rK=0)

