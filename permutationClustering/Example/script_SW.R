source("permutationClusteringSW.R") 
library(rCASC)
system("tar -zxvf testSCumi_mm10.tar.gz")
file=paste(getwd(),"testSCumi_mm10.csv",sep="/")
permutationClusteringSW(group="docker", scratch.folder=getwd(), file=file,nCluster=3, separator=",", logTen=0, clustering="tSne", perplexity=10 , seed=1111, rK=0)

