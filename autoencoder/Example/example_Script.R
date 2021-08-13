library(rCASC)
source("skeleton.R")
file=paste(getwd(),"setA.csv",sep="/")
dir.create(paste(getwd(),"scratch",sep="/"))
autoencoderA(group=c("sudo"))

