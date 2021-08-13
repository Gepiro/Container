array <- commandArgs(trailingOnly = TRUE)
library(rCASC)
file=array[1]
gtf=array[2]
scannobyGtf(group="docker", file=paste(getwd(),file,sep="/"), gtf.name=gtf, biotype="protein_coding",mt=TRUE, ribo.proteins=TRUE, umiXgene=3, riboStart.percentage=0,riboEnd.percentage=100, mitoStart.percentage=0, mitoEnd.percentage=100, thresholdGenes=100)
