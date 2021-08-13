args = commandArgs(trailingOnly=TRUE)

setwd("/data")

matrixName=args[1]
format=args[2]
separator=args[3]

if(separator=="tab"){separator2="\t"}else{separator2=separator}

mtx <- read.table(file=paste(matrixName,".",format,sep=""), sep=separator2, header=TRUE,row.names=1)
tmp = rownames(mtx)
i=1
while(i<=length(tmp)){
  tmp[i]=strsplit(tmp[i],":")[[1]][2]
  if(tmp[i]=="NA"){
    tmp=tmp[-c(i)]
    mtx=mtx[-c(i),]
  }else{i=i+1}
}
rownames(mtx)=make.unique(tmp,sep="_")
write.table(mtx,"output.csv",sep=separator2,col.names=NA)
