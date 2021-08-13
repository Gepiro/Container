source("skeleton_permAnalysis.R")
 library(rCASC)
permAnalysisGriphSW(group="docker",scratch.folder=getwd(),file=paste(getwd(),"testSCumi_mm10.csv",sep="/"),nCluster=3, 
                          clustering.output=paste(getwd(),"testSCumi_mm10_clustering.output.csv",sep="/"), 
                          clusterP=paste(getwd(),"testSCumi_mm10_final_clusterP.csv",sep="/"),
                          kill=paste(getwd(),"testSCumi_mm10_final_killedCell.csv",sep="/"), sep=",",sp=0.8)
