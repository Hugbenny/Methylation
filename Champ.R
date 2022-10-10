{library(ChAMP)
library(stringr)
print("Where are the .idat located ?")
# setwd(readLines(con="stdin",1))
setwd(readline())
# setwd("")
print("Which analysis do you want to do ? (450k/EPIC)")
# analysis<-readLines(con="stdin",1)
analysis<-readline()
# analysis<-"EPIC"
print("Do you want to exclude XY probes from analysis ? (TRUE/FALSE)")
# XY<-readLines(con="stdin",1)
XY<-readline()
# XY<-TRUE
myLoad<-champ.load(autoimpute=FALSE,detPcut=0.01,filterXY=XY,arraytype=analysis)
myLoad$pd$Slide<-as.character(as.factor(myLoad$pd$Slide))
dir.create(path="../Results")
dir.create(path="../Results/ChAMP")
save(myLoad, file = file.path("../Results/ChAMP/myLoad.rda"))
champ.QC(beta=myLoad$beta,pheno=myLoad$pd$Sample_Group,resultsDir="../Results/ChAMP/CHAMP_QCimages_Sample_Group/")
if(file.exists(file.path("../Results/ChAMP/myNorm.rda"))){
  print("Do you want to rerun normalisation ? (Yes/No)")
  # normali<-readLines(con="stdin",1)
  normali<-readline()
  if(normali == "Yes"){
    myNorm<-champ.norm(beta=myLoad$beta,arraytype=analysis,core=detectCores(),resultsDir="../Results/ChAMP/CHAMP_Normalization/")
    save(myNorm, file = file.path("../Results/ChAMP/myNorm.rda"))
    write.table(cbind(rownames(myNorm),data.frame(myNorm,row.names=NULL)),file = file.path("../Results/ChAMP/myNorm.csv"),sep = ",",row.names = F)
    betavalplot<-boxplot(myNorm,cex.axis = 0.8, las =2, ylab = "BetaValue")
    save(betavalplot,file = file.path("../Results/ChAMP/betavalueplot.png"))
    champ.SVD(beta = as.data.frame(myNorm),resultsDir="../Results/ChAMP/CHAMP_SVDimages/")
  if(normali == "No"){
    load(file.path("../Results/ChAMP/myNorm.rda"))
    }
  }
}
if(!(file.exists(file.path("../Results/ChAMP/myNorm.rda")))){
  myNorm<-champ.norm(beta=myLoad$beta,arraytype=analysis,core=detectCores(),resultsDir="../Results/ChAMP/CHAMP_Normalization/")
  save(myNorm, file = file.path("../Results/ChAMP/myNorm.rda"))
  write.table(cbind(rownames(myNorm),data.frame(myNorm,row.names=NULL)),file = file.path("../Results/ChAMP/myNorm.csv"),sep = ",",row.names = F)
  betavalplot<-boxplot(myNorm,cex.axis = 0.8, las =2, ylab = "BetaValue")
  save(betavalplot,file = file.path("../Results/ChAMP/betavalueplot.png"))
  champ.SVD(beta = as.data.frame(myNorm),resultsDir="../Results/ChAMP/CHAMP_SVDimages/")
}
print("Do you want to run batch correction ? (Yes/No)")
# batchcorr<-readLines(con="stdin",1)
batchcorr<-readline()
if(batchcorr == "Yes"){
  myCombat<-champ.runCombat()
  save(myCombat, file = file.path("../Results/ChAMP/myCombat.rda"))
  myNorm<-myCombat
  write.table(cbind(rownames(myNorm),data.frame(myNorm,row.names=NULL)),file = file.path("../Results/ChAMP/myNorm.csv"),sep = ",",row.names = F)
}
gc()

print("Do you want to run DMP/DMR analysis ? (Yes/No)")
# dmrdmpon<-readLines(con="stdin",1)
dmrdmpon<-readline()
while(dmrdmpon == "Yes"){
  
#Put name of controls from Sample_Groups
# print("How many different controls do you have ?")
# nbct<-readLines(con="stdin",1)
print("Put the name of the controls as in Sample_Group, separated by a space")
# controllist<-as.list(readLines(con="stdin",nbct))
controllist<-unlist(strsplit(readline()," "))
# controllist<-list("0mCT","0nCT")

#Put logFC threshold
print("What logFC threshold do you want ?")
# logFC<-readLines(con="stdin",1)
logFC<-readline()
# logFC<-0.05
dir.create(path = paste("../Results/ChAMP/DMP/DMP",logFC,sep="_"),recursive=T)
dir.create(path = paste("../Results/ChAMP/DMR/DMR",logFC,sep="_"),recursive=T)
DMPdir<-"../Results/ChAMP/DMP"
DMRdir<-"../Results/ChAMP/DMR"
DMPdirlog<-paste("../Results/ChAMP/DMP/DMP",logFC,sep="_")
DMRdirlog<-paste("../Results/ChAMP/DMR/DMR",logFC,sep="_")

#Loop to produce DMP and DMR for each patient
for(ct in controllist){
  for(i in myLoad$pd$Sample_Group){
    if(i!=ct && !(i %in% controllist)){
      myDMP<-champ.DMP(compare.group=c(i,ct),adjPVal = 0.05,arraytype=analysis)
      save(myDMP, file = file.path(DMPdir,paste("myDMP-",i,"-",ct,".rda",sep="")))
      j<-myDMP[[1]]
      j$cgid<-rownames(j)
      j<-j[which(abs(j$deltaBeta) >= logFC),]
      write.table(j,file = file.path(DMPdirlog,paste(i,"-",ct,sep="")),sep = ",",row.names = F)
      myDMR<-champ.DMR(compare.group=c(i,ct),cores=detectCores(),arraytype=analysis,method="ProbeLasso",minProbes=3,meanLassoRadius=375,
                       minDmrSep=300,minDmrSize=10,resultsDir="../Results/ChAMP/CHAMP_ProbeLasso/")
      k<-myDMR$ProbeLassoDMR
      k$Diff<-unlist(k[19]-k[18])
      k<-k[which(abs(k$Diff) >= logFC),]
      save(myDMR, file = file.path(DMRdir,paste("myDMR-",i,"-",ct,".rda",sep="")))
      write.table(k,file = file.path(DMRdirlog,paste(i,"-",ct,sep="")), sep = ",", row.names = F)
      # DMP.GUI(pheno=myLoad$pd$Sample_Group,cutgroupnumber=6)
      # DMR.GUI(pheno=myLoad$pd$Sample_Group,runDMP=FALSE,compare.group=c("0CT",i),arraytype=analysis)
      # CpG.GUI(arraytype=analysis)
      # QC.GUI(arraytype=analysis)
      gc()}}}

#DMP and DMR generation according to disease
# print("How many different phenotype do you have ?")
# nbph<-readLines(con="stdin",1)
print("Put the phenotype as in Sample_Pheno, separated by a space")
# phenolist<-as.list(readLines(con="stdin",nbph))
phenolist<-unlist(strsplit(readline()," "))
# phenolist<-list("FSHD2","BAMS")
print("What is the name of the control group as in Sample_Pheno ?")
# ct<-readLines(con="stdin",1)
ct<-readline()
# ct<-"0CT"


for(i in phenolist){
  if(i!=ct){
    myDMP<-champ.DMP(pheno=myLoad$pd$Sample_Pheno,compare.group=c(i,ct),adjPVal=0.05,arraytype=analysis)
    save(myDMP,file=file.path(DMPdir,paste("myDMP-",i,".rda",sep="")))
    j<-myDMP[[1]]
    j$cgid<-rownames(j)
    j<-j[which(abs(j$deltaBeta)>=logFC),]
    write.table(j,file=file.path(DMPdirlog,paste(i,"-CT",sep="")),sep=",",row.names=F)
    myDMR<-champ.DMR(pheno=myLoad$pd$Sample_Pheno,compare.group=c(i,ct),cores=detectCores(),arraytype=analysis,method="ProbeLasso",minProbes=3,meanLassoRadius=375,
                     minDmrSep=300,minDmrSize=10,resultsDir="../Results/ChAMP/CHAMP_ProbeLasso/")
    k<-myDMR$ProbeLassoDMR
    k$Diff<-unlist(k[19]-k[18])
    k<-k[which(abs(k$Diff) >= logFC),]
    save(myDMR,file=file.path(DMRdir,paste("myDMR-",i,".rda",sep="")))
    write.table(k,file=file.path(DMRdirlog,paste(i,"-CT",sep="")),sep=",",row.names=F)
    DMP.GUI(pheno=myLoad$pd$Sample_Pheno,cutgroupnumber=2)
    DMR.GUI(pheno=myLoad$pd$Sample_Pheno,runDMP=FALSE,compare.group=c("0CT",i),arraytype=analysis)
    CpG.GUI(arraytype=analysis)
    QC.GUI(arraytype=analysis)
    gc()}}

print("Do you want to rerun DMP/DMR analysis ? (Yes/No")
# dmrdmpon<-readLines(con="stdin",1)
dmrdmpon<-readline()
}
#Quick DMP analysis for all patients
# myDMP<-champ.DMP(arraytype=analysis)
# DMP.GUI(pheno=myLoad$pd$Sample_Group,cutgroupnumber=6)
# 
# myGSEA<-champ.GSEA(metho="gometh",arraytype=analysis,cores=detectCores())
# DMR.GUI(pheno=myLoad$pd$Sample_Group,runDMP=FALSE,compare.group=c("0CT",i),arraytype=analysis)
}
