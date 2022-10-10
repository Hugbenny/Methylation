{
#Converting fastq to fasta
library(seqinr)
library(ShortRead)
print("Where are the fastq located ?")
fastqfiles<-readline()
input <- dir(fastqfiles, pattern = "fastq.gz$")
output <- paste0(fastqfiles, "\\",  gsub("fastq.gz","fa", input))
for(i in seq(along = input)) writeFasta(readFastq(fastqfiles, input[i]), output[i])

# Launching BiQ Analyzer
system('java -jar "C:\\Program Files\\BiQ Analyzer HiMod\\BiQ5HiMod.jar"')

# Selection of best quality alignments
print("Where is the BiQ directory ?")
setwd(readline())
for(directory in list.dirs(path=getwd(),recursive=F)){
  for(result in list.files(paste(directory,"DR1",sep="\\"),pattern="results_Bisulfite.tsv")){
    file.copy(file.path(paste(directory,"DR1",result,sep="\\")),getwd())
    file.rename("results_Bisulfite.tsv",paste(basename(directory),".tsv",sep=""))
  }
}
total<-data.frame()
for(result in list.files(getwd(),pattern=".tsv",full.names=T)){
  temp<-read.csv(result,sep="\t")
  temp<-temp[temp$Sequence_Identity > 0.9 & temp$Conversion > 0.9,]
  total<-rbind(temp,total)
}
id<-data.frame(total$ID)
colnames(id)[1]<-"ID"
cpg<-data.frame(total$Mean_5mC...5hmC.pattern)
cpg$total.Mean_5mC...5hmC.pattern<-as.character(cpg$total.Mean_5mC...5hmC.pattern)
colnames(cpg)[1]<-"CpG"
sample<-data.frame(total$Sample)
colnames(sample)[1]<-"Samples"

# Calculating methylation levels and coverage
library(tidyr)
lol<-separate(data=cpg, col=CpG, into = c("CG_1", "CG_2", "CG_3", "CG_4",
                                              "CG_5" ,"CG_6", "CG_7", "CG_8", "CG_9",
                                              "CG_10", "CG_11", "CG_12", "CG_13", "CG_14", "CG_15",
                                              "CG_16", "CG_17", "CG_18", "CG_19", "CG_20", "CG_21",
                                              "CG_22", "CG_23", "CG_24", "CG_25", "CG_26", "CG_27",
                                              "CG_28", "CG_29", "CG_30","CG_31"),
              sep = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),
              remove = T,
              convert = F, extra = "warn", fill = "warn")
tot<-cbind(id,lol,sample)
tot$Reference<-"DR1"
tot<-tot[,c(1:32,34,33)]
write.table(tot, "DR1.tsv", sep = "\t", row.names = F)
print("Put the name of all samples separated by a space")
ech<-unlist(strsplit(readline()," "))
Sequence="DR1"

# Variables generales
{resultat.methylation<-matrix(nrow=1, ncol=length(ech), dimnames = list(Sequence, ech))
  resultat.couverture<-matrix(nrow=1, ncol=length(ech), dimnames = list(Sequence, ech))
  resultat.total<-matrix(nrow=2, ncol=length(ech), dimnames = list(c("Methylation", "Couverture") , ech))
  valeurs<-c("0", "1", "x")
  # Chargement des donnees
  fichier<-paste(Sequence,".tsv", sep="")
  print(paste("Chargement des donnees du fichier :", fichier))
  DATA<-read.table(fichier, sep="\t", header=TRUE)
  DATA.percentage<-matrix(nrow=length(ech), ncol=(length(DATA)-3), byrow = TRUE, dimnames = list(ech, 1:(length(DATA)-3)))
  DATA.couverture<-matrix(nrow=length(ech), ncol=(length(DATA)-3), byrow = TRUE, dimnames = list(ech, 1:(length(DATA)-3)))
  
  for (j in 1:length(ech)) {
    #Numeration des Methyl - Non Methyl - Non aligne
    DATA.table<-DATA[DATA$Sample==ech[j],2:(length(DATA)-2)]
    DATA.level<-table(factor(DATA.table[,1], levels=valeurs))
    for (i in 2:ncol(DATA.table)) {
      DATA.level<-rbind(DATA.level,table(factor(DATA.table[,i], levels=valeurs)))
    }
    #Pourcentage de methyl-C
    DATA.percentage[j,]<-DATA.level[,2]/(DATA.level[,1]+DATA.level[,2])
    resultat.methylation[Sequence,ech[j]]<-sum(DATA.level[,2])/sum(DATA.level[,1:2])
    #Pourcentage de non aligne
    DATA.couverture[j,]<-(DATA.level[,1]+DATA.level[,2])/(DATA.level[,1]+DATA.level[,2]+DATA.level[,3])
    resultat.couverture[Sequence,ech[j]]<-sum(DATA.level[,1:2])/sum(DATA.level[,1:3])
    val<-round(resultat.methylation[Sequence,j]*100,2)
    titre=paste(Sequence," - Barcode ", ech[j], " - Methylation Moyenne ", val,"%", sep = "")
    # print(titre)
  }
  
  resultat.total["Methylation",]<-resultat.methylation[Sequence,]
  resultat.total["Couverture",]<-resultat.couverture[Sequence,]
  
  write.table(resultat.total, paste("CpG-Global-",Sequence,".csv", sep=""), sep="\t",dec=",", na="NA")
  write.table(DATA.couverture, paste("CpG-Couverture-",Sequence,".csv", sep=""), sep="\t",dec=",", na="NA")
  write.table(DATA.percentage, paste("CpG-Meth-",Sequence,".csv", sep=""), sep="\t",dec=",", na="NA")
}


{fichier<-paste("Global-CpG-",Sequence,".pdf",sep = "")
  
  cent<-matrix(1,nrow=1,ncol=length(ech))
  pdf(file=fichier, width=(length(ech)+2), height=5)
  Fenetre<-par(mfrow=c(1,2))
  titre<-paste(Sequence," - % Methylation", sep = "")
  barplot(cent, col="white", ylim=c(0,1), space=0.1, las=1, cex.axis=1.5, axisnames = FALSE, sub=titre, cex.sub = 1)
  barplot(resultat.methylation[Sequence,], col="black", ylim=c(0,1), add=TRUE, space=0.1, las=1, cex.axis=1.5, axisnames = FALSE)
  titre<-paste(Sequence," - % Couverture", sep = "")
  barplot(cent, col="white", ylim=c(0,1), space=0.1, las=1, cex.axis=1.5, axisnames = FALSE, sub=titre, cex.sub = 1)
  barplot(resultat.couverture[Sequence,], col="black", ylim=c(0,1), add=TRUE, space=0.1, las=1, cex.axis=1.5, axisnames = FALSE)
  par(Fenetre)
  dev.off()
  
  
  #Graph Pourcentage de methyl-C
  fichier<-paste("Meth-CpG-",Sequence,".pdf",sep = "")
  
  pdf(file=fichier, width=((length(DATA.percentage[1,])+2)), height=length(ech)*5)
  Fenetre<-par(mfrow=c(length(ech),1))
  cent<-matrix(1,nrow=length(DATA.percentage[1,]),ncol=1)
  for (j in 1:length(ech)) {
    val<-round(resultat.methylation[Sequence,ech[j]]*100, 2)
    titre=paste(Sequence," - Barcode ", ech[j], " - Methylation Moyenne ", val,"%", sep = "")
    # print(titre)
    barplot(cent[,1], col="white", ylim=c(0,1), space=0.1, las=1, cex.axis=1, sub=titre, cex.sub = 2)
    barplot(DATA.percentage[ech[j],], col="black", ylim=c(0,1), add=TRUE, space=0.1, las=1, cex.axis=1, axisnames = FALSE)
  }
  par(Fenetre)
  dev.off()
  
  #Graph Pourcentage de couverture
  fichier<-paste("Couverture-CpG-",Sequence,".pdf",sep = "")
  
  pdf(file=fichier, width=((length(DATA.percentage[1,])+2)/2), height=length(ech)*5)
  Fenetre<-par(mfrow=c(length(ech),1))
  cent<-matrix(1,nrow=length(DATA.percentage[1,]),ncol=1)
  for (j in 1:length(ech)) {
    val<-round(resultat.couverture[Sequence,ech[j]]*100, 2)
    titre=paste(Sequence," - Barcode ", ech[j], " - Couverture Moyenne ", val,"%", sep = "")
    print(titre)
    barplot(cent[,1], col="white", ylim=c(0,1), space=0.5, las=1, cex.axis=1, sub=titre, cex.sub = 2)
    barplot(DATA.couverture[ech[j],], col="black", ylim=c(0,1), add=TRUE, space=0.5, las=1, cex.axis=1, axisnames = FALSE)
  }
  par(mar=c(1,1,1,1))
  dev.off()
  
  fichier<-paste("Histo-CpG-",Sequence,".pdf",sep = "")
  
  pdf(file=fichier, width=5, height=length(ech)*5)
  Fenetre<-par(mfrow=c(length(ech),1))
  cent<-matrix(1,nrow=length(DATA.percentage[1,]),ncol=1)
  for (j in 1:length(ech)) {
    val<-round(resultat.methylation[Sequence,ech[j]]*100, 2)
    titre=paste(Sequence," - Barcode ", ech[j], " - Methylation Moyenne ", val,"%", sep = "")
    print(titre)
    hist(DATA.percentage[ech[j],], freq=FALSE, main="", xlab="", ylab="", sub=titre, cex.axis=1.5, las=1, cex.lab=1, xlim=c(0,1), cex.sub = 1)
    abline(v=resultat.methylation[Sequence,ech[j]], lty="longdash")
    }
  par(Fenetre)
  dev.off() 
  }
}
