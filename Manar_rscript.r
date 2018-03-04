library(gplots)
memory.limit(size = 800000)
hmcols<-colorRampPalette(c("cyan4","tomato3","black","cyan3","darkred"))(5)
hmcols<-colorRampPalette(c("blue2","tomato3","black"))(3)
#setwd("C:/Users/Owner/Google Drive/RamLab/script/CS/Hagar/analysed_to_clusters")
setwd("C:/Users/Owner/Google Drive/RamLab/script/CS/Hagar/ManarMiseq")
f=list.files()
h=read.table(f[11])
h$pcG=h[,3]/(h[,2]+h[,3])
h$pcC=h[,5]/(h[,4]+h[,5])
hh=h[h[,1]>91,]
par(mfrow=c(2,1))
hist(hh$pcG,breaks=50,col="blue",xlim=c(0,1))
hist(hh$pcC,breaks=50,col="red",xlim=c(0,1))

hist(hh[,3],breaks=50,col="blue")
hist(hh[,5],breaks=50,col="red")


G=hh[hh[,3]>0 & hh[,5]==0,]
C=hh[hh[,5]>0 & hh[,3]==0,]
N=hh[hh[,5]==0 & hh[,3]==0,]
H=hh[hh[,5]>0 & hh[,5]<11 & hh[,3]>0 & hh[,3]<11,]

nrow(H)+nrow(G)+nrow(C)+nrow(N)
nrow(hh)

ref="GATTGGACATTCGGAAGAGGGCCCGCCTTCCCTGGGGAATCTCTGCGCACGCGCAGAACGCTTCGACCAATGAAAACACAGGAAGCCGTCCGCGCAACCGCGTTGCGTCACTTCTGCCGCCCCTGTTTCAAGGGATAAGAAACCCTGCGACAAAACCTCCTCCTTTTCCAAGCGGCTGCCGAAGATGGCGGAGGTGCAGGTATGGGCTCCGCGCGGGCCGGGGCGGCAAGGGGCCGGGTGGGATCCAGG"
temp=unlist(strsplit(as.vector(ref), split=""))  
cc=rep(0,length(temp))
for (n in c(1:length(temp))){
  if (temp[n]=="G"){cc[n]=2}
  if (temp[n]=="C"){cc[n]=3}
}

G1=G$V6
mat=matrix(0,nrow=length(G1),ncol=249)
for (n in c(1:length(G1))){
  temp=as.numeric(unlist(strsplit(as.vector(G1[n]), split=",")))  
  mat[n,temp]=1
}

temp=sample(c(1:length(G1)),10000)
matSamp=mat[temp,]
mmg=matSamp[,cc==2]
hmcols<-colorRampPalette(c("blue2","tomato3"))(2)
hmcols<-colorRampPalette(c("grey26","tomato3"))(2)
#pdf(paste(f[n],".1693AllAntiblue.clustered.pdf",sep="")) 
c=heatmap.2(mmg,col = hmcols,scale="none",trace="none",Colv = F,dendrogram = "row")
setwd("C:/Users/Owner/Google Drive/RamLab/script/CS/Hagar/ManarMiseq/240118")
#write.table(mmg,file="Sense_RPL13a_B200.txt",quote=F,col.names=F,row.names=F,sep="\t")   
mmg=as.matrix(read.table("Sense_RPL13a_B200.txt"))


C1=C$V7
matC=matrix(0,nrow=length(C1),ncol=249)
for (n in c(1:length(C1))){
  temp=as.numeric(unlist(strsplit(as.vector(C1[n]), split=",")))  
  matC[n,temp]=1
}

temp=sample(c(1:length(C1)),10000)
matSampC=matC[temp,]
mmc=matSampC[,cc==3]
hmcols<-colorRampPalette(c("blue2","tomato3"))(2)
#pdf(paste(f[n],".1693AllAntiblue.clustered.pdf",sep="")) 
c=heatmap.2(mmc,col = hmcols,scale="none",trace="none",Colv = F,dendrogram = "row")



mat=matrix(0,nrow=nrow(hh),ncol=249)
for (n in c(1:nrow(hh))){
  tempG=as.numeric(unlist(strsplit(as.vector(hh[n,"V6"]), split=","))) 
  tempC=as.numeric(unlist(strsplit(as.vector(hh[n,"V7"]), split=","))) 
  mat[n,tempG]=2
  mat[n,tempC]=3
}

temp=sample(c(1:nrow(hh)),10000)
mat=mat[temp,]
mm=mat[,cc>0]
hmcols<-colorRampPalette(c("grey","tomato3","green"))(3)
#pdf(paste(f[n],".1693AllAntiblue.clustered.pdf",sep="")) 
c=heatmap.2(mm,col = hmcols,scale="none",trace="none",Colv = F,dendrogram = "row")






for (n in c(1:5)){
  h=read.table(f[n])
  h=as.matrix(h)
  print(dim(h))
  z=apply(h,2,sum)
  h=h[,z>0]
  print(dim(h))
  h[h==9]=3
  h[h==7]=3
  count=apply(h,1,function(x){
    length(x[x==3])
  })
  h=h[count<10,]
  #h[h==2]=1
  #z=apply(h,2,min)
  #h=h[,z>1]
  bad=apply(h,1,function(x){
    table(x==2)[2]>0 & table(x==5)[2]>0
  })
  hh=h[is.na(bad),]
  sense=apply(hh,1,function(x){
    table(x==5)[2]>0
  })
  table(sense)
  anti=apply(hh,1,function(x){
    table(x==2)[2]>0 
  })
  table(anti)
  S=hh[!is.na(sense),]
  A=hh[!is.na(anti),]
  print(dim(S))
  print(dim(A))
  #hh=h[,c(77:1)]
  #print(dim(hh))
  #count=apply(hh,1,function(x){
  #  length(x[x==3])
  #})
  #hh=hh[count<10,]
  count=apply(A,1,function(x){
    length(x[x==6])
  })
  A=A[count<10,]
  count=apply(S,1,function(x){
    length(x[x==6])
  })
  S=S[count<10,]
  #print(dim(hh))
  #nn=nrow(hh)
  #if(nn>10000){
  #  int=sample(nn,10000)
  #  hh=hh[int,]
  #}
  #print(dim(hh))
  z=apply(A,2,min)
  A=A[,z!=4]
  A[A==6]=3
  z=apply(S,2,min)
  S=S[,z!=1]
  S[S==3]=6
  setwd("C:/Users/Owner/Google Drive/RamLab/script/CS/Hagar/analysed_to_clusters")
  pdf(paste(f[n],".9150antisense.clustered.pdf",sep="")) 
  a=heatmap.2(S[,1:100],col = hmcols,scale="none",trace="none",Colv = F,dendrogram = "row")
  dev.off()
  hhh=S[rev(a$rowInd),]
  write.table(hhh,file=paste(f[n],".9150antisense.clustered.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")   
}

h=read.table("Hagar2_sense8_conv_tab.txt.361.clustered.txt")
SS=rbind(as.matrix(h),S[,77:1])
h=read.table("Hagar2_anti8_conv_tab.txt.420.clustered.txt")
AA=rbind(as.matrix(h),A[,100:1])

hmcols<-colorRampPalette(c("blue","tomato3","black"))(3)
hmcols<-colorRampPalette(c("bisque4","tomato3","black"))(3)

pdf(paste(f[n],".1693AllAntiblue.clustered.pdf",sep="")) 
a=heatmap.2(AA,col = hmcols,scale="none",trace="none",Colv = F,dendrogram = "row")
dev.off()
hhh=AA[rev(a$rowInd),]
write.table(hhh,file=paste(f[n],".1693Allanti.clustered.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")   


h=read.table("Hagar_sense_conv_tab.txt.941.clustered.txt")

hmcols<-colorRampPalette(c("grey26","tomato3","black"))(3)
pdf("Mod_Hagar_sense_conv_tab.txt.341grey.clustered.pdf") 
a=heatmap.2(h[,80:1],col = hmcols,scale="none",trace="none",Colv = F,dendrogram = "row")
dev.off()
hhh=h[rev(a$rowInd),80:1]
write.table(hhh,file=paste(f[n],".341grey.clustered.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")   


hc=hclust(dist(1-hh),method= "average")
heatmap.2(hh[hc$order,],col = hmcols,scale="none",trace="none",Colv = F,Rowv = F,dendrogram = "none")

