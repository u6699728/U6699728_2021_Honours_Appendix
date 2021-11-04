library(gplots)
library(RColorBrewer)
library(dplyr)
library(corrgram)
library(Rsubread)
library(limma)
library(Glimma)
library(edgeR)


mydist<-function(c) {
  
  require(amap)
  Dist(c,method="pearson")
}
myclust<-function(c) { hclust(c,method='ward.D') }


makeRLEplot<-function(M,main,ylim,cols,cex){
  medians=apply(M,1,median)
  newmatrix=M-medians
  boxplot(newmatrix,main=main, col=cols, las=2,cex.name=cex,range=0,ylim=ylim,cex.axis=0.8)
  abline(h=0,col="grey",lty=1)
  #  medians
}

RNAseq_file = c("C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/Samtools/L1/L1.sorted.bam",
                "C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/Samtools/L2/L2.sorted.bam",
                "C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/Samtools/L3/L3.sorted.bam",
                "C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/Samtools/WT1/WT1.sorted.bam",
                "C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/Samtools/WT2/WT2.sorted.bam",
                "C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/Samtools/WT3/WT3.sorted.bam")

annotation_file = "C:/Users/u6699728/work/RNAseq_analysis/Attempts/Fourth_analysis_24_6_21/29_9_21_novel_analysis/Merged_gtf.gtf"


#featureCounts -T 36 -s 0 -M -P -a ${ANNOTATION}.gff -o ${OUTPUT} ${INPUT}.sorted.bam

Counted_files = featureCounts(files=RNAseq_file, annot.ext=annotation_file, isPairedEnd = TRUE, nthreads = 36, isGTFAnnotationFile = TRUE)
Counted_files

counted_information1 = Counted_files[1]
counted_information1
write.csv(counted_information1, "Counted_information1.csv")

counted_information2 = Counted_files[2]
counted_information2
write.csv(counted_information2, "Counted_information2.csv")

count.stats = Counted_files[4]
count.stats
write.csv(count.stats, "Count.stats.csv")

#########
s.counts1 <- read.delim("Count.stats.txt", header = T, sep="\t")
head(s.counts1)

dim(s.counts1)

map_stat1<-s.counts1[rowSums(s.counts1[,2:7])!=0,]
dim(map_stat1)

rownames(map_stat1)<-s.counts1$Status[rowSums(s.counts1[,2:7])!=0]
colnames(map_stat1)<-gsub(".sorted.bam","", colnames(s.counts1))
colnames(map_stat1)<-gsub("bam.","", colnames(map_stat1))


par(mfrow=c(1,1),mar=c(12,5,3,1))
barplot2(as.matrix(map_stat1[1:3,2:7]),las=2,cex.axis=0.6,cex.names=0.6,col=c("steelblue","steelblue4","violetred"),main="Sequencing reads summary")
legend("top",legend =rownames(map_stat1)[1:3],cex=0.6, text.col=c("steelblue","steelblue4","violetred"))

##########

RNA_data <- read.csv("Merged_information1.csv")
head(RNA_data)

dim(RNA_data)

g.inf <- RNA_data[,c(1,6)]
head(g.inf)

counts1<-RNA_data[,7:12]
colnames(counts1)<-c("L1", "L2", "L3", "WT1", "WT2", "WT3")
colnames(counts1)

s.idx <-order(substr(colnames(counts1),start=2,stop=2))
counts1.o <-counts1[,s.idx]
lCPM<-cpm(counts1.o,log=TRUE)
cols<-c(rep("red",3),rep("purple",3))

par(mfrow=c(1,1))        
plotMDS(lCPM,col=cols,main="MDS plot log2(CPM)")
legend("bottomleft",legend = c("1: MIM159","2: WT"),text.col=c("red","purple"),cex=0.6)

thresh <- lCPM > 0.5                                                                     #retain CPM>0.5
head(thresh)

colSums(thresh)

keep <- rowSums(thresh) >= 3
summary(keep) 

counts.keep <-counts1.o[keep,]
head(counts.keep)

dim(counts.keep)

g.inf.keep<-g.inf[keep,]
dim(g.inf.keep)

lCPM.keep<-cpm(counts.keep,log=TRUE)
par(mfrow=c(1,2),mar=c(5,5,3,1))
dim(lCPM)

plot(density(lCPM[,1]), col=cols[1], lwd=2, ylim=c(0,0.5), xlim=c(-3,15),las=2, main="", xlab="")
title(main="Density plot of log2(CPM)_all_Genes", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:6){
  den <- density(lCPM[,i])
  lines(den$x, den$y, col=cols[i], lwd=2)
}
legend("topright", colnames(lCPM), text.col=cols, bty="n",cex=0.8)

plot(density(lCPM.keep[,1]), col=cols[1], lwd=2, ylim=c(0,0.2), xlim=c(-7,15),las=2, main="", xlab="")
title(main="log2(CPM)  lcmp>0.5 & >=3 samples", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:6){
  den <- density(lCPM.keep[,i])
  lines(den$x, den$y, col=cols[i], lwd=2)
}
legend("topright", colnames(lCPM), text.col=cols, bty="n",cex=0.8)

par(mfrow=c(1,2),mar=c(6,5,3,1))
boxplot(lCPM.keep,col=cols,main="Boxplot of lCPM (counts.keep)",las=2,cex.axis=0.8)
abline(h=median(lCPM.keep),col="grey")
makeRLEplot(lCPM.keep,col=cols,main="RLE of logCPM (counts.keep)",ylim=c(-0.8,0.8),cex=0.8)

groups<-(c("MIM159", "MIM159", "MIM159",
                       "Wildtype", "Wildtype", "Wildtype"))
design<-model.matrix(~ 0 + groups)
table(groups)

design

y0 <- DGEList(counts.keep)
y0$samples

y <- calcNormFactors(y0,method="TMM")
y$samples

v <- voomWithQualityWeights(y,design,plot = TRUE,col=cols)

dim(v$E)
head(v$E)

out<-cbind(g.inf.keep,v$E)
head(out)

write.csv(out,file="tobacco_kept_gene_expression_2021.csv")

par(mfrow=c(1,2),mar=c(7,5,3,1))

boxplot(v$E, xlab="", ylab="v$E",las=2,main="Boxplot after TMM normalization",col=cols,cex.axis=0.8)
##  add a grey horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="grey")

makeRLEplot(v$E,main="RLE plot of v$E",ylim=c(-0.8,0.8),cols=cols,cex=0.6)

mycol=colorRampPalette(c("red","royalblue"))(10)
corrgram(v$E,upper.panel= panel.cor,labels=substr(colnames(v$E),start=1,stop=2),
         col.regions = colorRampPalette(c("red4","red", "salmon","white","red","purple", "royalblue"))
)

par(mfrow=c(1,1),mar=c(6,5,3,1))
plotMDS(v$E,col=cols,main="MDS of six samples")
legend("bottomright",legend = c("1:MIM159","2:WT"),text.col=c("red","purple"),cex=0.6)

v <- voomWithQualityWeights(y,design,plot = FALSE)
rownames(v)<-g.inf.keep$Symbol

fit <- lmFit(v)

cont.matrix <- makeContrasts(MIM159_vs_Wildtype=groupsMIM159-groupsWildtype,
                             levels=design)
cont.matrix

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont

fit.cont <- eBayes(fit.cont)
fit.cont

dim(fit.cont)

summa.fit <- decideTests(fit.cont,adjust.method = "BH")
summary(summa.fit)

options(digits = 2)

fit.cont$gene<-g.inf.keep[,c(1,2)]
for (i in 1) { 
  out<-topTable(fit.cont,coef=i,sort.by="p",genelist=fit.cont$gene,n=nrow(fit))
  print("+++++++++++++++++++++++++++++++++++")
  print(paste0("Number of DEGs for ",colnames(cont.matrix)[i]," is: " ))
  print(summary(summa.fit)[,i])
  print("top10 DEGs are:")
  print(topTable(fit.cont,coef=i,sort.by="p",genelist=fit.cont$gene))
  write.csv(as.matrix(out),file=paste0("topTable_",colnames(cont.matrix)[i],"_2021.csv"))
}

nice.col <- brewer.pal(6,name="Dark2")
for (i in c(1)) { 
  glXYPlot(x=fit.cont$coefficients[,i], y=-log10(fit.cont$p.value[,i]),
           xlab="logFC", ylab="-log10(p.value)", main=colnames(cont.matrix)[i],
           counts=y$counts, groups=groups, status=summa.fit[,i],
           #         anno=fit.cont$gene, side.main="Sym", folder="volcano_Tissue_vs_Bone_withqualityWeights")
           anno=fit.cont$gene, side.main="X", folder=paste0("volcano_2021",colnames(cont.matrix)[i]))
}
write.csv(fit.cont$coefficients, "logFC.csv")
for (i in c(1,2)){
  idx<-c(summa.fit[,i]=='-1' | summa.fit[,i]=='1')
  plot(x=fit.cont$coefficients[,i], y=-log10(fit.cont$p.value[,i]),pch=20,col="darkgrey",
       xlab="log2(FC)", ylab="-log10(adj.p.value)", main=colnames(cont.matrix)[i], cex=0.75)
  points(x=fit.cont$coefficients[summa.fit[,i]=='-1',i], y=-log10(fit.cont$p.value[summa.fit[,i]=='-1',i]),col="#8795E8",cex=0.75, pch=20)
  points(x=fit.cont$coefficients[summa.fit[,i]=='1',i], y=-log10(fit.cont$p.value[summa.fit[,i]=='1',i]),col="#FF6AD5",cex=0.75, pch=20)
  # label top 12 DEGs
  g.lab<-rownames(topTable(fit.cont,coef=i,sort.by="p",genelist=fit.cont$gene,n=12))
  idx.lab<-rownames(fit.cont)%in%g.lab
  text(x=fit.cont$coefficients[idx.lab,i], y=(-log10(fit.cont$p.value[idx.lab,i])+0.15),col="purple",cex=0.7,labels =fit.cont$gene$X[idx.lab])
}

fit.treat <- treat(fit.cont,lfc=log2(1.5))
summa.treat <- decideTests(fit.treat)
sum.treat<-summary(summa.treat)
sum.treat

options(digits = 2)
topTreat(fit.treat,coef=1,sort.by="p",genelist=fit.cont$gene)

for (i in 1){   
  out<-topTreat(fit.treat,coef=i,sort.by="p",n=nrow(fit),genelist=fit.treat$gene)
  print("+++++++++++++++++++++++++++++++++++")
  print(paste0("Number of treat DEGs for cmparison ",i," : ",colnames(cont.matrix)[i]," is: " ))
  print(sum.treat[,i])
  print("top10 treat DEGs are:")
  print(topTreat(fit.treat,coef=i,sort.by="p",genelist=fit.treat$gene))
  write.csv(as.matrix(out),file=paste0("topTreat_",colnames(cont.matrix)[i],"_2021.csv"))
}

for (i in 1){ 
  glMDPlot(fit.treat, coef=i, counts=y$counts, groups=groups,
           status=summa.treat[,i],anno=fit.cont$gene, side.main="X", main=colnames(cont.matrix)[i],
           folder=paste0("MD_",colnames(cont.matrix)[i]))
}


top70DE.sym<-topTreat(fit.treat,coef=1,sort.by="p",genelist=fit.treat$gene,n=70)$X
idx<-match(top70DE.sym,fit.treat$gene$X)
heatmap.2(x=as.matrix(v$E[idx,]),distfun=mydist,hclustfun=myclust,scale="row",labRow=fit.treat$gene$X[idx],srtCol=30,
          #        dendrogram="row",Colv=FALSE,
          ColSideColors=cols, col=bluered(32),cexRow=0.8,cexCol=1.5,trace="none",las=2,margins=c(10,8),lhei=c(0.4,2))

pdf("Rplot_heatmap_top100_comparison26_topTreat.pdf",height=10)
for (i in c(1)) { 
  top100DE.sym<-topTreat(fit.treat,coef=i,sort.by="p",genelist=fit.treat$gene,n=100)$X
  #  g.idx<-groups%in%unlist(strsplit(colnames(cont.matrix)[i],split="_vs_"))
  
  idx<-match(top100DE.sym,fit.treat$gene$X)
  heatmap.2(x=as.matrix(v$E[idx,]),distfun=mydist,hclustfun=myclust,scale="row",labRow=fit.cont$gene$X[idx],     srtCol=20,ColSideColors=cols,
            #        dendrogram="row",Colv=FALSE,
            main=colnames(cont.matrix)[i],
            col=bluered(32),cexRow=0.4,cexCol=0.8,trace="none",las=2,margins=c(8,10),lhei=c(0.4,2))
  
}
dev.off()

pdf("Rplot_heatmap_top100_comparisons26_topTable.pdf",height=9)
for (i in c(1)) { 
  top100DE.sym<-topTable(fit.cont,coef=i,sort.by="p",genelist=fit.cont$gene,n=100)$X
  # g.idx<-groups%in%unlist(strsplit(colnames(cont.matrix)[i],split="_vs_"))
  
  idx<-match(top100DE.sym,fit.cont$gene$X)
  heatmap.2(x=as.matrix(v$E[idx,]),distfun=mydist,hclustfun=myclust,scale="row",labRow=fit.cont$gene$X[idx],srtCol=30,
            #        dendrogram="row",Colv=FALSE,
            main=colnames(cont.matrix)[i],ColSideColors=cols,
            col=bluered(32),cexRow=0.4,cexCol=0.8,trace="none",las=2,margins=c(8,8),lhei=c(0.4,2))
}
dev.off()

head(summa.fit)

out<-cbind(fit.cont$gene,summa.fit)
rownames(out)<-NULL
head(out)

write.csv(as.matrix(out),file="decisionTest_allComparisonsTopTable.csv")

table(rowSums(summa.fit==0)==2)

delist<-out[rowSums(summa.fit==0)!=6,]
dim(delist)

out<-cbind(fit.cont$gene,summa.treat)
rownames(out)<-NULL
head(out)

write.csv(as.matrix(out),file="decisionTest_allComparisonsTopTreat.csv")


table(rowSums(summa.treat==0)==6)

delist<-summa.treat[rowSums(summa.treat==0)!=6,]
dim(delist)

par(mfrow=c(1,2),mar=c(7,3,3,1))
#g<-c(1:3)
g<-c(1)
a<-vennCounts(summa.fit[,g], include="both")
vennDiagram(a, include="both", names=NULL, mar=rep(1,4), cex=c(1,1.2,1.2), lwd=1,
            circle.col=c( "darkblue","red4","red","blue","darkgreen"),main="topTable", counts.col=c("blue","red"))


vennDiagram(summa.fit[,g], include=c("up","down"), names=NULL, mar=rep(1,2), cex=c(1,1.2,1.2), lwd=1,main="topTable",
            circle.col=c( "darkblue","red4","red","blue","darkgreen"), counts.col=c("blue","red"))
