#
#
#
# Script to calculate the differential expression of pair experimental designs.
#
#
# Tested in R 2.15
#
#
#

library(edgeR)
library(knitr)
#library(R2HTML)
#library(gplots)



# Obtaining the paths and the groups information.

cmd <- commandArgs(TRUE)

raw_expression <- cmd[1]
raw_expression
out <- cmd[2]
out
samples <- make.names(unlist(strsplit(cmd[3],";")))
samples
groups <- unlist(strsplit(cmd[4],";"))
groups
pairs <- unlist(strsplit(cmd[5],";"))
pairs
expname <- cmd[6]
expname

# Creating the table of raw expression.

expression  <- read.table( raw_expression , header=FALSE) 
row.names(expression) <- (expression$V1)
expression$V1 <- NULL
colnames(expression) <- samples 

#knit('differential_ex_nonpair.R')


# Creation of the edgeR object
expression.dglist <- DGEList(expression, group=groups)

# Only counting miRNA if there are more than 5 counts in at least half of the samples.
expression.dglist <- expression.dglist[rowSums(cpm(expression.dglist)> 5) > (length(expression)-2)/2,]

# Calculation of the normalization factors which correct for the different compositions of the samples. 
expression.dglist <- calcNormFactors(expression.dglist)

#pdf(paste(out,"/",expname,"_gra.pdf",sep=""), width=11,height=8.5)

# Plotting Quality Control of Differential Expression
png(paste(out,"/QUA_NOR_plot1.png",sep=""),width=1280,height=960,res=150)
par( mfrow=c(1,1), mar=c(1,1,1,1), oma=c(3,1,0,0) )
boxplot(log2(expression[rowSums(expression[length(expression)]) != 0,1:length(expression)]),las=2,xlab="",ylab="Raw log2 Counts",main="Raw Expression Counts")
normcnt <- cpm(expression.dglist, normalized.lib.sizes=T)
dev.off()

# display normalized, filtered counts
png(paste(out,"/QUA_NOR_plot2.png",sep=""),width=1280,height=960,res=150)
boxplot(log2(normcnt),xlab="",ylab="Normalized log2 Counts",main="Normalized Counts: \n miRNAs with > 5 counts in more than half of the samples",las=2)
dev.off()
par( mfrow=c(1,1), mar=c(4,5,1,1), oma=c(3,1,2,1))

#########################################

# ploting Multidimensional Scaling
png(paste(out,"/MDS_plot.png",sep=""),width=1280,height=960,res=150)
plotMDS(expression.dglist, main=paste("MDS plot"))
dev.off()

# Design matrix for paired groups
design <- model.matrix(~pairs+expression.dglist$samples$group)

rownames(design) <- rownames(expression.dglist$samples)
colnames(design)[length(colnames(design))] <- expname

## Generalized Linear Model (GLM) is the statistical method used in this pipeline.

# Estimation of the Common Dispersion. Each miRNA gets assigned the same dispersion estimate. 
expression.dglist <- estimateGLMCommonDisp(expression.dglist, design)

# Estimation of the Tagwise. Each miRNA will get its own unique dispersion estimate.
expression.dglist.tagwdis <- estimateGLMTagwiseDisp(expression.dglist, design)

# Summary
names(expression.dglist.tagwdis)
head(expression.dglist.tagwdis$tagwise.dispersion)
summary(expression.dglist.tagwdis$tagwise.dispersion)



# Using GLM model looking for differences between expressed miRNAs.
glm.expression.dglist.tagdis <-glmFit(expression.dglist.tagwdis, design, dispersion=expression.dglist.tagwdis$tagwise.dispersion)
lrt.expression.dglist.tagdis <- glmLRT(glm.expression.dglist.tagdis)

# Normalized expression.
normcnt <- round(cpm(expression.dglist.tagwdis, normalized.lib.sizes=T))
nmexpr.eRout.tagdis <- merge(normcnt, topTags(lrt.expression.dglist.tagdis, n=10000), by.x="row.names", by.y="row.names")
names(nmexpr.eRout.tagdis)[1] <- "miRNA"
nmexpr.eRout.tagdis <- nmexpr.eRout.tagdis[order(nmexpr.eRout.tagdis$PValue),]
write.table(nmexpr.eRout.tagdis, paste(out,"/differential_expression.xls",sep=""),sep="\t", quote=F,row.names=F)

## Plotting Results

# Plotting the relationship between concentration and fold-change across the genes. 
pvalue <- nmexpr.eRout.tagdis[nmexpr.eRout.tagdis[,"PValue"] < 0.05,]
de.genes.tgw <- rownames(pvalue)

png(paste(out,"/SMEAR_plot.png",sep=""),width=1280,height=960,res=150)
plotSmear(expression.dglist, de.tags=de.genes.tgw , main="Poisson", cex= .5, xlab="Log Concentration", ylab="log Fold-Change")
dev.off()

#  Plotting FDR Volcano Plot
png(paste(out,"/FDR.volcano_plot.png",sep=""),width=1280,height=960,res=150)
plot(nmexpr.eRout.tagdis$logFC, -log10(nmexpr.eRout.tagdis$FDR), col=ifelse(nmexpr.eRout.tagdis$FDR<0.05,"red","black"), main="FDR volcano plot", xlab="log2FC", ylab="-log10(FDR)")
dev.off()

#  Plotting P-Value Volcano Plot
png(paste(out,"/PValue_plot.png",sep=""),width=1280,height=960,res=150)
plot(nmexpr.eRout.tagdis$logFC, -log10(nmexpr.eRout.tagdis$PValue),col=ifelse(nmexpr.eRout.tagdis$PValue<0.01,"red","black"),main="P-value volcano plot",xlab="log2FC",ylab="-log10(Pvalue)")
dev.off()


#  Plotting Histogram of P-Value distribution.
png(paste(out,"/PValue_dis_plot.png",sep=""),width=1280,height=960,res=150)
hist(nmexpr.eRout.tagdis$PValue, breaks=20,xlab="P Value", ylab="Frequency",main="p-value distribution")
dev.off()

#  Plotting barplot with miRNA PValue > 0.05 and logFC > |1|.
Change <- pvalue[pvalue[,"logFC"] > 1 | pvalue[,"logFC"] < -1 | pvalue[,"FDR"] < 0.01,]
Change.sort <- Change[order(Change$logFC),]


png(paste(out,"/barplot_plot.png",sep=""),width=1280,height=960,res=150)
barplot(Change.sort $logFC, names.arg=Change.sort$miRNA, cex.names= 0.6, las=2, col=rainbow(100), axes=TRUE, ylab='logFC',  main='miRNA that |logFC| > 1 , P-value < 0.05 , FDR < 0.01' )
dev.off()


#  Writing table with miRNA PValue > 0.05 and logFC > |1|
write.table(Change.sort , paste(out,"/dif_exp_sum.xls",sep=""),sep="\t", quote=F,row.names=F)

lon <- length(samples) + 1

heat <-  Change.sort[,c(2:lon)] 
rownames(heat) <- Change.sort[,1]

png(paste(out,"/heatmap_plot.png",sep=""),width=800,res=150)
heatmap(as.matrix(heat))
dev.off()



#heatmap.2(as.matrix(heat), dendrogram="col", cellnote=as.matrix(PVal),
#          notecol="black",col=redgreen(75),scale="none",key=TRUE, keysize=1.5,
#          density.info="none", trace="none", cexRow=0.7,cexCol=1.2)

#############################################################################################
#
#
# Comparison. If you have the real data in a table with a format similar to the table .dif_exp_sum.xls you can look for correlation.  data2 <- Data to compare.
#
#

#data1 <- Change.sort[,c("miRNA","logFC")]
#rownames(data1) <- data1[,"miRNA"]
#data1[,"miRNA"] <- NULL

#data2 <- read.table('../cap.csv',header=FALSE)
#rownames(data2) <- data2[,"V1"]
#data2[,"V1"] <- NULL

#dflist <- list(data1,data2)
#idx <-Reduce(intersect,lapply(dflist, rownames))
#idx.graph <- data.frame(cbind(data1[idx,],data2[idx,]))

#plot(idx.graph$X1, xlab='miRNA' , ylab='logFC', col="red", main="miRNA correlation", sub="Black points= CAP-miRSEQ Pipeline,             Red points= Onyx Pipeline")
#points(idx.graph$X2, )

#print('miRNA found for onyx but not present in the other database')
#data2.only <- setdiff(rownames(data2),idx)
#data2.only

#print('miRNA found in the other data but not present in onyx')
#data1.only <- setdiff(rownames(data1),idx)
#data1.only

#############################################################################################


#dev.off()

#Print HTML 

#html = HTMLInitFile( out , filename=expname, CSSFile="style.css",  Title="Onyx | miRNA Differential Expression Analysis Pipeline")

#HTML("<head>", file=html)
#HTML("<title>Onyx | miRNA Differential Expression Analysis Pipeline</title>", file=html)
#HTML("</head>", file=html)
#HTML("<body>", file=html)
#HTML(Change.sort, file=html)
#HTML("</body>", file=html)

#HTMLEndFile()
