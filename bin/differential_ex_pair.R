library(edgeR)

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


expression  <- read.table( raw_expression , header=FALSE) 
row.names(expression) <- (expression$V1)
expression$V1 <- NULL
colnames(expression) <- samples 



expression.dglist <- DGEList(expression, group=groups)
expression.dglist <- expression.dglist[rowSums(cpm(expression.dglist)> 5) > (length(expression)-2)/2,]
expression.dglist <- calcNormFactors(expression.dglist)

pdf(paste(out,"/",expname,"_gra.pdf",sep=""), width=11,height=8.5)
plotMDS(expression.dglist, main=paste("MDS plot"))

design <- model.matrix(~pairs+expression.dglist$samples$group)

rownames(design) <- rownames(expression.dglist$samples)
colnames(design)[length(colnames(design))] <- expname

expression.dglist <- estimateGLMCommonDisp(expression.dglist, design)
expression.dglist.tagwdis <- estimateGLMTagwiseDisp(expression.dglist, design)
names(expression.dglist.tagwdis)
head(expression.dglist.tagwdis$tagwise.dispersion)
summary(expression.dglist.tagwdis$tagwise.dispersion)

glm.expression.dglist.tagdis <-glmFit(expression.dglist.tagwdis, design, dispersion=expression.dglist.tagwdis$tagwise.dispersion)
lrt.expression.dglist.tagdis <- glmLRT(glm.expression.dglist.tagdis)

normcnt <- round(cpm(expression.dglist.tagwdis, normalized.lib.sizes=T))
nmexpr.eRout.tagdis <- merge(normcnt, topTags(lrt.expression.dglist.tagdis, n=10000), by.x="row.names", by.y="row.names")
names(nmexpr.eRout.tagdis)[1] <- "miRNA"
nmexpr.eRout.tagdis <- nmexpr.eRout.tagdis[order(nmexpr.eRout.tagdis$PValue),]
write.table(nmexpr.eRout.tagdis, paste(out,"/",expname,".differential_expression.xls",sep=""),sep="\t", quote=F,row.names=F)

plot(nmexpr.eRout.tagdis$logFC, -log10(nmexpr.eRout.tagdis$FDR), col=ifelse(nmexpr.eRout.tagdis$FDR<0.05,"red","black"), main="FDR volcano plot", xlab="log2FC", ylab="-log10(Pvalue)")

hist(nmexpr.eRout.tagdis$PValue, breaks=20,xlab="P Value", ylab="Frequency",main="p-value distribution")


pvalue <- nmexpr.eRout.tagdis[nmexpr.eRout.tagdis[,9] < 0.05,]
Change <- pvalue[pvalue[,6] > 1 | pvalue[,6] < -1,]
Change.sort <- Change[order(Change$logFC),]

barplot(Change.sort $logFC, names.arg=Change.sort$miRNA, cex.names= 0.6, las=2, col=rainbow(100), axes=TRUE, ylab='logFC', xlab='miRNA', main='miRNA that |logFC| > 1 and P-value < 0.05' )

write.table(Change.sort , paste(out,"/",expname,".dif_exp_sum.xls",sep=""),sep="\t", quote=F,row.names=F)


# Comp


data1 <- Change.sort[,c("miRNA","logFC")]
rownames(data1) <- data1[,"miRNA"]
data1[,"miRNA"] <- NULL

data2 <- read.table('../cap.csv',header=FALSE)
rownames(data2) <- data2[,"V1"]
data2[,"V1"] <- NULL

dflist <- list(data1,data2)
idx <-Reduce(intersect,lapply(dflist, rownames))
idx.graph <- data.frame(cbind(data1[idx,],data2[idx,]))

plot(idx.graph$X1, xlab='miRNA' , ylab='logFC', col="red", main="miRNA correlation", sub="Black points= CAP-miRSEQ Pipeline,             Red points= Onyx Pipeline")
points(idx.graph$X2, )

data2.only <- setdiff(rownames(data2),idx)
data2.only
data1.only <- setdiff(rownames(data1),idx)
data1.only
dev.off()
