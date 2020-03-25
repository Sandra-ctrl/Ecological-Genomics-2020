countsTable <- read.table("RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable)
head(countsTableRound)
conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)
grep("10", names(countsTableRound), value = TRUE)
day10countstable <- subset(countsTableRound, grep("10", names(countsTableRound), value = TRUE)) #doesn't work has to be logical

day10countstable <- countsTableRound %>% select(contains("10"))
dim(day10countstable)

conds10<- subset(conds, day=="10")
dim(conds10)
head(conds10)

colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)

rowSums(countsTableRound)
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))

apply(countsTableRound,2,mean)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ climate + day + treatment)
dim(dds)
dds <- dds[rowSums(counts(dds)) > 76]
dim(dds)

dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds, alpha = 0.05)
res <- res[order(res$padj),]
head(res)

summary(res)

res_treatCD <- results(dds, name="treatment_D_vs_C", alpha=0.05)
res_treatCD <- res_treatCD[order(res_treatCD$padj),]
head(res_treatCD)
summary(res_treatCD)

plotMA(res_treatCD,ylim=c(-3,3))

vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd,intgroup=c("climate","treatment","day"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))
data$day <- factor(data$day, levels=c("0","5","10"), labels = c("0","5","10"))

ggplot(data, aes(PC1, PC2, color=day, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

d <-plotCounts(dds, gene="MA_7017g0010", intgroup = (c("treatment","climate","day")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=day, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p

p <-ggplot(d, aes(x=treatment, y=count, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p

library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate","day")])
pheatmap(mat, annotation_col=df)



