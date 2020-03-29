library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

## Import the counts matrix
countsTable <- read.table("RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)

## Import the samples description table - links each sample to factors of the experimental design.
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)

############ Try with only Day 10 data

# grep("10", names(countsTableRound), value = TRUE)
# day10countstable <- subset(countsTableRound, grep("10", names(countsTableRound), value = TRUE)) #doesn't work has to be logical

day10countstable <- countsTableRound %>% select(contains("10"))
dim(day10countstable)

conds10<- subset(conds, day=="10")
dim(conds10)
head(conds10)

## Let's see how many reads we have from each sample:
colSums(day10countstable)
mean(colSums(day10countstable))
barplot(colSums(day10countstable), las=3, cex.names=0.5,names.arg = substring(colnames(day10countstable),1,13))
abline(h=mean(colSums(day10countstable)), col="blue", lwd =2)

# What's the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))
# wow! This shows dispersion across genes - differences in magnitude of expression

# What's the average number of counts per gene
rowSums(day10countstable)
mean(rowSums(day10countstable))
median(rowSums(day10countstable))

# What's the average number of counts per gene per sample
apply(countsTableRound,2,mean)
apply(day10countstable,2,mean)

## Create a DESeq object and define the experimental design here with the tilde

dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)
dim(dds)

# Filter out genes with few reads 

dds <- dds[rowSums(counts(dds)) > 30]
dim(dds)
# 24300    30
## Run the DESeq model to test for differential gene expression: 
# 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 
# 3) run negative binomial glm
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# Running the model: design = ~ climate + treatment + climate:treatment

# [1] "Intercept"            "climate_HD_vs_CW"     "treatment_D_vs_C"    
# [4] "treatment_H_vs_C"     "climateHD.treatmentD" "climateHD.treatmentH"

# Order and list and summarize results from specific contrasts
# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
res <- results(dds, name="treatment_D_vs_C", alpha = 0.05)
res <- res[order(res$padj),]
head(res)
summary(res)

res <- results(dds, name="treatment_H_vs_C", alpha = 0.05)
res <- res[order(res$padj),]
head(res)
summary(res)

res <- results(dds, alpha = 0.05)
res <- res[order(res$padj),]
head(res)

dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)
summary(dds)
dds <- dds[rowSums(counts(dds))>76]
dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)
summary(dds)
dim(dds)
dds <- dds[rowSums(counts(dds))>30]
dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)
summary(dds)
dim(dds)
dds <- dds[rowSums(counts(dds)) > 30]
dim(dds)
dds <- DESeq(dds)
resultsNames(dds)

res_treatCD <- results(dds, name="treatment_D_vs_C", alpha=0.05)
res_treatCD <- res_treatCD[order(res_treatCD$padj),]
head(res_treatCD)
summary(res_treatCD)

plotMA(res_treatCD,ylim=c(-3,3))

res_treatCH <- results(dds, name="treatment_H_vs_C", alpha = 0.05)
res_treatCH <- res_treatCH[order(res_treatCH$padj),]
head(res_treatCH)
summary(res_treatCH)

plotMA(res_treatCH,ylim=c(-3,3))

res_interClimTreat <- results(dds, name="climateHD.treatmentD", alpha=0.05)
res_interClimTreat <- res_interClimTreat[order(res_interClimTreat$padj),]
head(res_interClimTreat)
summary(res_interClimTreat)
plotMA(res_interClimTreat,ylim=c(-3,3))

res_interClimTreatH <- results(dds, name="climateHD.treatmentH", alpha=0.05)
res_interClimTreatH <- res_interClimTreatH[order(res_interClimTreatH$padj),]
head(res_interClimTreatH)
summary(res_interClimTreatH)
plotMA(res_interClimTreatH,ylim=c(-3,3))

res_ClimTreat <- results(dds, name="climate_HD_vs_CW", alpha=0.05)
res_ClimTreat <- res_ClimTreat[order(res_ClimTreat$padj),]
head(res_ClimTreat)
summary(res_ClimTreat)
plotMA(res_ClimTreat,ylim=c(-3,3))

res_Int <- results(dds, name="intercept", alpha = 0.05)
res_Int <- res_Int[order(res_Int$padj),]
head(res_Int)
summary(res_Int)

# PCA
vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd,intgroup=c("climate","treatment"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("Control","Hot","Dry+Hot"))
data$climate <- factor(data$climate, levels=c("HD","CW"), labels = c("Hot-Dry","Cold-Wet"))


ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Counts of specific top gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="MA_133272g0010", intgroup = (c("treatment","climate")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p

p <-ggplot(d, aes(x=treatment, y=count, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate")])
pheatmap(mat, annotation_col=df)

