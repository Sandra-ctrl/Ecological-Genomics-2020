library(methylKit)
library(tidyverse)
library(ggplot2)
library(pheatmap)

getwd() 
# first, we want to read in the raw methylation calls with methylkit

# set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)
dir <- "/Users/sandr/OneDrive/Documents/GitHub/Ecological-Genomics-2020/Epigenetics_data"

# read in the sample ids
samples <- read.table("~/Github/Ecological-Genomics-2020/Epigenetics_data/sample_id.txt", header=FALSE)

# now point to coverage files
files <- file.path(dir, samples$V1)
all(file.exists(files))

# convert to list
file.list <- as.list(files)

# get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))

# use methRead to read in the coverage files
myobj <- methRead(location = file.list,
                  sample.id =   nmlist,
                  assembly = "atonsa",
                  dbtype = "tabix",
                  context = "CpG",
                  resolution = "base",
                  mincov = 20,
                  treatment = 
                    c(0,0,0,0,
                      1,1,1,1,
                      2,2,2,2,
                      3,3,3,3,
                      4,4,4,4),
                  pipeline = "bismarkCoverage",
                  dbdir = "~/Documents/GitHub/Ecological-Genomics-2020/Epigenetics_data")
getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
par(mfrow=c(1,1))

par(mfrow=c(5,4))
for(i in 1:length(myobj)) 
  getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)

meth <- methylKit:::readMethylBaseDB(
  dbpath = "/Users/sandr/OneDrive/Documents/GitHub/Ecological-Genomics-2020/Epigenetics_data/methylBase_united.txt.bgz",
  dbtype = "tabix",
  sample.id =   unlist(nmlist),
  assembly = "atonsa", 
  context = "CpG",
  resolution = "base",
  treatment = c(0,0,0,0,
                1,1,1,1,
                2,2,2,2,
                3,3,3,3,
                4,4,4,4),
  destrand = FALSE)
meth
pm <- percMethylation(meth)
ggplot(gather(as.data.frame(pm)), aes(value)) + 
  geom_histogram(bins = 10, color="black", fill="grey") + 
  facet_wrap(~key)

sp.means <- colMeans(pm)
dim(pm)
p.df <- data.frame(sample=names(sp.means),
                   group = substr(names(sp.means), 1,6),
                   methylation = sp.means)

ggplot(p.df, aes(x=group, y=methylation, color=group)) + 
  stat_summary(color="black") + geom_jitter(width=0.1, size=3) 

clusterSamples(meth, dist="correlation", method="ward.D", plot=TRUE)

PCASamples(meth, screeplot=TRUE)
PCASamples(meth, screeplot=FALSE)

# subset with reorganize()

meth_sub <- reorganize(meth,  sample.ids= (c("AA_F25_1","AA_F25_2","AA_F25_3", "AA_F25_4",
                                             "HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4")), 
                       treatment=c(0,0,0,0,1,1,1,1),
                       save.db=FALSE)

# calculate differential methylation

myDiff=calculateDiffMeth(meth_sub,
                         overdispersion="MN",
                         mc.cores=1,
                         suffix = "AA_HH", adjust="qvalue",test="Chisq")
myDiff

# where MN corrects for overdispersion
# fit a logistic regression to methylation values where explanatory variable is the treatment (case or control). 
# and we compare the fit of the model with explanatory variable vs the null model (no explanatory variable) 
# and ask if the fit is better using a Chisq test. 
# the methylation proportions are weighted by their coverage, as in a typical logistic regression. Note that in theory you could enter these as two column success and failure data frame, which is common in logistic regressions.

# use overdispersion: Chisq without overdispersion finds more true positives, but many more false positives. good compromise is overdispersion with Chisq. reduced true pos, but really reduces false pos rate.

# get all differentially methylated bases
myDiff <-getMethylDiff(myDiff,difference=10,qvalue=0.05)
myDiff

# we can visualize the changes in methylation frequencies quickly.
hist(getData(myDiff)$meth.diff)

# get hyper methylated bases
hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper")
#
# get hypo methylated bases
hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo")

#heatmaps first

# get percent methylation matrix
pm <- percMethylation(meth_sub)

# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]

# add snp, chr, start, stop

din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, getData(myDiff)$start, sep=":"), din, pm.sig)
colnames(df.out) <- c("snp", colnames(din), colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))

####
# heatmap
####

my_heatmap <- pheatmap(pm.sig,
                       show_rownames = FALSE)

ctrmean <- rowMeans(pm.sig[,1:4])

h.norm <- (pm.sig-ctrmean)

my_heatmap <- pheatmap(h.norm,
                       show_rownames = FALSE)

##### if you want to change colors. only because I don't love the default colors.

paletteLength <- 50
myColor <- colorRampPalette(c("cyan1", "black", "yellow1"))(paletteLength)
myBreaks <- c(seq(min(h.norm), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(h.norm)/paletteLength, max(h.norm), length.out=floor(paletteLength/2)))

my_heatmap <- pheatmap(h.norm,
                       color=myColor, 
                       breaks=myBreaks,
                       show_rownames = FALSE)

#####
#let's look at methylation of specific gene or snp
####

df.out
df.plot <- df.out[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")
df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)

# looking at snp LS049205.1:248
# if you choose a different snp, you can create different plots.

df.plot %>% filter(snp=="LS042748.1:7742") %>% 
  ggplot(., aes(x=group, y=methylation, color=group, fill=group)) +
  stat_summary(fun.data = "mean_se", size = 2) +
  geom_jitter(width = 0.1, size=3, pch=21, color="black")

## write bed file for intersection with genome annotation

write.table(file = "/Users/sandr/OneDrive/Documents/GitHub/Ecological-Genomics-2020/Epigenetics_data/diffmeth.bed",
            data.frame(chr= df.out$chr, start = df.out$start, end = df.out$end),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
