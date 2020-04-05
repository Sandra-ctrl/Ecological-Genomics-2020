library(methylKit)
library(tidyverse)
library(ggplot2)
library(pheatmap)

library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE)
detach("ggplot2")
detach("package:ggplot2", unload=T)
# first, we want to read in the raw methylation calls with methylkit

# set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)

# read in the sample ids
samples <- read.table("~/Github/Ecological-Genomics-2020/sample_id.txt", header=FALSE)

# now point to coverage files
files <- file.path(dir, samples$V1)
getwd()
dir <- "/Users/sandr/OneDrive/Documents/GitHub/Ecological-Genomics-2020"
files <- file.path(dir, samples$V1)
all(file.exists(files))
# convert to list
file.list <- as.list(files)
# get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))
# use methRead to read in the coverage files
myobj <- methRead()
myobj <- methRead(location= file.list,
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
                  dbdir = "~/Documents/GitHub/Ecological-Genomics-2020")

######
# visualize coverage and filter
######

# We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats()

# and can plot all of our samples at once to compare.

# filter samples by depth with filterByCoverage()
filtered.myobj <- filterByCoverage()

######
# merge samples
######

#Note! This takes a while and we're skipping it

# use unite() to merge all the samples. We will require sites to be present in each sample or else will drop it

meth <- unite()