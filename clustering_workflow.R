# clustering_exploration.R workflow with example cluster/exploratory data analyses
# RNA-Seq Workshop, University of Zurich, February 2017
# Laura Carroll, lmc297@cornell.edu

# Step 0: set our working directory to our bayseq output directory
#setwd("/path/to/normalized_counts.txt/directory")
library(vegan)
library(ggplot2)

# Step 1: load in our data (normalized read counts from bayseq_workflow.R)
infile <- read.delim("normalized_counts.txt")
# let's get columns 2 through 9 of this file
countfile <- infile[,2:9]
rownames(countfile) <- infile[,1]
# how many differentially-expressed genes do we have?
max(which((infile$FDR.DE<0.05 & infile$FC >= 2.5)|(infile$FDR.DE<0.05 & infile$FC <= 0.4)))
# let's make a heatmap with our differentially-expressed genes
top329 <- as.matrix(countfile[1:329,])
heatmap(top329,distfun = function(x) vegdist(x,method = "bray"))


# Step 2: K-Means
# let's pick a maximum number of clusters to test
maxclus=10
# loop through each value of k, saving the total within-cluster sum of squares
kvec<-c()
for(k in 1:maxclus){
  kout<-kmeans(countfile,centers=k,iter.max = 100, nstart=10)
  testk<-kout$tot.withinss
  kvec<-append(kvec,testk)
}
# let's plot our k on the 
plot(x=c(1:maxclus),y=kvec,xlab = "K", ylab="Total Within SS")
# Where is our "elbow?" Use the value you want for "centers=" below
kout<-kmeans(countfile,centers=2,iter.max = 100, nstart=10)
library(cluster)
# let's plot our k-means output
clusplot(countfile, kout$cluster,color = TRUE,shade = TRUE,labels=2,lines=0)

# Step 3: PCA
# let's run pca on our genes
pca<-prcomp(countfile,scale = TRUE,center = TRUE)
summary(pca)
# let's color our data by our k-means cluster designations
covar<-as.character(kout$cluster)
names(covar)<-names(kout$cluster)
data <- as.data.frame(cbind(covar,pca$x[,1],pca$x[,2],pca$x[,3]))

# Make base plot
g <- ggplot(data, aes(x=pca$x[, 1], y=pca$x[, 2], size=pca$x[, 3]))

# Plot colored by k-means cluster
g + geom_point(aes(color=covar)) +
  theme_bw() +
  scale_color_manual(values = c("blue",
                                "magenta"))

# now let's run pca on our samples!
# first, transpose the matrix
tk<-t(countfile)
pca<-prcomp(tk,scale = TRUE,center = TRUE)
summary(pca)
# let's color it by condition (bhi or css)
covar<-ifelse(grepl("bhi",rownames(tk)),"bhi","css")
data <- as.data.frame(cbind(covar,pca$x[,1],pca$x[,2],pca$x[,3]))

# Make base plot
g <- ggplot(data, aes(x=pca$x[, 1], y=pca$x[, 2], size=pca$x[, 3]))

# Plot colored by condition
g + geom_point(aes(color=covar)) +
  theme_bw() +
  scale_color_manual(values = c("purple",
                                "green"))
