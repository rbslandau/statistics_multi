###################################################################################
# Script for Session 7 of the course: Applied multivariate Statistics with R  	  #
#						by Ralf B. Schäfer,	WS 2018/19							  #
# 						     Cluster analysis									  #
###################################################################################

# we first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine here
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

library(vegan)
werra_sp <- read.table("https://raw.githubusercontent.com/rbslandau/statistics_multi/master/Data/River_sp.csv", sep = ';', header = TRUE)
werra_env <- read.table("https://raw.githubusercontent.com/rbslandau/statistics_multi/master/Data/River_env.csv", sep = ';', header=TRUE)
# we load some species data from governmental stream monitoring and 
# related grouping information

###########################
# Data preparation        #
# Reduce weight of        #
# dominant taxa			  #
###########################	

### Transformation required?
# check species maxima
apply(werra_sp[ ,-1], 2, max)
# remove species that are absent in all sites
# see script for RDA
# for rationale of following lines
werra_pa <- decostand(werra_sp[ ,-1], "pa")
# calculate sum per species
werra_sum <- apply(werra_pa, 2, sum)
sort(werra_sum)
# remove species that occur in less than 2 sites
werra_fin <- werra_sp[ ,-1][ , ! werra_sum <3]

range(apply(werra_fin, 2, max))
# maxima differ strongly
range(apply(werra_fin^0.5, 2, max))
# also after square root transformation
range(apply(werra_fin^0.25, 2, max)) 
# double square root leads to stronger downweighing of abundant taxa
# is required in this case
werra_sp_t <- werra_fin^0.25

### Create Bray-Curtis distance matrix
werra_dist <- vegdist(werra_sp_t)
# if no distance measure is specified, vegdist uses the Bray Curtis distance
# Be careful, the dist function uses Euclidean as default!
# Note that the distance measure influences the cluster result

werra_clus <- hclust(werra_dist, method="average")
# cluster analysis with average linkage - 
# see lecture and R help for other methods available
# we also use Wards method for comparison
werra_clus2 <- hclust(werra_dist, method="ward.D2")

############################################
# Check preservation of initial distances  #
############################################

# we can check the relationship between the cophenetic matrix
# and the initial distance matrix to see how the distances 
# between the initial objects are preserved

# distance matrix before clustering
werra_dist
# distance matrix after clustering for average linkage
cophenetic(werra_clus)
# correlation of matrices for average linkage
cor(cophenetic(werra_clus), werra_dist)
# correlation for wards method
cor(cophenetic(werra_clus2), werra_dist)
# correlation is lower
# Compute Stress1 (analogous to NMDS) for average linkage
sqrt(sum((werra_dist - cophenetic(werra_clus))^2) / sum(werra_dist^2) ) 
# Compute Stress1 (analogous to NMDS) for Wards method
sqrt(sum((werra_dist - cophenetic(werra_clus2))^2) / sum(werra_dist^2) ) 
# Stress 1 much higher for Wards method

#############################
# Plot cluster dendrograms  #
#############################
# plot results for both methods
par(mfrow= c(1,2), cex = 1.5)
plot(werra_clus, ann = TRUE)
plot(werra_clus2, ann = TRUE)
# group 2-4 is placed in different cluster
# for Wards method

# check changes more elegantly
library(dendextend) 
# converting to dendogram objects as dendextend works with dendogram objects 
werr_clus_dend <- as.dendrogram(werra_clus)
werr_clus_dend2 <- as.dendrogram(werra_clus2)
par(mfrow= c(1,1), cex = 1.5)
tanglegram(werr_clus_dend, werr_clus_dend2, sort = TRUE, margin_outer= 3.5, main_left = "Average linkage", main_right = "Wards method", edge.lwd = 2)
# dashed lines indicate different clustering
# coloured lines connect same objects

###################################
# Using cluster validity indices  #
###################################
# how many groups are optimal?

# internal validation
# we have no prior information on group/cluster membership

# we compare the results for 2 to 4 groups
# for the average linkage method

# first we visualise the groups
par(mfrow= c(1,1), cex=1.5)
plot(werra_clus, ann=TRUE)
rect.hclust(werra_clus, k = 2, border = 'steelblue')
rect.hclust(werra_clus, k = 3, border = 'darkgreen')
rect.hclust(werra_clus, k = 4, border = 'orange')

# and extract the different grouping vectors
cut_2 <- cutree(werra_clus, k = 2)
cut_3 <- cutree(werra_clus, k = 3)
cut_4 <- cutree(werra_clus, k = 4)

# The plclust function  can be employed 
# to assign the cluster group to the dendrogram. 
# Check help of the function for further options
plot(werra_clus, labels = cut_2)
plot(werra_clus, labels = cut_4)

# the package factoextra also allows several options for nice dendrograms
library("factoextra")
fviz_dend(werra_clus)
# plot clusters for k = 4
fviz_dend(werra_clus, k = 4, # Cut in four groups
		k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
		color_labels_by_k = TRUE, # color labels by groups
		rect = TRUE, # Add rectangle around groups
		rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
		rect_fill = TRUE, lower_rect = -0.05
		)
# see Kassambara A. (2017) Practical guide to cluster analysis in R:
# unsupervised machine learning, Edition 1. STHDA.
# Page 84ff for more information on plot editing and different plot types
# such as phylogenetic plots	


# we can also visualise cluster results using nmds
werra_mds <- metaMDS(werra_dist)
# note the limitations of NMDS regarding
# precise representation of distances
# 
# Plot NMDS for cluster solution
par(cex=1.5)
ordiplot(werra_mds, type = "n")
points(werra_mds, pch = 16)
# next we add our cluster solution
ordispider(werra_mds, factor(cut_2), label = TRUE)
ordihull(werra_mds, factor(cut_2), lty = "dotted")

# check for 4 groups
ordiplot(werra_mds, type = "n")
points(werra_mds, pch = 16)
# next we add our cluster solution
ordispider(werra_mds, factor(cut_4), label = TRUE)
ordihull(werra_mds, factor(cut_4), lty = "dotted")
# 2 groups looks more convincing
# beside visual checking, several indices have been developed
# to evaluate the quality of the cluster solution 

##### computation of internal validation indices ######

# now we calculate several indices for the different groups
library(fpc)
indi_2 <- cluster.stats(werra_dist, cut_2)
indi_3 <- cluster.stats(werra_dist, cut_3)
indi_4 <- cluster.stats(werra_dist, cut_4)
# the function calculates all indices automatically

# note that the R package clusterSim
# allows for the comparison of several
# choices for distance measures and cluster methods
# NbClust is another R package that provides many 
# cluster validity indices. It also features automatic
# computation of indices for all clusters
# from minimum to maximum number of clusters defined

indi_2
# range of indices provided
# we restrict the comparison to a few indices

# Calinski and Harabasz index
indi_2$ch
indi_3$ch
indi_4$ch
# the higher the better
# the highest CH for 2 groups

# average silhouette width
indi_2$avg.silwidth
indi_3$avg.silwidth
indi_4$avg.silwidth
# highest avg. silhouette width for 2 groups
# avg. sildwidth closer to 0.2 (no cluster structure)
# than to 0.5 (reasonable cluster structure)

# calculate individual silhouette values s(x)
library(cluster)
sil <- silhouette(cut_2, werra_dist)
print(sil)
# shows that observation 8 seems misclassified
plot(sil)
# plot for individual s(x)

# GAP index
# for hierarchical clustering, we need to prepare 
# a helper function that directly returns the cluster vector
hclusCut <- function(x, k, d.meas = "bray", clus.meth="average" )
   list(cluster = cutree(hclust(vegdist(x, method=d.meas), method = clus.meth), k=k))
# GAP function
clusGap(werra_sp_t, hclusCut, K.max = 5, spaceH0="original")
# the SSQ are actually only for four clusters lower than from a 
# reference distribution

# for visualisation of the indices, see also Kassambara 2017: 
# chapters 12 to 14

##### computation of stability index ######
boot_hc1 <- clusterboot(werra_dist, distances=TRUE, clustermethod=disthclustCBI, method = "average", k = 2)
boot_hc2 <- clusterboot(werra_dist, distances=TRUE, clustermethod=disthclustCBI, method = "average", k = 4)
# run for k = 2 and k = 4 from above
print(boot_hc1)
# high stability for both clusters
# "dissolved" indicates no of times the bootstrapped cluster 
# is smaller or equal than the value set for dissolve (default = 0.5)
# "recovered" indicates no of times the bootstrapped cluster 
# is higher or equal than the value set for recover (default = 0.75)
# plot frequency distribution of bootstrap results
plot(boot_hc1)
print(boot_hc2)
# similar stability for 3 of the 4 clusters

# overall, 2 or 4 clusters could be justified based on the indices
# the data actually comes more or less from two groups:
# upstream and downstream of a salinity discharge point
# we can use this information to compute external validity

##### computation of external validation index: Rand index ######
extindi_2 <- cluster.stats(werra_dist, clustering = cut_2, alt.clustering=as.numeric(werra_env$Position))
extindi_2$corrected.rand
# reasonably high (in terms of Jaccard index) 
cut_2
as.numeric(werra_env$Position)
# only 1 case misclassified (it does not matter that 1 and 2 are interchanged)
# however, the upstream/downstream position is not necessarily the "truth"
# the cluster solution can be regarded as a different truth
# in this sense, external validation is not required

extindi_4 <- cluster.stats(werra_dist, clustering = cut_4, alt.clustering=as.numeric(werra_env$Position))
extindi_4$corrected.rand
# index logically much lower, given that external data consists of 2 groups

# visual comparison of external data
# using nmds results from above
# set colors for up-/downstream
cols <- c("darkred", "steelblue")
# x11() on windows
quartz()
ordiplot(werra_mds, type = "n")
points(werra_mds, col = cols[werra_env$Position], pch = 16)
legend("bottomleft", pch = 16, 
	col = cols, 
	legend = c('downstream', 'upstream'))

# next we add our cluster solution
ordispider(werra_mds, factor(cut_2), label = TRUE)
ordihull(werra_mds, factor(cut_2), lty = "dotted")
# one point is misclassified! But see comment above:
# upstream/downstream position is not necessarily 
# truth that determines ecological similarity


######################################################
# 	Extra Section: Don't in statistical analysis 	   #
#    Check cluster analysis with ANOVA 	             #
###################################################### 

# Although advocated in some text books,
# you must not use the group vector in a MANOVA/PERMANOVA 
# To check for significant differences in the cluster solution 
# see warning example below

# translated from Matlab code
# https://gist.github.com/mrkrause/2b315222abd00c902a1d

# This is a warning example that statistical tests on 
# results from the cluster analysis are significant even for 
# complete random data

# we construct clusters from completely random data
nrep <- 1000
sam_size <- 1000;
# sample size
groups_k <- 5; 
 
p_random <- c(NA)
# construct vector with NA
p_cluster <- c(NA)
# construct vector with NA
 
for(i in 1:nrep)
		{ 
		data <- runif(sam_size)	
		# sample from random uniform distribution
		cluster_group <- kmeans(data, groups_k)	
		# k means clustering
		random_group <- sample(1:groups_k, 1000, replace=TRUE)
		# assign random groups	
		p_random[i] <- summary(aov(data ~ random_group))[[1]][["Pr(>F)"]][1]
		# p value for random group
		p_cluster[i] <- summary(aov(data ~ cluster_group$cluster))[[1]][["Pr(>F)"]][1]
		# p value for clusters from random data
			}
#
  
par(mfrow=c(1,2))
hist(p_random, breaks=c(seq(0,1,0.025)), main = "Random data", ylab ="No of runs", xlab="p value from ANOVA")
# 
sum(p_random < 0.05)/1000
# as expected about 5% of random groups are significant  
 
hist(p_cluster, breaks=c(seq(0,1,0.025)), main = "Cluster data", ylab ="No of runs", xlab="p value from ANOVA")
sum(p_cluster < 0.05)/1000
# by contrast, more than 90% of clusters from complete random data are significant


##########################
#  End	Extra Section    #
#						 #
##########################

###################
#   Exercise      # 
###################
# a) Conduct a cluster analysis for the glass data set. 
# You should standardize the environmental variables
# before analysis using the scale() function. 
# Compare the results for "complete", "single" and "average" linkage. 
# Check the correlations of the cophenetic matrix with the initial distance matrices. 
# Evaluate the number of clusters for average linkage using the GAP index for up to k = 20

library(chemometrics)
data(glass)

# b) We have groups for the glass data, which can by accessed through:
data(glass.grp)
glass.grp

# Validate the cluster solution exernally, i.e. check the match between the
# given groups and the groups from cluster analysis
# Extract the cluster grouping vector for average linkage
# and check consistency with the given groups using the adjusted Rand index

######################
# k-means Clustering #
######################

library(chemometrics)
data(glass)
sc_glass <- scale(glass)
# scale glass data

# conduct k means clustering
set.seed(1000)
kclus <- kmeans(sc_glass, centers = 4, iter.max = 1000, nstart = 100)
# n start gives the number of random assignments at the start
# and iter.max the maximum number of iterations. 
# Exercise: Change the values by orders of magnitude and check what happens

# compute euclidean distance based nmds
dist_glass <- dist(sc_glass)
glass_mds <- metaMDS(dist_glass)

# Plot NMDS for cluster solution using the colours based on k-means
groups <- levels(factor(kclus$cluster))
ordiplot(glass_mds, type = "n")
cols <- c("steelblue", "darkred", "darkgreen", "pink")
for(i in seq_along(groups)){
  points(glass_mds, select = factor(kclus$cluster) == groups[i], col = cols[i], pch = 16)
}
ordispider(glass_mds, factor(kclus$cluster), label = TRUE)
ordihull(glass_mds, factor(kclus$cluster), lty = "dotted")

# calculate internal validation indices for different ks
# the function cascadeKM is a wrapper for k-means
# and computes the Calinski and Harabasz index

# using the calinski criterion
mult_k_clus <- cascadeKM(sc_glass, inf.gr=2, sup.gr=10, iter=50, criterion="ca")
# runs k means for different numbers of groups
# we reduce the number of iterations to save time
# inf.gr = min k; inf.gr = max.k

par(cex=1.5, las=1)
plot(mult_k_clus)
# would suggest 2 groups 

# you can access the partitioning into the different clusters with
mult_k_clus$partition[ , 1]
# gives the clustering, i.e. partitioning for two groups

colnames(mult_k_clus$partition)[1]
# gives the name for the partitioning

mult_k_clus$partition[ , 2]
# gives the clustering, i.e. partitioning for three groups

colnames(mult_k_clus$partition)[2]
# gives the name for the partitioning

################
#   Exercise   #
################
# Compare the solutions of 2 and 4 groups
# using two other CVIs (note: use the clustering!) and visualise for 2 and 4 groups
# Which solution would you regard as the most appropriate?

