###################################################################################
# Script for Session 5 of the course: Applied multivariate Statistics with R  	  #
# 						by Ralf B. Schäfer,	WS 2018/19							  #
# 			Similarity measures and Non-metric multidimensional scaling			  #
###################################################################################

# we first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine here
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

library(vegan)
# calculation of association measures
# we use an example matrix with 6 species observed at 5 sites
# and calculate some distance measures for this matrix
mat1 <- matrix(data = c(40, 5, 12, 0, 0, 40, 10, 1, 18, 15, 22, 3, 100, 5, 22, 16, 200, 50, 2, 1, 1, 0, 0, 0, 40, 10, 0, 0, 0, 20), nrow = 5, byrow = TRUE)
mat1


##############
# Exercise  #
##############
# a) Compute the Bray-Curtis, Euclidean and Jaccard dissimilarities/distances.
# What are the differences between the coefficients regarding the relationship
# of the highest to lowest dissimilarity/distance? What else do you observe?
# Use the function vegdist() to calculate dissimilarity/distance and check, which arguments you need to provide.
# Be careful when calculating Jaccard dissimilarity (you need to set a further argument)!

# b) How does standardisation of the data affect the results?
# Use decostand(yourmatrix, method="max") and recalculate the dissimilarities/distances.
# This will divide each observation by the maximum value of each variable (species).

############
# Solution #
############

######################################
# Nonmetric multidimensional scaling #
######################################

# function to check for number of axis in NMDS
# NMDS is repeated 10 times for each dimension to
# check for the influence of the random initial configurations
# x must be a data-object which can be coerced into a matrix;
# max is the maximum number of dimensions
# dist_meas is the distance measure
scree_values <- function(x, max, dist_meas) {
  xx <- as.matrix(x)
  scree.points <- NULL
  scree.temp <- 1
  for (i in 1:max)
  {
    sol.mds <- NULL
    sol.mds <- replicate(
      10, metaMDS(xx, k = i, trymax = 30, distance = dist_meas)$stress,
      simplify = TRUE
    )
    scree.points <- append(scree.points, sol.mds)
  }
  return(scree.points)
}
# see ?metaMDS to change according to your needs

# We use the Dune vegetation data
library(vegan)
data(dune)
head(dune)
summary(dune)
# data in similar order of magnitude
# -> no standardisation required

# calculate relationship between dimensions and stress
scree.res <- scree_values(dune, max = 10, dist_meas = "bray")

# plot dimensions vs. stress
par(cex = 1.5, las = 1)
plot_vec <- as.vector(sapply(1:(length(scree.res) / 10), function(x) rep(x, 10), simplify = "vector"))
# prerequisite for plotting replicated runs of nmds
plot(plot_vec, scree.res, ylab = "Stress1", xlab = "dimensions")
lines(1:(length(scree.res) / 10), as.vector((by(scree.res, plot_vec, mean))))
abline(h = 0.10, col = "red")
# red line indicates threshold for good Stress1 value

# we select two dimensions as Stress level is around 10
# if distance measure is not specified bray curtis is used
specnmds <- metaMDS(dune, k = 2)

# plot Shepard diagram
stressplot(specnmds)
# Non-metric fit is based on stress and calculated as R^2 = sqrt(1- S^2)
# Linear R^2 = explained variation of regression of observed distances against ordination distances

# preparation for plotting
sumcols <- colSums(dune)
# calculate sum of columns. If two objects are plotted
# above each other in the ordination, select the object that has higher column sums
par(cex = 1.5)
plot(specnmds, dis = "sp", type = "n")
nmdsplot <- orditorp(specnmds, display = "sp", priority = sumcols, pcol = "gray20", pch = "+", cex = 0.8)
# species scores are obtained by weighted averaging as for RDA (using WA scores)

# fitting of environmental variables on the NMDS ordination results
data(dune.env)
head(dune.env)

# we transform some factor variables into continuous numeric variables
dune.env$Moisture <- as.numeric(as.character(dune.env$Moisture))
dune.env$Manure <- as.numeric(as.character(dune.env$Manure))

# we combine some factor levels
combine <- c("BF", "HF")
levels(dune.env$Management)[levels(dune.env$Management) %in% combine] <- paste(abbreviate(combine, 5), collapse = "&")

# fitting of environmental variables to result
ef <- envfit(specnmds, dune.env[, c(1, 2, 3, 5)], permu = 500)
ef
# envfit() fits variables to gradients in the ordination
# The position of arrows and centroids is given, when calling the object
# The R^2-values result from the regression of variables against ordination scores
# The p-values result from permutation: randomly permuted vectors or factors are generated
# and for each a R^2-value is calculated. These R^2 values are subsequently used
# to calculate the percentage of permuted vectors/factors exhibiting higher R^2-values
# than the vector/factor under scrutiny.
# Example:
# After 500 permutations, 5% of the R^2-values of a randomly permuted vector/factor (i.e. 25)
# have a R^2 ≥ 0.4. This means that vectors/factors with R^2 > 0.4 are significant at the level of 0.05.

plot(ef, p.max = 0.05, axis = TRUE, col = "red")
# scaling of arrows by correlation coefficient

# Note that the response is not necessarily linear
# as the ordination is unconstrained from environmental variables
par(cex = 1.5)
ordisurf(specnmds ~ Moisture, data = dune.env, plot = T, bubble = 4)
# fits surface for Moisture
# Moreover, the vectors/factors are fitted to the ordination distances
# that are not the true distances between objects
# -> significant relationship may not necessarily be ecologically meaningful

ef <- envfit(specnmds, dune.env[, c(2)], permu = 500)
plot(ef)
# compare with simple vector, can give misleading impression

## There is also a function to explore 3 dimensional ordination plots
# in the library vegan3d
library(vegan3d)

# in case you have not installed this library, run
# install.packages("vegan3d")
# and then run
# library(vegan3d)

specnmds2 <- metaMDS(dune, k = 3)

par(cex = 2)
ordirgl(specnmds2, type = "t")

#############
# Exercise   #
#############
# Conduct a NMDS for the lichen data (varespec (data is in vegan package))
# and examine differences in the results for the euclidean and bray-curtis distance
# (see metaMDS function on how to change the distance measure).
# Fit environmental variables to the resulting ordination plot.
# Inspect intercorrelation between variables before fitting and
# remove highly correlated variables based on expert knowledge (or by gut feeling :-)).

################################################
#   Extra Section: Influential points		   #
#   Warning example provided by Jari Oksanen   #
################################################

library(labdsv)
# bryce canyon vegetation data
data(bryceveg)
nrow(bryceveg)
# 160 observations

# environmental data for the sites
data(brycesite)

# conduct NMDS
m <- metaMDS(bryceveg)
par(cex = 1.5)
ordisurf(m ~ slope, brycesite, knots = 10, bubble = 4)
# fits trend surface and slope

(ef_2 <- envfit(m ~ slope, brycesite))
# highly significant

plot(ef_2)
# direction seems to be influenced by a few points

# we repeat the analysis and removee four sites of 160
brycesite$slope
# a couple of sites have very high values for slope

m2 <- metaMDS(bryceveg[!brycesite$slope > 34, ])
# remove all sites with slope > 34

# conduct NMDS
par(cex = 1.5)
ordisurf(m2 ~ slope, brycesite[!brycesite$slope > 34, ], knots = 10, bubble = 4)
(ef_3 <- envfit(m2 ~ slope, brycesite[!brycesite$slope > 34, ]))
# not significant anymore after removal of 4 sites
plot(ef_3)
# trend not convincing anymore

###########################
# End Extra Section		  #
###########################
