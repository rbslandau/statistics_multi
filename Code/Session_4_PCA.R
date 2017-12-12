###################################################################################
# Script for Session 4 of the course: Applied multivariate Statistics with R  	  #
# 						by Ralf B. Schäfer,	WS 2017/18							  #
# 						Principal Component Analysis						      #
###################################################################################

# # code for figures on lecture slides - uncomment to run
# ### create normal distribution of sensitivity
# par(cex= 2, las=1, mar=c(5,4,1,1))
# norm_dist <- dnorm(seq(1,200,1), 50, 20)
# plot(norm_dist, type = "l", xlim = c(0,200), ylab = "Species density", xlab = "Environmental variable", lwd=2)
# # density plot

# # for plot with multimodal distribution
# norm_dist2 <- dnorm(seq(1,200,1), 60, 150)
# norm_dist3 <- dnorm(seq(1,200,1), 150, 40)
# plot(norm_dist+norm_dist2+norm_dist3 , type = "l", xlim = c(0,200), ylab = "", xlab = "", lwd=2, axes = F)
# axis(1, labels = F)
# axis(2, labels = F)
# mtext("Environmental gradient", side = 1, cex =2, line = 1)
# mtext("Population density", side = 2, las = 0, cex = 2, line = 1)

# # for plot with polynomial
# library(polynom)
# p <- polynomial(c(91194.69,-1417,1 ))
# # define polynomial equation for lambda
# solve(p)
# par(cex=2, las=1)
# plot(p, axes=F, lwd=1.5, xlab="", ylab=expression(paste("F(", lambda, ")")))
# axis(1, pos=0)
# axis(2, pos=0)
# mtext(expression(lambda), side=1, cex=2)

# let us first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine here
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

library(vegan)
library(mvoutlier)
# the vegan package contains several functions
# for multivariate analysis of ecological data
# and represents a standard package
data(varechem)
head(varechem)
# environmental data from soils

# check for multivariat normality
chisq.plot(varechem, quan = 1, ask = FALSE)
# looks relatively normally distributed, some deviation
# you could screen univariate distributions to check
# whether a particular variable contributes to the deviation

# principal component analysis
(va_pca <- rda(varechem, scale = TRUE))
# you should always scale unless the descriptors have the same order of magnitude
summary(va_pca, display = NULL)

# first two axis represent 60% of variation
par(cex = 1.5)
plot(va_pca, display = c("sites"), type = "text", scaling = 1) #
# we can interpret the distances between sites in this plot

# we check the overall variance of all variables
sum(apply(scale(varechem), 2, var))
# overall variance is 14 - because we have 14 environmental variables

# now check the eigenvalues
va_pca$CA$eig -> ev
# extraction of eigenvalues
# variance is preserved in our new eigenvalues:
sum(ev)
# sum of eigenvalues is 14

###########################
# How many axes are needed?#
###########################

# Apply sum criterion
alpha <- 0.7
summary(va_pca)$cont$importance[3, ] > alpha
# 2 axes needed in order to explain 70% of variance

# Apply screeplot
par(mfrow = c(1, 1))
screeplot(va_pca, type = "lines")
# in a screeplot we look for a separation between important and non-important components
# as indicated by an elbow.
# The name derives from the illustration of a mountain from which stones drop
# and we search for the end of the mountain and the beginning of the scree
# for the given data, this is not very obvious, but the reduction from PC5 to
# PC6 is minor and consequently 5 PCs could be selected based on the plot

# Apply broken stick criterion
screeplot(va_pca, type = "lines", bstick = TRUE)
# third component slightly lower than broken stick criterion
# broken stick for individual PCs and cumulative variance provided in the following library:
library(BiodiversityR)
PCAsignificance(va_pca)
# PCA axes with a higher percentage of (accumulated) variance than the broken-stick variance are regarded as significant.
# Only the first two axes fulfill broken stick criterion and eight axes the accumulated variance criterion

###########################
#     Extra Section 	  #
###########################

## illustration of broken stick model in R
# we look at the case for p = 2, i.e. the stick is broken into 2 pieces
p <- 2
# we determine the length of the longest piece y for a random breaking of a stick of unit length,
# hence we sample from a uniform distribution between 0 and 1.
# If the selected value y is larger than 0.5, we keep it, otherwise 1-y would be the longest piece
target_vec <- c(1:10000)
# repeat 10000 times
# set seed for reproducible example
set.seed(length(target_vec))
for (i in 1:length(target_vec)) {
  y <- runif(p - 1)
  target_vec[i] <- ifelse(y >= 0.5, y, 1 - y)
}

# if else rule as outlined above
mean(target_vec)
# the model predicts a length of 0.75 (equals the variance) for the longest piece (first component) in the case of two variables

###########################
#  End Extra Section 	  #
###########################

# Apply k-fold cross validation to select number of PCs
# see help and Josse J. & Husson F. (2012) Computational Statistics & Data Analysis 56, 1869–1879.
# for details
library(missMDA)
set.seed(100)
estim_ncpPCA(varechem, method.cv = "Kfold", nbsim = 100, pNA = 0.25)
# pNA is 1/k, i.e. 4-fold cross validation in this case
# minimum MSEP for 2 PCs, though only slightly lower than 3 PCs

# note that per default only the first 5 PCs (argument ncp.max =5) are considered
# you can set a different ncp.max
# note that this function can also be used for data sets with missing values
# see function "imputePCA" for details

# the explained variance in cross validation
# may be preferred as metric if most variance should be preserved
# for example, you could select the number of variables within 1 standard error
# from the maximum explained variance

# the function pcaCV provides the explained variance and conducts
# casewise-CV in contrast to the function ncpPCA
# that randomly distributes missing values across the matrix.
# Using the functions of this packages requires
# a PCA object as fitted with the function princomp
X.rpc <- princomp(varechem, cor = TRUE)
library(chemometrics)
pca_cv <- pcaCV(varechem, scale = TRUE, segments = 4, plot.opt = TRUE)
# the argument segment sets k in k-fold cross-validation
# 8 PCs required to explain 80% of variance
# note the difference to the sum criterion, which is overly optimistic
# because it does not consider prediction to new data

###########################
#    Outlier checking	  #
###########################

par(mfrow = c(1, 2), cex = 2)
res <- pcaDiagplot(varechem, X.rpc, a = 2, quantile = 0.975)
# we set the number of PCs to 2 based on MSEP above (argument a =2)

# we look for orthogonal distances and score distances
which(res$ODist > 4)
# gives rowname and the row of outlier
# (site 11 in row 23 of varechem would be an outlier regarding orthonormal distance)
# since this observation is only an outlier regarding the distance and not score
# (and not extreme compared to the other points), we do not take action

################
#    Biplot	   #
################

# we can display the sites and variables simultaneously
par(mfrow = c(1, 3), cex = 1.2)
biplot(va_pca, scaling = 3, display = c("sp", "site"))
# the scores of the vectors on the axes indicate
# the relationship with the component axis
# higher score means higher relationship
# the biplot displays the following relationships:
# between sites, between sites and environmental variables
# (orthogonal projection on variables)
# and between environmental variables (correlation and maximum rate of change associated with variable)
# read chapter 3 in the following document
# to understand the impact of axis scaling
vegandocs("decision-vegan")
# Borcard (2011, p.121) states: "Bottom line: if the main interest of the analysis is
# to interpret the relationships among objects [here sites, RBS ], choose scaling 1.
# If the main interest focuses on the relationships among descriptors [here env. variables, RBS ], choose scaling 2"
# scaling = 3 represents a compromise, but note that both objects and descriptors are not presented with original distances/correlations
# see also Legendre & Legendre 2012: 639-641
# For mathematical background of scaling in biplots see Jolliffe (2002): 90ff
# Jolliffe I.T. (2002) Principal component analysis, 2nd edn. Springer, New York.

#########################################################################
# How well are the variables captured by the selected number of axes?  #
#########################################################################
# How well are the variables captured by the axes?
pcaVarexpl(varechem, a = 2)
# this function calculates the explained variance of each variable
# for the number of PCs
pcaVarexpl(varechem, a = 14)
# of course, for 14 PCs (= no. of vars), we obtain an R2 of 1 for all variables

# loadings provide the weight of the original variables
# in the computation of PCs
(loadings <- X.rpc$loadings)
# same result for rda object
va_pca$CA$v

# a * a^t results in the identity matrix
va_pca$CA$v %*% t(va_pca$CA$v)
# note that due to rounding errors
# the values do not yield to 0 above and below
# the diagonal

# variance of individual variables is obtained by
apply(va_pca$CA$v ^ 2, 1, sum)

# Correlation loadings can aid in the interpretation
# of the relationship between the original variables
# and the PCs
# they are obtained by multiplying the eigenvectors with
# the standard deviation of the individual PCs
sdev <- X.rpc$sdev
var_cor_func <- function(var.loadings, comp.sdev) {
  var.loadings * comp.sdev
}
var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
var.cor
# several variables exhibit high correlations with first and subsequent axes

# same result as when correlating the original variables with the new axes
pca_axes <- vegan::scores(va_pca, disp = "sites", choices = c(1:ncol(varechem)), scaling = 0)
# further details on the extraction of PC axes is given below
cor(scale(varechem), pca_axes)

########################################################
# Simplify variable interpretation through sparse PCA  #
########################################################
# we note that most variables load at least to a minor extent
# on all PCs. This complicates interpretation
# Sparsity methods limit the number of variables that load on a PC

library(pcaPP)
# we estimate the optimal penalty term first
k.max <- 2
##  k.max = max number of considered sparse PCs
oTPO <- opt.TPO(scale(varechem), k.max = k.max, method = "sd")

oTPO$pc
# the model selected by opt. TPO
oTPO$pc$load
# and related sparse loadings

##  Tradeoff Curves: Explained Variance vs. sparseness
par(mfrow = c(1, k.max))
for (i in 1:k.max) plot(oTPO, k = i)
# L0 gives number of variables with zero loadings on PC
# lambda opt is the estimated optimal penalty
# ECV1 is the empirical cumulated variance of the PCs

##  Tradeoff Curves: Explained Variance vs. lambda
par(mfrow = c(1, k.max))
for (i in 1:k.max) plot(oTPO, k = i, f.x = "lambda")
# lambda represents penalty term

#  Objective function vs. lambda
# we do not go into details here, but if you are interested
# you could look at the lambda against the TPO function
# by uncommenting the next two lines
# par (mfrow = c (1, k.max))
# for (i in 1:k.max)        objplot (oTPO, k = i)

# loadings
oTPO$pc$load
# calculation of correlation loadings
t(apply(oTPO$pc$load, 1, var_cor_func, oTPO$pc$sdev))
# eigenvalues of PCs
oTPO$pc$sdev[1:2] ^ 2

# use optimized lambdas to compute sparsePCA
spc <- sPCAgrid(scale(varechem), k = k.max, lambda = oTPO$pc.noord$lambda, method = "sd")

par(mfrow = c(1, 1))
biplot(spc)
# all variables now pushed on direction of PCs

##############################
# Extraction of PC scores    #
##############################

# scores can be used in PC regression
# here we select 4 axes
pca_axes <- vegan::scores(va_pca, disp = "sites", choices = c(1:ncol(varechem)), scaling = 0)
# scaling depends on whether we focus on sites or variables, 0 means raw scares (not scaled)
# we can address the axes as usually
pca_axes[, 1:4]

# The first PCs preserve what is common in the data, not necessarily what predicts
# a gradient in an external data set (or what differentiates groups)
# Especially when the aim is prediction, this may lead to problems, if only
# the first axes are included in the analysis as later axes may still be important
# see discussion in Jolliffe (2002), p173f. An example is discussed where the PCs 8 and 10
# bear almost no variance (3 and < 1 %) but are among the most important predictors
# for a response variable

# also have a look at this blog:
# http://www.r-bloggers.com/pca-or-polluting-your-clever-analysis/

# extraction of sparse PCA scores via:
spc$scores
# the scores could subequently be used in PC regression

########################################
#     		Extra Section 	 		   #
# Influence of non-linearity on PCA    #
# 									   #
########################################

# in this example we inspect the PCA results
# if the gradient is non-linear, i.e. for species data
# setup species with unimodal gradient
hss <- c(1, 2, 4, 7, 8, 7, 4, 2, 1)
# unimodal gradient
spec1 <- c(hss, rep(0, 10)) #
spec2 <- c(rep(0, 5), hss, rep(0, 5))
spec3 <- c(rep(0, 10), hss)
# display three species against sites
par(mfrow = c(1, 2), cex = 1.6)
plot(spec1, col = "blue", type = "o", ylab = "Species abundance", xlab = "Site", main = "Three species against environmental gradient")
points(spec2, col = "red")
lines(spec2, col = "red")
points(spec3, col = "magenta")
lines(spec3, col = "magenta")
# conduct PCA
data <- cbind(spec1, spec2, spec3)
pca <- rda(data)
biplot(pca)
# suggests that there are 2 or 3 gradients in the data
# however, we have constructed our data to have only one gradient;
# compare Legendre 2012 (p. 483)
par(mfrow = c(1, 1), cex = 1.6)
screeplot(pca, type = "lines", bstick = TRUE)
# first axis not higher than broken stick model

par(mfrow = c(1, 2), cex = 1.6)
pcaVarexpl(data, a = 1)
# first axis tries to capture species one and three
pcaVarexpl(data, a = 2)
# second axis tries to capture species 2 and remaining variance of other species

# to understand how this happens we inspect the bivariate case:
data2 <- cbind(spec1, spec2)
pca2 <- rda(data2)
biplot(pca2)
plot(spec1, spec2, type = "n")
text(spec1, spec2, labels = 1:19)
# pattern cannot be catched linearly

########################################
#     	End	Extra Section 	 		   #
# 									   #
########################################


##################
#   Exercise     #
##################
# Conduct a complete PCA for the glass data set.
# How many principal components are needed to visualise the data?
# Extract the number of meaningful pcs and check the intercorrelation of these axes.
# What do you observe?
# Which variables are most important regarding the first axis? Explore with a sparse PCA.
library(chemometrics)
data(glass)
