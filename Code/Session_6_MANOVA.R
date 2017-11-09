###################################################################################
# Script for Session 6 of the course: Applied multivariate Statistics with R  	  #
#						by Ralf B. Schäfer,	WS 2017/18							  #
# 						Hotelling T2 test and MANOVA					  		  #
###################################################################################

# let us first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine here
# to simplify the identification of your path, you can use the following function
# file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

tibet_data <- read.table("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/tibet_skull/at_download/file", sep = ";", dec = ".", header=TRUE)
head(tibet_data)
str(tibet_data)
# taken from Everitt 2005
# Tibetan skull data: Group one (type I) found in graves in Sikkim and neighboring areas of Tibet; 
# group two (type II) consisting of the remaining 15 skulls picked up on battlefield in the Lhasa district
# Hypothesis (in former times): Tibetans from Khans might be survivors of a particular fundamental human type
# unrelated to the Mongolian and Indian types that surrounded them.
skull_dat <- tibet_data[ , 1:5]

##################################################
# check model assumptions for Hotelling T2 test  #
##################################################

# check assumption of multivariate normality 
library(mvoutlier)
chisq.plot(skull_dat[tibet_data$Type == 1, ], quan = 1, ask= FALSE)
# looks ok, given the small sample size
chisq.plot(skull_dat[tibet_data$Type == 2, ], quan = 1, ask= FALSE)
# looks ok, given the small sample size

# check assumption of homogeneity of within group covariance matrices 
# using a hypohtesis test
library(vegan)
skull_dist <- dist(skull_dat)
# computes euclidean distance matrix
(skull_mhv <- betadisper(skull_dist, tibet_data$Type))
# using a homogeneity of variance test introduced by Anderson (2006)
set.seed(111)
permutest(skull_mhv)
# no indication of differences using permutational test

############################ 
#    Visual checking 	   #
############################

# visual method to check for homoscedasticity of covariances
library(pcaPP)
par(cex=1.5)
plotcov(cov(skull_dat[tibet_data$Type == 1, ]), cov(skull_dat[tibet_data$Type == 2, ]), method1 = "Group 1", method2 = "Group 2")
# shows the covariances of both samples in a plot. 
# except for one or two covariances, all look very similar

# check homoscedasticity of variances
# as for anova

# for first variable
# put on same median to check homogeneity of variance
med <- tapply(tibet_data$Length, tibet_data$Type, median)
# compute median of groups
w <- tibet_data$Length - med[tibet_data$Type]
# compute distance to median and plot
plot(w ~ Type, data = tibet_data)

# for second variable
# put on same median to check homogeneity of variance
med <- tapply(tibet_data$Breadth, tibet_data$Type, median)
# compute median of groups
w <- tibet_data$Breadth - med[tibet_data$Type]
# compute distance to median and plot
plot(w ~ Type, data = tibet_data)

# would be required for all variables

##########################
#  End	Visual checking  #
#						 #
##########################

#########################
#   Hotelling T2 test   #
#########################

# the function for the Hotelling T2 test is directly available in a R package
library(ICSNP)
HotellingsT2(as.matrix(skull_dat) ~ tibet_data$Type)
# skull_dat is a data.frame (can be checked with class(skull_dat)
# but function requires matrix as input data

manov_skull <- manova(as.matrix(skull_dat) ~ tibet_data$Type)
summary.manova(manov_skull, test = "Hotelling-Lawley") 
# you obtain the same result for a manova, when selecting the Hotelling-Lawley test statistics

### in the case of strong deviations from the multivariate normal distribution
# an alternative non-parametric tests can be used
library(cramer)
cramer.test(as.matrix(skull_dat[tibet_data$Type == 1, ]), as.matrix(skull_dat[tibet_data$Type == 2, ]))
# a new test for multivariate two sample problem, introduced by
# Baringhaus, L. and Franz, C. (2004) On a new multivariate two-sample test, 
# Journal of Multivariate Analysis, 88, p. 190-206

##################
## MANOVA section#
##################

# conduct manova on soil data from the following study:
# Holmgren et al. (2015) Positive shrub-tree interactions 
# facilitate woody encroachment in boreal peatlands. 
# Journal of Ecology, 103(1): 58-66
# You can download the data here:
# http://datadryad.org/resource/doi:10.5061/dryad.jf2n3/2
# download the file Table 1_environ...
# and simplify filename to Table_1.xls and
# move file into your working directory
# open in a spreadsheet program and save as .csv format

soil_dat <- read.csv("https://raw.githubusercontent.com/rbslandau/statistics_multi/master/Data/Table_1.csv", sep = ";", dec = ".")
head(soil_dat)
# many columns empty
summary(soil_dat)
levels(soil_dat$Hummock)
# restrict to levels A, B, C and D
soil_dat2 <- soil_dat[soil_dat$Hummock %in% c("A", "B", "C", "D"), ]
# check data
soil_dat2
# we remove non-needed columns and rows and rename columns 
soil_dat3 <- soil_dat2[ , -c(2:5, 9: ncol(soil_dat2))]
soil_dat3 
names(soil_dat3) <- c("Cat_", "N", "P", "K")
# remove NAs in variables
soil_dat4 <- soil_dat3[! is.na(soil_dat3$N), ]

# extract quantitative variables and scale
# otherwise results dominated by variable with highest variance
soil_fin <- scale(soil_dat4[ ,2:4])

# check assumption of multivariate normality 
chisq.plot(soil_fin[soil_dat4$Cat_ == "A", ], quan = 1, ask= FALSE)
# Looks acceptable
chisq.plot(soil_fin[soil_dat4$Cat_ == "B", ], quan = 1, ask= FALSE)
# Looks acceptable
chisq.plot(soil_fin[soil_dat4$Cat_ == "C", ], quan = 1, ask= FALSE)
# Looks acceptable
chisq.plot(soil_fin[soil_dat4$Cat_ == "D", ], quan = 1, ask= TRUE)
# observation 1 and 13 are outlier
# looks acceptable after removal of outlier

soil_fin[soil_dat4$Cat_ == "D", ]
# observation with rowname 31 and 113 is outlier
soil_fin <- soil_fin[! row.names(soil_fin) %in% c(31,113), ]
soil_dat4 <- soil_dat4[! row.names(soil_dat4) %in% c(31,113),  ]
# check again
chisq.plot(soil_fin[soil_dat4$Cat_ == "D", ], quan = 1, ask= FALSE)
# looks better now

# check assumption of homogeneity of within group covariance matrices 
soil_dist <- dist(soil_fin)
# computes euclidean distance matrix
(soil_mhv <- betadisper(soil_dist, soil_dat4$Cat_))
# using a homogeneity of variance test introduced by Anderson (2006)
set.seed(111)
permutest(soil_mhv)
# generally lower dispersion of C and D compared to A and B
# this is statistically significant -> deviation from homogeneity of covariance
# we have to consider this in interpretation

# look at collinearity
par(cex=1.3)
plotcov(cov(soil_fin[soil_dat4$Cat_ == "B", ]), cov(soil_fin[soil_dat4$Cat_ == "C", ]), method1 = "Group B", method2 = "Group C")
# similar covariances, N and P higher correlation for group B
plotcov(cov(soil_fin[soil_dat4$Cat_ == "A", ]), cov(soil_fin[soil_dat4$Cat_ == "D", ]), method1 = "Group A", method2 = "Group D")
# similar covariances, high correlations for N and P
# Despite some differences, relatively similar covariance structure

# manova command
manov <- manova(as.matrix(soil_fin) ~ soil_dat4$Cat_)
# requires matrix as input data
summary.manova(manov) 
# gives the manova results, the test can be selected using 
# test = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy") argument
# Pillai is the default
summary.manova(manov, test = c("Wilks"))
summary.manova(manov, test = c("Roy"))
# Pillai and Wilks not significant
# Roy with p = 0.03 significance
# but Roy is strongest influenced by deviations from 
# homogeneity of covariance matrices

# Note that we have a highly balanced design here
summary(soil_dat4$Cat_)
# hence deviations from assumptions have minor influence
# on Pillai

# for Pillais trace we obtain an average R2 per descriptor
# by dividing Pillais trace by the number of eigenvalues (i.e. variables)
man_eigen <- summary.manova(manov)$Eigenvalues
l_eigen <- length(man_eigen)
pill_man <- summary.manova(manov)$stats[1,2]
pill_man/l_eigen
# average R2 per variable is 0.054

summary.aov(manov)
# gives the results for the single anovas - significant for N
# note that the p-values are not corrected for multiple testing here.
# in addition, the relationships between variables are ignored

### non parametric MANOVA in the case of strong deviations from normality
library(ICSNP)
rank.ctest(as.matrix(soil_fin) ~ soil_dat4$Cat_)
# yields to same results

## robust MANOVA 
# see http://www.statistik.tuwien.ac.at/public/filz/papers/09CSDAmanova.pdf 
# for a paper on robust MANOVA by Peter Filzmoser. 
# In general, Peter Filzmoser provides many papers on robust statistical methods using R
# http://file.statistik.tuwien.ac.at/filz/papers/09CSDAmanova.pdf
# the classical, robust and rank-based MANOVA with the Wilk test statistic
# can be conducted using:
# Wilks.test {rrcov}

## MANOVA using RDA
# first step would be checking for homogeneity of covariance matrices
# has already been done above
manov_rda <- rda(soil_fin ~ soil_dat4$Cat_)
summary(manov_rda)
# explained variance of factor is 6.2 %
anova(manov_rda, step =1000, perm.max=1000)
# result similar to Manova - but a permutation test is used
summary.manova(manov) # only for comparison

# the RDA context offers the opportunity of plotting
par(cex=1.5)
# plot variables
plot(manov_rda, scaling=3, display=c("sp"), xlim=c(-1.5,2), ylim=c(-1.5,2.2))
# plot centroids of factors
text(manov_rda, display="cn", scaling=3, labels= c("A", "B", "C", "D"), col = c("green", "red", "light blue", "orange"))
# add individual sites
wa_sc <- scores(manov_rda, display="wa", scaling=3)
text(wa_sc[soil_dat4$Cat_=="A" , ], col="green", labels = row.names(wa_sc[soil_dat4$Cat_=="A", ]))
text(wa_sc[soil_dat4$Cat_=="B" , ], col="red", labels = row.names(wa_sc[soil_dat4$Cat_=="B", ]))
text(wa_sc[soil_dat4$Cat_=="C" , ], col="light blue", labels = row.names(wa_sc[soil_dat4$Cat_=="C", ]))
text(wa_sc[soil_dat4$Cat_=="D" , ], col="orange", labels = row.names(wa_sc[soil_dat4$Cat_=="D", ]))
# shows that centroids of all factors are close to overall centre
# spread of A and B higher

##############
#  Exercise  #
##############
# Conduct a MANOVA for the glass data set.
# Do the glass vessels originate from different manufacturers?
# use ?glass to obtain more information on the data
# conduct a RDA on the groups as well and visualise the results.

library(chemometrics)
## data preparation
data(glass)
data(glass.grp)
# the data is extremely uneven in terms of sample size
# we only choose the groups with relatively even sample size
# and only a few variables

glass_dat <- glass[glass.grp != 1, 1:5]
glass_grp <- factor(glass.grp[glass.grp != 1])
# you have to define glass.grp as factor
summary(glass_grp)
