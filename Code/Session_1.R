#################################################################################
# Script for Session 1 of the course: Applied multivariate Statistics with R		#
#					by Ralf B. Schäfer											#
# 					WS 2017/18													#
# The script describes simulation approaches for basic statistical methods		#
# as well as basic exploratory techniques that are used in model diagnostics 	#
#################################################################################

# let us first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

########################################
# Permutation: Inference via shuffling #
########################################

### we create two populations with a different distribution
# Rationale: We will use these data in permutational hypothesis testing
set.seed(20)
pop1 <- rnorm(1000, mean = 100, sd = 4)
# draw samples from a normal distribution
set.seed(20)
pop2 <- runif(1000, min = 50, max = 135)
# draw samples from a uniform distribution
# they have different means
mean(pop1)
mean(pop2)
# and clearly a very different distribution
par(mfrow = c(1,2))
hist(pop1)
hist(pop2)

# we sample with a relatively small size from the populations
# this translates, for example, to a monitoring study 
# with biometric measurements of a sample of the population 
no_sam1 <- 30
no_sam2 <- 30
set.seed(80)
samp1 <- sample(pop1, no_sam1)
mean(samp1)
set.seed(55)
samp2 <- sample(pop2, no_sam2)
mean(samp2)
par(mfrow=c(1,2))
hist(samp1, breaks = 15)
hist(samp2, breaks = 15)
# we create a categorical variable that codes for the two populations
categ <- c(rep("Samp1", no_sam1), rep("Samp2", no_sam2))
# and combine the two populations into one vector
data_vec <- c(samp1,samp2)
# create data.frame
datafr <- data.frame(data_vec, categ)
names(datafr) <- c("Size", "Pop_type")

# Imagine we have sampled from the two populations
# and now want to compare their body size
# we have the hypothesis that the body sizes differ between the two populations

# we create a function
meanDif <- function(x, grp) 
						{
						mean(x[grp == "Samp1"]) - mean(x[grp == "Samp2"])
						}
with(datafr, meanDif(Size, Pop_type))

# we use permutation for hypothesis testing
library(permute)
perm_vec <- numeric(length = 10000)
# Legendre & Legendre (2012) suggest to use at least 10,000 permutations
# for inference, especially when aiming at publishing the results 
N <- nrow(datafr)
set.seed(300)
for(i in seq_len(length(perm_vec) - 1)) # loop runs 9999 times
	{
	  perm <- shuffle(N)
	  perm_vec[i] <- with(datafr, meanDif(Size, Pop_type[perm]))
	}
# results in 9999 shuffled data sets
# add non-shuffled data	because we test the null hypothesis
# i.e. assume that the unpermuted statistic could be obtained under H0
perm_vec[10000] <- with(datafr, meanDif(Size, Pop_type))	  

# plot the results of the permutations
par(cex = 1.2, mfrow = c(1,1))
hist(perm_vec, main = "", breaks = 20, xlab = expression("Mean difference (Pop1 - Pop2)"))
rug(perm_vec[10000], col = "red", lwd = 2, ticksize = 0.5)

# compute permutation p-value
Dbig <- sum(perm_vec >= perm_vec[10000])
Dbig/length(perm_vec)
# permutation test identifies statistical significant difference

##############################
# Classical approach: t-test #
##############################

# classically you would use a t-test (you should be familiar with this test)
# we check the assumptions to decide whether we can assume variance equality

###################################################
## Boxplots for checking homogeneity of variance  #
###################################################

plot(Size ~ Pop_type, data = datafr)
# put on same median - code taken from plot.hov function
# this is only done to ease visual comparison
med <- tapply(datafr$Size, datafr$Pop_type, median)
# compute median of groups
w <- datafr$Size - med[datafr$Pop_type]
# compute distance to median and plot
# i.e. we move both data sets to 0 to ease comparison of the variance
plot(w ~ Pop_type, data = datafr)
# the figure displays the boxplots shifted to a median of zero
# variances clearly differ

# although we know that the assumption of normality is violated,
# we run a t-test anyways, accounting for different variances
(body_t2 <- t.test(Size ~ Pop_type, data = datafr, var.equal = FALSE))
# no statistical significance, p-value two-fold higher than for permutation

# how would we check for normal distribution?

##########################################################
# QQplots and boxplots to check for normal distribution  #
##########################################################

# we need to check for each sample
# quantile-quantile plot
qqnorm(datafr$Size[datafr$Pop_type =="Samp1" ], datax = TRUE)
# data fit a line relatively well, though some deviation is visible
qqline(datafr$Size[datafr$Pop_type =="Samp1" ], datax = TRUE)
# the function qreference in the DAAG package produces reference plots 
# to aid in the evaluation whether the data are normally distributed
library(DAAG)
qreference(datafr$Size[datafr$Pop_type =="Samp1" ], nrep = 8) 
# nrep controls the number of reference plots
# data does look perfectly normal

# similarly, you could try to select your qqplot from 
# a range of qq plots to see whether it deviates markably.
# The code is provided here:
# http://www.r-bloggers.com/checking-glm-model-assumptions-in-r/

qqnorm(datafr$Size[datafr$Pop_type =="Samp2" ], datax = TRUE)
qqline(datafr$Size[datafr$Pop_type =="Samp2" ], datax = TRUE)
# slightly stronger deviation from line
qreference(datafr$Size[datafr$Pop_type =="Samp2" ], nrep = 8) 
# despite we sampled from a uniform distribution,
# the deviation seems still moderate

# Another option to check normality is the histogram
# histograms can also be used to check for
# asymmetry in the distribution of variables (i.e. skewness)

hist(datafr$Size[datafr$Pop_type =="Samp2" ], breaks = 15) 
# shows frequency, breaks determine number of categories
# Note the number and position of breaks can change the outlook!
hist(datafr$Size[datafr$Pop_type =="Samp2" ], breaks = 5)
# deviation from normal distribution looks smaller 

hist(datafr$Size[datafr$Pop_type =="Samp2" ], breaks = 15, probability = TRUE, xlim =c(40, 140)) 
# shows probability densities
# we can add density lines of the distribution
dens <- density(datafr$Size[datafr$Pop_type =="Samp2" ])  
# estimates the density using a given smoother
lines(dens) 
# and of a normal distribution with the same mean and s
dens2 <- dnorm(0:140, mean = mean(datafr$Size[datafr$Pop_type =="Samp2" ]), sd = sd(datafr$Size[datafr$Pop_type =="Samp2" ]))
lines(dens2, col = "red") 
# the empirical cumulative distribution function would be calculated and plotted with
# ecdf_tot <- ecdf(datafr$Size[datafr$Pop_type =="Samp2" ])
# plot(ecdf_tot)

# provides mean, median and quartiles
summary(datafr$Size[datafr$Pop_type =="Samp2" ])

dev.off()

###########################
#     Extra Section 	  #
###########################

# Why are we not using a test to check for normal distribution?
# The reason is that for small sample sizes normality tests often do not reject 
# the null hypothesis, though it is wrong (H0: data originates from normal distribution)
# This example is taken from a blog: http://blog.fellstat.com/?p=61

# set seed for random number generator. This makes the example reproducible
set.seed(3300)
# We draw 15 observations from a binomial distribution with each of 5 trials and the probability of success
x <- rbinom(15, 5, 0.6) 
shapiro.test(x) 
# null hypothesis is not rejected
# Also the shapiro.test depends highly on sample size (rejecting H0 to often with big sample sizes)

qreference(x, nrep = 8) 
#deviation is quite obvious

###########################
#   End of Extra Section  #
###########################

########################
# Bootstrapping of CI  #
########################
library(boot)
dasboot <- boot(datafr[datafr$Pop_type=="Samp1", 1], function(x, i){mean(x[i])}, 10000,
                           parallel="multicore", ncpus=8)
plot(dasboot)
par(cex = 1.4)
hist(dasboot$t, breaks = 100, main = "", xlab = "t*")
boot.ci(dasboot, type = "bca")
# adjusted bootstrap percentile (BCa)
# add confidence intervals to plot



# directly implemented
library(coin)
pvalue(oneway_test(Size ~ Pop_type, data = datafr, distribution=approximate(B=9999)))

###################################
# Simple linear regression model  #
###################################



# now we load a file from the course website
read.table("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/possum.csv")
# does not look good

# we have to specify different options: 
# header = TRUE 
# sep (for separator) = ";" 
# dec (for decimal point) = "."
# We also omit the row.names.
# See ?read.table for the various options
pos_dat <- read.table(url("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/possum.csv"),  dec=".", sep=";", header=TRUE, row.names = NULL)
close(url("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/possum.csv"))

pos_dat 
# looks ok
head(pos_dat) 
# first 6 rows (for long data.frames)
str(pos_dat)  
# structure of object (data.frame with 14 variables and 104 rows, also class of variables)

###########################################
# file information
###########################################

# This file was taken from the DAAG package which is maintained by John Maindonald 
# and W. John Braun and contains information on possums that were published in: 
#
# Lindenmayer, D. B., Viggers, K. L., Cunningham, R. B., and Donnelly, C. F. 1995. 
# Morphological variation among columns of the mountain brushtail possum, 
# Trichosurus caninus Ogilby (Phalangeridae: Marsupiala). 
# Australian Journal of Zoology 43: 449-458
#
# You can find more details on the variables in the DAAD package (if installed):
# library(DAAG)
# data(possum)
# ?possum

# you should save the file for later use on your local hard drive
save(pos_dat, file = "Possum.Rdata")
# if you start a new session you can just load the saved data with:
# load("Possum.Rdata")

# Finally, you could conventionally download the file to your local harddrive and read it from there:
# read.csv2("/Users/ralfs/Downloads/Possum.csv") # would also work if adjusted to your path
# read.csv and read.csv2 are different functions for reading data frames 
# that have predefined options for csv files
# The most general function is read.table, were you have to define the format of your data
# (header, separator, decimal point, ...)

##########################
### Data handling
##########################




##########################
### Exploratory analysis #
##########################

###############################################
# Imagine that we would be interested whether #
# possums differ in their total length		  #
# between Victoria and other states		 	  #
###############################################

# we first inspect the total length variable
plot(pos_dat$totlngth) 
# overview

# polish plot
par(cex = 1.4)
# set parameters for following plot, in this case increase font size
# please check ?par() for explanation on the commands used
plot(pos_dat$totlngth, las = 1, cex.lab = 1.3, 
     ylab = "Total length [cm]", main = "Data overview") 
     # scatter plot

## Boxplot of Total Length
# can be used to screen for outliers
boxplot(pos_dat$totlngth, las = 1, cex.lab = 1.3, 
        ylab = "Total length [cm]", main = "Data overview")

# how can I save a graph? either use the graphical user interface or use jpeg, png, tiff...
jpeg("Boxplot.jpeg", quality = 100)
  par(cex = 1.4)
  boxplot(pos_dat$totlngth, las = 1, cex.lab = 1.3, 
        ylab = "Total length [cm]", main = "Data overview")
dev.off()
# switches off the device (in this case saves the contents to file in wd)

# although the boxplot is widely used there are interesting alternatives such as the beanplot
# for an introduction see the article see:
# 
# Peter Kampstra 2008: Beanplot: A Boxplot Alternative for Visual Comparison of Distributions. 
# Journal of Statistical Software, Code Snippets. 28 (1): 1-9. 
# Freely available at http://www.jstatsoft.org/v28/c01/
# 
# you need to install the package beanplot in order to use beanplots:
# install.packages("beanplot")
# load package

library(beanplot)
# use beanplot
beanplot(pos_dat$totlngth) 
# gives single observations, median and shows the estimated density function

# if we want to send the output to a new the graphics device we use:
quartz() 
# please note: on Windows OS use windows(), on Linux x11()

## Now we check for normal distribution and homogeneity of variances
# we consider homogeneity of variances first, because
# deviations are more relevant

#############################################
## 1. Checking for homogeneity of variance  #
#############################################

plot(totlngth ~ Pop, data = pos_dat)
# put on same median - code taken from plot.hov function
# this is only done to ease visual comparison
med <- tapply(pos_dat$totlngth, pos_dat$Pop, median)
# compute median of groups
w <- pos_dat$totlngth - med[pos_dat$Pop]
# compute distance to median and plot
# i.e. we move both data sets to 0 to ease comparison of the variance
plot(w ~ Pop, data = pos_dat)
# the figure displays the boxplots shifted to a median of zero
# variance looks a bit higher in Victoria, but is still acceptable

# checking for equality of variances with the F test
var.test(pos_dat$totlngth[pos_dat$Pop == "Vic"],  pos_dat$totlngth[pos_dat$Pop == "other"])
var.test(totlngth ~ Pop, data = pos_dat) # equivalent using formula notation
# both are equal, first one is using the formula notation
# null hypothesis is not rejected

# overall, we can assume that the variance is homogeneous
# in case that we conducted a t-test (or other test)
# and the p value would be close to our alpha
# we should compare the results to those of a methods
# that does not require homogeneity of variance

###################################################
# 2. check for normal distribution using qqplot   #
###################################################

# we need to check for each sample (i.e. from Victoria and other states)
# quantile-quantile plot
qqnorm(pos_dat$totlngth[pos_dat$Pop =="other" ], datax = TRUE)
# data fit the line relatively well, though some deviation is visible
qqline(pos_dat$totlngth[pos_dat$Pop =="other" ], datax = TRUE)
# a function qreference in the DAAG package that produces reference plots 
# to aid in the evaluation whether the data are normally distributed
library(DAAG)
qreference(pos_dat$totlngth[pos_dat$Pop =="other" ], nrep = 8) 
# nrep controls the number of reference plots
# data does not look conspicuous

# similarly, you could try to select your qqplot from a range of qq plots to see whether it deviates
# markably. The code is provided here:
# http://www.r-bloggers.com/checking-glm-model-assumptions-in-r/

qqnorm(pos_dat$totlngth[pos_dat$Pop =="Vic" ], datax = TRUE)
# data fit the line relatively well
qqline(pos_dat$totlngth[pos_dat$Pop =="Vic" ], datax = TRUE)
qreference(pos_dat$totlngth[pos_dat$Pop =="Vic" ], nrep = 8) 
# again some deviation is visible, but given that
# the t-test is relatively robust against violation 
# of this assumption, we don't need to worry

# Another option to check normality is the histogram
# but we will rather stick to qq plots in the future
# histograms can also be used to check for
# asymmetry in the distribution of variables

hist(pos_dat$totlngth[pos_dat$Pop =="other" ]) 
# shows frequency
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], probability = TRUE) 
# shows probability densities

# add density lines
dens <- density(pos_dat$totlngth[pos_dat$Pop =="other" ])  
# estimates the density using a given smoother
lines(dens) 
# the empirical cumulative distribution function would be calculated and plotted with
# ecdf_tot <- ecdf(pos_dat$totlngth[pos_dat$Pop =="other" ])
# plot(ecdf_tot)

# #Note that different breaks change the outlook!

xlim <- range(dens$x)
ylim <- range(dens$y) * 1.2 
# extract minimum and maximum values of density
# the plot suggests a normal distribution
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], breaks = 72.5 + (0:5) * 5, probability = TRUE, 
     xlim = xlim, ylim = ylim, 
     xlab = "Total length (cm)", main = "") 
lines(dens) 

# with other breaks it looks rather skewed
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], breaks = 75 + (0:5) * 5, probability = TRUE, 
     xlim = xlim, ylim = ylim, 
     xlab = "Total length (cm)", main = "") 
lines(dens) 

# effect of the number of breaks
par(mfrow = c(2,2))
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], breaks = 5, main = "5 breaks")
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], breaks = 10, main = "10 breaks")
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], breaks = 20, main = "20 breaks")
hist(pos_dat$totlngth[pos_dat$Pop =="other" ], breaks = 50, main = "50 breaks")

# gives a numeric representation of the median and the quartiles
summary(pos_dat$totlngth[pos_dat$Pop =="other" ])

dev.off()
# Why are we not using a test to check for normal distribution?
# The reason is that for small sample sizes normality tests often do not reject 
# the null hypothesis, though it is wrong (H0: data originates from normal distribution)
# This example is taken from a blog: http://blog.fellstat.com/?p=61

# set seed for random number generator. This makes the example reproducible
set.seed(3300)
# We draw 15 observations from a binomial distribution with each of 5 trials and the probability of success
x <- rbinom(15, 5, 0.6) 
shapiro.test(x) 
# null hypothesis is not rejected
# Also the shapiro.test depends highly on sample size (rejecting H0 to often with big sample sizes)

qreference(x, nrep = 8) 
#deviation is quite obvious

# after we have learnt how to check assumptions, you should be able to 
# complete the exercises below 

# # The following section gives the code for the graphs on the slides - 
# # if you want to understand a bite more of R, try to understand what is happening, 
# # also using the help pages. You could also modify some of the values and see what happens
# # 
# ## Scatterplot
# par(mfrow = c(1,3))
# par(cex = 1.8) # because mfrow resets cex
# plot(totlngth~hdlngth, data = pos_dat, las = 1, cex.lab = 1.4, 
     # ylab = "Total length [cm]", xlab = "Head length [mm]", main = "Scatterplot")

# ## boxplot
# temp <- boxplot(append(pos_dat$totlngth, 70), las = 1, cex.lab = 1.4 , xlim = c(0.8, 1.6),
        # ylab = "Total length [cm]", main = "Boxplot")
# # the function append() appends values to a vector

# # add annotations to plot
# text(1.23, temp$stats[5], "upper whisker", cex = 0.8, adj = 0)
# text(1.23, temp$stats[4], "up. quartile (75%)", cex = 0.8, adj = 0)
# text(1.23, temp$stats[3], "median", cex = 0.8, adj = 0)
# text(1.23, temp$stats[2], "low. quartile (25%)", cex = 0.8, adj = 0)
# text(1.23, temp$stats[1], "lower whisker", cex = 0.8, adj = 0)
# text(1.23, temp$out,"outlier (> 1.5 IQR)", cex = 0.8, adj = 0)

# ## Histgramm
# hist(pos_dat$totlngth, probability = TRUE, las = 1, xlim = c(71, 100), ylim = c(0, 0.12),
     # xlab = "Total length [cm]", main = "Histogram with density curve \n and normal distribution")
# dens <- density(pos_dat$totlngth)
# lines(dens, lwd = 2) 
# normdens <- dnorm(seq(70, 100, 0.1), mean = mean(pos_dat$totlngth), sd=sd(pos_dat$totlngth))
# lines(seq(70, 100, 0.1), normdens, lwd = 2, col = "red")

################################################################################
# Exercise 1: Check whether male and female possums have the same total length #
################################################################################

## Checking for homogeneity of variance
plot(totlngth ~ sex, data = pos_dat)
# put on same median - code taken from plot.hov function
med <- tapply(pos_dat$totlngth, pos_dat$sex, median)
# compute median of groups
w <- pos_dat$totlngth - med[pos_dat$sex]
# compute distance to median and plot
plot(w ~ sex, data = pos_dat)
# looks very similar

# checking for equality of variances with the F test
var.test(totlngth ~ sex, data = pos_dat)
# null hypothesis is not rejected

# check for normal distribution using qqplot
# quantile quantile plot 
qqnorm(pos_dat$totlngth[pos_dat$sex == "m"], datax = TRUE)
qqline(pos_dat$totlngth[pos_dat$sex == "m"], datax = TRUE)
qreference(pos_dat$totlngth[pos_dat$sex == "m"], nrep = 8) 
# nothing to worry
qqnorm(pos_dat$totlngth[pos_dat$sex == "f"], datax = TRUE)
qqline(pos_dat$totlngth[pos_dat$sex == "f"], datax = TRUE)
qreference(pos_dat$totlngth[pos_dat$sex == "f"], nrep = 8) 
# looks also ok

# we assume normal distribution and homogeneity of variance and use the t-test
t.test(pos_dat$totlngth[pos_dat$sex == "m"],  pos_dat$totlngth[pos_dat$sex == "f"], var.equal = TRUE)
# t.test(totlngth ~ sex, data = pos_dat, var.equal = TRUE)
# equivalent using formula notation

######################
# New section: ANOVA #
######################

# This part of the script generates the plots for the presentation
# Run this part of the sript and try to understand the idea behind these plots - 
# use the R help to get more information on the code snippets

# a <- sample(1:50, 20, replace = TRUE)
# b <- sample(60:70, 20, replace = TRUE)
# c <- sample(80:150, 20, replace = TRUE)
# new_rand <- c(a, b, c)
# dis_dat <- data.frame(Value = new_rand, 
#                       Group = c(rep("a",20),rep("b",20),rep("c",20)))
# order_vec <- sample(1:60 ,60)
# dis_dat_rn <- dis_dat[order_vec, ]

# png("SSY.png", bg = "transparent")
#   par(cex = 1.5)
#   plot(dis_dat_rn[ , 1], ylab = "y", xlab = "order", las = 1)
#   abline(mean(dis_dat_rn[ , 1]), 0, lwd = 2)
#   for(i in 1:nrow(dis_dat_rn)) {
#     lines(c(i, i), c(mean(dis_dat_rn[ , 1]), dis_dat_rn[i, 1]), lwd = 1.5)
#   }
#   mtext(expression(SSY == sum((y-bar(y)))^2), side = 3, line = 1.5, cex = 2)
# dev.off()

# png("SSE.png", bg="transparent")
  # par(cex = 1.5)
  # plot(dis_dat_rn[ , 1], type = "n", ylab = "y", xlab = "order", las = 1)

  # points(seq(1, 60, 3), dis_dat_rn[dis_dat_rn[ , 2] == "a", 1], pch = 16)
  # points(seq(2, 60, 3), dis_dat_rn[dis_dat_rn[ , 2] == "b", 1], pch = 18 )
  # points(seq(3, 60, 3), dis_dat_rn[dis_dat_rn[ , 2] == "c", 1], pch = 21)

  # abline(mean(dis_dat_rn[dis_dat_rn[ , 2] == "a", 1]), 0, lwd = 2)
  # abline(mean(dis_dat_rn[dis_dat_rn[ , 2] == "b", 1]), 0, lwd = 2, lty = 2)
  # abline(mean(dis_dat_rn[dis_dat_rn[ , 2] == "c", 1]), 0, lwd = 2, lty = 4)
  
  # k <- -2
  # for (i in 1:20){
    # k <- k+3
    # lines(c(k, k), c(mean(dis_dat_rn[dis_dat_rn[ , 2] == "a", 1]), dis_dat_rn[dis_dat_rn[ , 2] == "a", 1][i]), lwd = 2)
    # lines(c(k+1, k+1), c(mean(dis_dat_rn[dis_dat_rn[ , 2]== "b", 1]), dis_dat_rn[dis_dat_rn[ , 2] == "b", 1][i]), lwd = 2)
    # lines(c(k+2, k+2), c(mean(dis_dat_rn[dis_dat_rn[ , 2]== "c", 1]), dis_dat_rn[dis_dat_rn[ , 2] == "c", 1][i]), lwd = 2)
  # }
  # mtext(expression(SSE == sum(sum((y-bar(y[j])))^2)), side = 3, line = 1.5, cex = 2)
# dev.off()
# end of plot


####################################################################################
# Exercise 2: Divide the possums into 4 age groups of approximately similar sample size 
# and examine whether the head lengths differ
####################################################################################

# see linear regression model below for diagnostics

## Solution
# Use the cut function to create 4 groups of similar sample size and 
# add this as column to our data.frame
pos_dat$age_groups <- cut(pos_dat$age, c(0, 2.25, 3, 5, 10))
head(pos_dat)

# put on same median to ease the checking of homogeneity of variance
med <- tapply(pos_dat$hdlngth, pos_dat$age_groups, median)
# compute median of groups
w <- pos_dat$hdlngth - med[pos_dat$age_groups]
# compute distance to median and plot
plot(w ~ age_groups, data = pos_dat)

# variation seems to decrease with age

# however, we can fit the anova and run model diagnostics afterwards to decide
# whether the assumptions of normality and homogeneity of variance are met

## conduct anova
summary(aov(hdlngth ~ age_groups, data = pos_dat))
# displays the anova table for the response variable (value) as explained by the factor (Group)
# p is slightly smaller than 0.05, and would be considered statistically significant

# this only tests the null hypothesis - in most cases the effects size for each factor is more interesting
summary.lm(model_1 <- aov(hdlngth ~ age_groups, data = pos_dat))
# note that the hypothesis tests displayed for summary.lm are treatment contrasts
# this means that the difference of each level to the first level is tested.
# you could set other contrasts, see the section on contrasts in Crawleys R Book

# how to interpret the parameters in the table?
mean_a <- mean(pos_dat$hdlngth[pos_dat$age_groups == "(0,2.25]"], na.rm = TRUE)
# this is the intercept
mean_a
mean_b <- mean(pos_dat$hdlngth[pos_dat$age_groups == "(2.25,3]"], na.rm = TRUE)
mean_b
mean_b - mean_a 
# this gives the second line in the summary table: the difference between mean for a and mean for b

# We could also remove the intercept from our model (by adding -1), then we have the estimates (means) for every group
summary.lm(model_1 <- aov(hdlngth ~ age_groups - 1, data = pos_dat))
# however, this makes not much sense, as then the model compares the estimate to 0
# note that this inflates the explained variance 

# now we turn to model diagnostics
par(mfrow = c(1, 2), cex = 1.5)
plot(model_1, which = 1:2)
# plot 1 allows checking for homoscedasticity (equal variances of residuals)
# plot 2 allows checking for normal distribution of the residuals.
# the homogeneity of variances can also be examined using hypothesis testing such 
# as Levenes Test (leveneTest {car}) or Bartlett`s test (bartlett.test), 
# but we will not do this for reasons discussed elsewhere (see slides and notes)
# we take a look at the reference plots to evaluate whether the residuals are normally distributed

library(DAAG)
quartz() # please note: on Windows OS use windows(), on Linux x11()
qreference(residuals(model_1), nrep = 4)
# nrep controls the number of reference plots
# there may be a slight departure from normality due to a few observations 

# given that variances do not look equal, we consider Welch's anova
oneway.test(hdlngth ~ age_groups, data = pos_dat, var.equal = FALSE)
# leads to a slightly higher p value (0.054 instead of 0.041)
# this example shows that it is questionable to 
# employ a binary decision such as significant or non-significant
# we clearly have a borderline case here
# In a paper, you would rather report 0.054 to be on the safe side
# and look at effect sizes

###################################
#Demonstration linear regression  #
###################################

library(vegan)
data("varechem")
# we use some soil data to demonstrate how to conduct a linear regression analysis
?varechem
head(varechem)
str(varechem)

# we check the relationship between K and S
# linear regression for data
mod_1 <-lm(S ~ K, data = varechem) 
# regression for air pollution by factories (as expressed with the ~)
summary.lm(mod_1)
# this function gives an overview on residuals and coefficients. 

# plot relationship
plot(S ~ K, data = varechem)
# and regression line
par(las = 1,cex = 1.8)
plot(S ~ K, data = varechem, ylab="S concentration", xlab="K concentration")
# we extract the lines from the regression model with the function abline()
abline(mod_1, lwd = 2)

##for publication style output use
library(memisc)
mtable(mod_1)
# there are various arguments that can be adjusted for this function

### model checking
# can also be used for anova
par(mfrow = c(1, 3))
plot(mod_1, which = 1:3)
# the plots 1 and 3 can be used to spot departure from the assumption of homoscedasticity
# plot 2 is the quantile-quantile plot and serves as a diagnosis tool to check for normality of residuals, 
# plot 3 is similar to plot 1 but uses standardized residuals. 
# Standardized residuals are obtained by dividing each residual by its standard deviation. 
# They show (1) how many standard deviations any point is away from the fitted regression model and 
# (2) have constant variance even if residuals have nonconstant variance (due to leverage). 
# In the absence of leverage points plot 1 and 3 are very similar. 
# Standardized residuals can also be used to check for outliers. 
# Points that are more than 2 standard deviations away from the regression line may 
# be considered as outlier (see Sheater 2009: p.60).

# for more exploration of residuals, see:
https://www.r-bloggers.com/visualising-residuals/

# if you feel insecure regarding the evaluation of deviation from normality, 
# the qreference function is helpful and allows comparison with data drawn from a normal distribution
qreference(residuals(mod_1), nrep = 8)
# nothing to worry

# Plots on influence of specific points
par(mfrow = c(1,3))
plot(mod_1, which = 4:6)
# these 3 plots display deliver information on influential observations.
# the second plot is usually the first to consider
# if any observation is inside the boundaries of Cooks distance
# we don't need to worry about influential points 

# calculation of two times the average hat value
# 2* p/n - see lecture
2*2/nrow(varechem) 
# gives double the average hat value and 
# plot 2 and 3 allow checking for influential points with a hat value larger than this average hat value. 

dfbeta(mod_1) 
# gives the the change in the model parameter when the respective point is omitted
# leave-one-out regression
# example
mod_2 <-lm(S ~ K, data = varechem[-1, ]) 
summary(mod_1)
summary(mod_2)
dfbeta(mod_1)[ 1,]

###########################
#   some helper functions #
###########################
### Extracting results from model
coefficients(mod_1)
# gives the coefficients of the model

fitted(mod_1)
# fitted values (predictions for the x)

confint(mod_1)
# confidence intervals for the parameters

residuals(mod_1)
# model residuals

rstandard(mod_1)
# standardized residuals

anova(mod_1)
# gives the analysis of variance table with the explained variance for the variables 
# and the unexplained variance as residuals

str(mod_1) 
# structure of fitted model object, shows what can be extracted from model

# plots for slides
# though this section is only relevant for creating the slides, you may want to check 
# the commands in case you are interested in graphical plotting
# 
# We load some new data that gives the biological growth of a species in dependence of tannin. 
# This data was provided by Crawley (2015), see pages 116 following.
# tan_dat <- read.table("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/tannin.txt", header=TRUE)
# summary.lm(th_mod <- lm(growth ~ tannin, data = tan_dat))#
# png("SSE_lm.png", bg = "transparent")
  # par(cex = 1.5, las = 1)
  # plot(growth ~ tannin, data = tan_dat, cex.main = 1.3,
       # xlab = "Var 1", ylab = "Var 2",        
       # main = c(expression(SSE == sum(italic((y[i]-b[0]-b[1]*x[i]))^2,i = 1,n))))
  # abline(th_mod, lwd = 1.5) # fits regression line, uses coefficients from object th_mod
  # fitted <- predict(th_mod) # prediction of y values based on regression equation
  # for (i in 1:9){ 
    # lines(c(tan_dat$tannin[i], tan_dat$tannin[i]), c(tan_dat$growth[i], fitted[i]), lwd = 1.5)
  # }
# dev.off()

# png("SSY_lm.png", bg="transparent")
  # par(cex = 1.5, las = 1)
  # plot(growth ~ tannin, data = tan_dat, cex.main = 1.3,
      # xlab = "Var 1", ylab = "Var 2", main = c(expression(SSE == sum((y-symbol(a)-symbol(b)*x)^2))))
  # abline(mean(tan_dat$growth), 0, lwd = 2)
  # for(i in 1:length(tan_dat$growth)) {
    # lines(c(tan_dat$tannin[i], tan_dat$tannin[i]), c(mean(tan_dat$growth), tan_dat$growth[i]), lwd = 1.5)
  # }
# dev.off()

# # additional figures
# growth2 <- c(tan_dat$growth, 0) # point added that represents a leverage point
# tannin2 <- c(tan_dat$tannin, 20) # point added that represents a leverage point
# th_mod2 <- lm(growth2 ~ tannin2)
# influence(th_mod2)$hat
# summary.lm(th_mod2)
# png("Hat_val_1.png")
#   par(las = 1, cex = 1.5)
#   plot(tannin2, growth2, main = "Hat value for x=20 is 0.81", xlab = "x", ylab = "y")
#   abline(th_mod2, lwd = 2)
# dev.off()
# 
# # and second graph...
# growth3 <- c(tan_dat$growth, -12) # point added that represents a leverage point
# tannin3 <- c(tan_dat$tannin, 20) #
# th_mod3 <- lm(growth3~tannin3)
# influence(th_mod3)$hat
# summary.lm(th_mod3)
# png("Hat_val_2.png")
#   par(las = 1, cex = 1.5)
#   plot(tannin3, growth3, main = "Hat value for x=20 is 0.81", xlab = "x", ylab = "y")
#   abline(th_mod3, lwd = 2)
# dev.off()

################################################################################################
# Exercise 3: Develop regression models for total length of possums predicted 
# by foot length and skull width as predicted by head length. 
# Would you recommend to predict the total length of the possum from traces of their feet in the snow? 
# And is an invasive measurement of the skull width necessary or can it be approximated with the head length?
################################################################################################

## Look at the data
# relationship between foot length and possums.
plot(totlngth~footlgth, data = pos_dat, 
     ylab = "Total length", xlab = "Foot length")
# slightly positive relationship between both variables

## model fitting
foot_mod <- lm(totlngth ~ footlgth, data = pos_dat)
summary.lm(foot_mod)
# only 20% of variation explained - model is not very good for predicting the total length of the possum

### model checking
par(mfrow = c(2, 2))
plot(foot_mod)
# the plots 1 and 3 can be used to spot departure from the assumption of homoscedasticity. 
# Both figures do indicate that heteroscedasticity is not a problem here
# In addition, plot 3 can be used to spot outliers. 
# In fact, observations 39 and 44 may be considered as outliers as they deviate > 2 sd from the fitted model. 
# Altough plot 4 confirms that both observations (39 and 44) are outliers, 
# it also suggests that they are not influential in terms of Cooks distance. 
# Plot 2 is the quantile-quantile plot and serves as a diagnosis tool to check for normality of residuals. 
# We use reference plots from a normal distribution to ease interpretation:
quartz()
# please note: on Windows OS use windows(), on Linux x11()
qreference(residuals(foot_mod), nrep = 8)
# no indication of non-normal distribution

##################################################################################################
# although not really necessary for this case, for purpose of demonstration we fit the 
# model again without both outliers to see whether the parameter estimates and model fit changes
foot_mod2 <-lm(totlngth ~ footlgth, data = pos_dat[-c(39, 44), ])
summary.lm(foot_mod2)
# indeed a slight increase in model fit!
par(mfrow = c(2, 2))
plot(foot_mod2)
# nothing to worry in model diagnostics
# let`s plot the model and see to which extent the outliers influence the model parameters
quartz()
par(las = 1,cex = 1.8)
plot(totlngth ~ footlgth, data = pos_dat, ylab="Total length", xlab="Foot length")
# we extract the lines from the regression model with the function abline()
abline(foot_mod, lwd = 2)
abline(foot_mod2, lwd = 2, lty = 2, col = "red")
# only slight change without outliers
# finally we mark the outliers
points(pos_dat$footlgth[c(39, 44)], pos_dat$totlngth[c(39,44)], 
       pch = 16, col = "blue")
# for a scientific publication you would also report the model results with and 
# without the outliers and mark the outliers in a similar manner (but not necessarily in blue :-)

