#################################################################################
# Script for Session 1 of the course: Applied multivariate Statistics with R		#
# 					by Ralf B. Schäfer											#
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
par(mfrow = c(1, 2), cex = 1.2)
# set parameters for following plot, in this case increase font size
# with cex and plot two figures beside each other with mfrow
# please check ?par() for explanation on the commands used
hist(pop1)
hist(pop2)

# we sample with a relatively small size from the populations
# this translates, for example, to a monitoring study
# with biometric measurements of a sample of the population
no_sam1 <- 30
no_sam2 <- 30
set.seed(50)
samp1 <- sample(pop1, no_sam1)
mean(samp1)
set.seed(30)
samp2 <- sample(pop2, no_sam2)
mean(samp2)
par(mfrow = c(1, 2))
hist(samp1, breaks = 15)
hist(samp2, breaks = 15)
# we create a categorical variable that codes for the two populations
categ <- c(rep("Samp1", no_sam1), rep("Samp2", no_sam2))
# and combine the two populations into one vector
data_vec <- c(samp1, samp2)
# create data.frame
datafr <- data.frame(data_vec, categ)
names(datafr) <- c("Size", "Pop_type")

# Imagine we have sampled from the two populations
# and now want to compare their body size
# we have the hypothesis that the body sizes differ between the two populations

# we create a function
meanDif <- function(x, grp) {
  mean(x[grp == "Samp1"]) - mean(x[grp == "Samp2"])
}
with(datafr, meanDif(Size, Pop_type))

# we use permutation for hypothesis testing
library(permute)
perm_vec <- numeric(length = 10000)
# Legendre & Legendre (2012) suggest to use at least 10,000 permutations
# for inference, especially when aiming at publishing the results
N <- nrow(datafr)
set.seed(400)
for (i in seq_len(length(perm_vec) - 1)) # loop runs 9999 times
{
  perm <- shuffle(N)
  perm_vec[i] <- with(datafr, meanDif(Size, Pop_type[perm]))
}
# results in 9999 shuffled data sets
# add non-shuffled data	because we test the null hypothesis
# i.e. assume that the unpermuted statistic could be obtained under H0
perm_vec[10000] <- with(datafr, meanDif(Size, Pop_type))

# plot the results of the permutations
par(cex = 1.2, mfrow = c(1, 1))
hist(perm_vec, main = "", breaks = 20, xlab = expression("Mean difference (Pop1 - Pop2)"))
rug(perm_vec[10000], col = "red", lwd = 2, ticksize = 0.5)

# compute permutation p-value for one sided test
Dbig <- sum(perm_vec >= perm_vec[10000])
Dbig / length(perm_vec)
# permutation test identifies statistical significant difference
# but only for one-sided test!
# for two-sided test
Dbig2 <- sum(abs(perm_vec) >= perm_vec[10000])
Dbig2 / length(perm_vec)
# two sided test not significant
# note that you should more or less always use a two-sided test
# see: Lombardi C.M. & Hurlbert S.H. (2009)
# Misprescription and misuse of one-tailed tests. Austral Ecology 34, 447–468.

# direct implementation of the permutational t-test
library(coin)
pvalue(oneway_test(Size ~ Pop_type, data = datafr, distribution = approximate(B = 9999)))

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

# if we want to send the output to a new the graphics device we use:
quartz()
# please note: on Windows OS use windows(), on Linux x11()

# although the boxplot is widely used there are interesting alternatives
# such as the beanplot
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
beanplot(w ~ Pop_type, data = datafr)
# gives single observations (white lines), median
# and shows the estimated density function

# although we know that the assumption of normality is violated,
# we run a t-test anyways, accounting for different variances
(body_t2 <- t.test(Size ~ Pop_type, data = datafr, var.equal = FALSE))
# no statistical significance, p-value slightly higher than for permutation

# how would we check for normal distribution?

##########################################################
# QQplots and boxplots to check for normal distribution  #
##########################################################

# we need to check for each sample
# quantile-quantile plot
qqnorm(datafr$Size[datafr$Pop_type == "Samp1" ], datax = TRUE)
# data fit a line relatively well, though some deviation is visible
qqline(datafr$Size[datafr$Pop_type == "Samp1" ], datax = TRUE)
# the function qreference in the DAAG package produces reference plots
# to aid in the evaluation whether the data are normally distributed
library(DAAG)
qreference(datafr$Size[datafr$Pop_type == "Samp1" ], nrep = 8)
# nrep controls the number of reference plots
# data does look perfectly normal

# similarly, you could try to select your qqplot from
# a range of qq plots to see whether it deviates markably.
# The code is provided here:
# http://www.r-bloggers.com/checking-glm-model-assumptions-in-r/

qqnorm(datafr$Size[datafr$Pop_type == "Samp2" ], datax = TRUE)
qqline(datafr$Size[datafr$Pop_type == "Samp2" ], datax = TRUE)
# stronger deviation from line
qreference(datafr$Size[datafr$Pop_type == "Samp2" ], nrep = 8)
# deviation is conspicuous compared to samples from normal distribution

# Another option to check normality is the histogram
# histograms can also be used to check for
# asymmetry in the distribution of variables (i.e. skewness)

hist(datafr$Size[datafr$Pop_type == "Samp2" ], breaks = 15)
# shows frequency, breaks determine number of categories
# Note the number and position of breaks can change the outlook!
hist(datafr$Size[datafr$Pop_type == "Samp2" ], breaks = 5)
# deviation from normal distribution looks smaller

hist(datafr$Size[datafr$Pop_type == "Samp2" ], breaks = 15, probability = TRUE, xlim = c(40, 140))
# shows probability densities
# we can add density lines of the distribution
dens <- density(datafr$Size[datafr$Pop_type == "Samp2" ])
# estimates the density using a given smoother
lines(dens)
# and of a normal distribution with the same mean and standard deviation
dens2 <- dnorm(0:140, mean = mean(datafr$Size[datafr$Pop_type == "Samp2" ]), sd = sd(datafr$Size[datafr$Pop_type == "Samp2" ]))
lines(dens2, col = "red")
# the empirical cumulative distribution function would be calculated and plotted with
# ecdf_tot <- ecdf(datafr$Size[datafr$Pop_type =="Samp2" ])
# plot(ecdf_tot)

# provides mean, median and quartiles
summary(datafr$Size[datafr$Pop_type == "Samp2" ])
# close graphics window
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
# We draw 15 observations from a binomial distribution
# with each of 5 trials and the probability of success
x <- rbinom(15, 5, 0.6)
shapiro.test(x)
# null hypothesis is not rejected
# Also the shapiro.test depends highly on sample size
# (rejecting H0 to often with big sample sizes)

qreference(x, nrep = 8)
# deviation is quite obvious

###########################
#   End of Extra Section  #
###########################

########################
# Bootstrapping of CI  #
########################

# finally, we boostrap confidence intervals for the
# first sample population
library(boot)
boot_samp1 <- boot(
  datafr[datafr$Pop_type == "Samp1", 1],
  function(x, i) {
    mean(x[i])
  },
  10000, parallel = "multicore", ncpus = 8
)
# on Windows set parallel ="no" and remove 'ncpus = 8'
# plot distribution of bootstrapped statistic and qqplot
plot(boot_samp1)
# normally distributed

# frequency plot
par(cex = 1.4)
hist(boot_samp1$t, breaks = 100, main = "", xlab = "t*")
# compute adjusted bootstrap confidence interval
(cis <- boot.ci(boot_samp1, type = "bca"))
# add to plot
lines(x = c(cis$bca[4], cis$bca[5]), y = c(300, 300), col = "red", lwd = 2)
mtext(text = c("95% Confidence interval"), side = 3, col = "red", cex = 1.2)

###########################################
# Demonstration Simple linear regression  #
###########################################

library(vegan)
data("varechem")
# we use some soil data to demonstrate how to conduct a linear regression analysis
?varechem
head(varechem)
str(varechem)

# we check the relationship between K and S
# linear regression for data
mod_1 <- lm(S ~ K, data = varechem)
# regression for S explained by K (as expressed with the ~)
summary.lm(mod_1)
# this function gives an overview on residuals and coefficients.

# plot relationship
plot(S ~ K, data = varechem)
# and regression line
par(las = 1, cex = 1.8)
plot(S ~ K, data = varechem, ylab = "S mg/kg", xlab = "K mg/kg")
# we extract the lines from the regression model with the function abline()
abline(mod_1, lwd = 2)

## for publication style output use
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
# https://www.r-bloggers.com/visualising-residuals/

# if you feel insecure regarding the evaluation of deviation from normality,
# the qreference function is helpful and allows comparison with data drawn from a normal distribution
qreference(residuals(mod_1), nrep = 8)
# nothing to worry

# Plots on influence of specific points
par(mfrow = c(1, 3))
plot(mod_1, which = 4:6)
# these 3 plots display deliver information on influential observations.
# the second plot is usually the first to consider
# if any observation is inside the boundaries of Cooks distance
# we don't need to worry about influential points

# calculation of two times the average hat value
# 2* p/n - see lecture
2 * 2 / nrow(varechem)
# gives double the average hat value and
# plot 2 and 3 allow checking for influential points with a hat value larger than this average hat value.

dfbeta(mod_1)
# gives the the change in the model parameter when the respective point is omitted
# leave-one-out regression
# example
mod_2 <- lm(S ~ K, data = varechem[-1, ])
summary(mod_1)
summary(mod_2)
dfbeta(mod_1)[ 1, ]
# dfbeta provides the difference in the coefficients
# between the models 1 and 2

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

###################################
# Bootstrapping regression models #
###################################

library(car)
# we use a convenient function to produce the confidence interval
# first for the residual method (fixed x resampling)
set.seed(30)
boot_mod_res <- Boot(mod_1, f = coef, R = 9999, method = c("residual"))
confint(boot_mod_res, level = .95, type = "bca")
confint(mod_1)
# bootstrapped confidence intervals for the parameters
# differ slightly for intercept, can be regarded as more exact

# bootstrapping cases (random x resampling)
set.seed(30)
boot_mod_case <- Boot(mod_1, f = coef, R = 9999, method = c("case"))
confint(boot_mod_case, level = .95, type = "bca")
confint(boot_mod_res, level = .95, type = "bca")
confint(mod_1)
# bootstrapping cases provides most narrow residuals

# see Fox (2015):647ff and an online Appendix to Fox & Weisberger (2010): https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Bootstrapping.pdf for further details

######################################
# Cross-validating regression models #
######################################

# we compute the cross-validated Mean square prediction error
# that can be used to compare models, may not be most interesting for
# the simple linear regression model
cv.lm(data = varechem, form.lm = formula(S ~ K), m = 5, seed = 30, plotit = "Observed")
# compare to mean square error from fitted model
mean(mod_1$residuals ^ 2)
# prediction error (not surprisingly) higher than model error

# we can also compute a cross-validated version of the R2
# function provided by Kabacoff 2011:214
# R2 is validated using cross validation
# execute all code below to set the function
# adjust k for k-fold cross-validation
shrinkage <- function(fit, k = 5) {
  library(bootstrap)
  theta.fit <- function(x, y) {
    lsfit(x, y)
  }
  theta.predict <- function(fit, x) {
    cbind(1, x) %*% fit$coef
  }
  x <- fit$model[, 2:ncol(fit$model)]
  y <- fit$model[, 1]
  results <- crossval(x, y, theta.fit, theta.predict, ngroup = k)
  r2 <- cor(y, fit$fitted.values) ^ 2
  r2cv <- cor(y, results$cv.fit) ^ 2
  cat("Original R-square =", r2, "\n")
  cat(k, "Fold Cross-Validated R-square =", r2cv, "\n")
  cat("Change =", r2 - r2cv, "\n")
}
#### end of function

# calculate for model 1
shrinkage(mod_1)
# r2 shrinked by about 5%

############################
# Load data for Exercises  #
############################

# now we load a file from the course website
read.table("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/possum.csv")
# does not look good

# we have to specify different options:
# header = TRUE
# sep (for separator) = ";"
# dec (for decimal point) = "."
# We also omit the row.names.
# See ?read.table for the various options
pos_dat <- read.table(url("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/possum.csv"), dec = ".", sep = ";", header = TRUE, row.names = NULL)
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
# read.csv2("/Users/ralfs/Downloads/Possum.csv")
# would also work if adjusted to your path
# read.csv and read.csv2 are different functions for reading data frames
# that have predefined options for csv files
# The most general function is read.table, were you have to define the format of your data
# (header, separator, decimal point, ...)

###############################################################
# Exercise: Imagine that we would be interested whether       #
# male and female possums differ in their total length		  #
# Choose a statistical approach that allows to      			  #
# answer this question and compare the results for the 		  #
# classical and simulation-based approach  					  #
# Inspect the data for outliers (and other issues you 	      #
# consider necessary) before analysis using exploratory tools #
###############################################################


#################################################################
# Exercise: Chasing possums can be laborious.   					#
# An easy way would be to predict their length from footprints.	#
# Would you recommend to predict the total length of the possum #
# from traces of their 	feet in the snow?   						#
# And is an invasive measurement of the skull width necessary   #
# or can it be approximated with the head length?				#
# Choose a statistical approach that allows to      			#
# answer this question and compare the results for the 		  	#
# classical and simulation-based approach  					  	#
# Conduct complete model diagnostics							#
#################################################################





####################################################################################
# # 						Code for graphs on slides
# # if you want to understand a bite more of R, try to understand what is happening,
# # also using the help pages. You could also modify some of the values and see what happens
# #
# ## Scatterplot
# par(mfrow = c(1,3))
# par(cex = 1.8) # because mfrow resets cex
# plot(totlngth ~ hdlngth, data = pos_dat, las = 1, cex.lab = 1.4,
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
##############################################################################
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
