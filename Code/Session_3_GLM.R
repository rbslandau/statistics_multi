#################################################################################
# Script for Session 3 of the course: Applied multivariate Statistics with R	#
# 						by Ralf B. Schäfer											#
# 						WS 2017/18														#
# 		The script describes Generalised linear models in R						#
#################################################################################

# let us first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine here
#
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory.
# Subsequently, copy the path the the folder without (!!!) the file
# reference into the setwd() function

## Demonstration - code modified from Maindonald 2010, chapter 8
# Research question: Which environmental variables are most relevant
# to explain the occurrence of the Southern Corroboree frog,
# a critically endangered species

## multiple regression for GLM
library(DAAG)
# display information on data set that will be used
?frogs

### analysis with lm shows violation of assumption
flaw_mod <- lm(frogs$pres.abs ~ frogs$distance)
summary(flaw_mod)
par(cex = 1.5, mfrow = c(2, 2))
plot(flaw_mod)

#############################################################
# Checking distribution of variables and transformation     #
#############################################################

# We have a look at the data, since the data is from a restricted area,
# geographical parameters are omitted from analysis

library(car)
spm(frogs[, c(4:10)])
# Distance, NoOfSites and NoOfpools look left skewed.

# we also check the data range
summary(frogs$distance)
summary(frogs$NoOfPools)
summary(frogs$NoOfSites)
# NoOfSites has a relatively narrow data range
# does not require transformation
# The two other variables span 2 orders of magnitude
# Transformation should be considered

# we check whether sqrt or log transformation improves the distribution
par(mfrow = c(2, 3))
for (nam in c("distance", "NoOfPools")) {
  y <- frogs[, nam]
  plot(density(y), main = "", xlab = nam)
  plot(density(sqrt(y)), main = "", xlab = nam)
  plot(density(log(y)), main = "", xlab = nam)
}

# log transformation selected for NoOfPools and distance,
frogs$logNoPools <- log(frogs$NoOfPools)
frogs$logdistance <- log(frogs$distance)

##############################
# Checking for collinearity  #
##############################

# we inspect the data for collinearity
# columns 5 and 6 are the untransformed variables
# which are not selected
spm(frogs[, c(4, 7:12)])
# very high correlations between meanmax, altitude and meanmin
# this function represent an alternative
# to the function introduced in the context of
# multiple linear regression

dev.off()

# check correlation coefficients
cor(frogs[, c(4, 7:12)])
# almost correlation of 1 for altitude, meanmax and meanmin

# compute the VIF
frog.glm <- glm(
  pres.abs ~ logdistance + logNoPools + NoOfSites + avrain
    + altitude + meanmax + meanmin, family = binomial(link = "logit"),
  data = frogs, na.action = "na.fail"
)
# as for the lm, we first fit the glm. Syntax is similar to that of
# the lm, except for the specification of the binomial distribution
# and logit link. See ?family for other distributions
vif(frog.glm)
# very high VIFs

# we omit altitude from the model, combine meanmin and meanmax and recheck VIF
frog.glm2 <- glm(
  pres.abs ~ logdistance + logNoPools + NoOfSites
    + avrain + I(meanmax + meanmin), family = binomial,
  data = frogs, na.action = "na.fail"
)
# combination of meanmax and meanmin via I
# is equivalent to including the mean of both variables
# this represents an alternative to omission of one of the variables

vif(frog.glm2)
# avrain and combination of meanmax and meanmin now below 4

#######################################################
# Computation of all possible models  				  #
#######################################################
# As for the linear model, we first compute
# compute all possible models
# in fact, application of this technique as well as the
# following ones is largely identical code-wise as for the lm

library(MuMIn)
# computation of all models
models <- dredge(frog.glm2)
models
# for delta AIC of 2 only 1 best fit model

### model averaging
cand_models <- get.models(models, subset = delta < 3)
# here a delta AIC of 3 is chosen
# as otherwise only 1 model is selected

summary(model.avg(cand_models))
# except for NoOfSites no difference between models

##################################################################
# stepwise backward model selection using hypothesis testing on variables
##################################################################

# we check the initial model fit using Wald z-test
summary(frog.glm2)
# Based on z value (normal distribution), which is
# related to the test of the parameter differing from 0,
# we should remove NoOfSites

frog.glm3 <- update(frog.glm2, ~ . - NoOfSites)
summary(frog.glm3)
# now all parameters significant

#################################################################
# stepwise backward model selection with Likelihood ratio test 	#
#################################################################

# instead of checking for individual coefficients, we could also
# test for significant differences in the likelihood between models
# this should generally be preferred, though both tests often
# lead to similar inferences

# as for the lm, we can use the anova test and specify
# the chisquare test for comparing the deviances
# technically, this is a Likelihood ratio test
anova(frog.glm3, frog.glm2, test = "Chisq")
# no significant difference, hence this test supports simplification
# you should use test = "F" for quasi-likelihood based models

# using anova for a single model
# gives sequential reduction in deviance,
anova(frog.glm3, test = "Chisq")
# when terms are added sequentially (Type I ANOVA)

# makes not much sense to draw inference from this, as the following example shows:
frog.glm4 <- glm(
  pres.abs ~  I(meanmax + meanmin) + avrain +
    logdistance + logNoPools, family = binomial,
  data = frogs, na.action = "na.fail"
)
anova(frog.glm4, test = "Chisq")
# now all terms are significant, only by changing order in formula!

# We should rather us Type II ANOVA
Anova(frog.glm4)
# is the same as checking whether one term
# can be dropped without reducing the model fit significantly:
drop1(frog.glm4, test = "Chisq")
# all deletions result in a significant reduction in deviance

###################################################################
# Stepwise backward selection using information theoretic approach#
###################################################################

# if using information theoretic approaches,
# the functions are again the same as for the linear model
# calculation of AIC
AIC(frog.glm2)
AIC(frog.glm3)
# AIC lower, also suggesting the removal of NoOfSites

# calculation of BIC
BIC(frog.glm2)
BIC(frog.glm3)
# BIC lower, also suggesting the removal of NoOfSites

##################################################
# automatic model building - stepwise modelling  #
##################################################

# we have set the full model above, we specify the null model
null <- glm(pres.abs ~ 1, family = binomial, data = frogs, na.action = "na.fail")
# nullmodel: no variables, only mean

# calculate sample size for BIC
n <- length(frogs$pres.abs)
step(
  frog.glm2, direction = "both", trace = 100,
  scope = list(upper = frog.glm2, lower = null), k = log(n)
)
# automatic forward and backward model building with the BIC
# for this data set all techniques result in the same model

#############################
# Post-selection shrinkage	#
#############################
# regression coefficients
library(shrink)
# re-fit the model with x = TRUE and y = TRUE, which returns
# information (model matrix and the response)
# required by the shrink function
frog.glm3_s <- glm(
  pres.abs ~ logdistance + logNoPools +
    +avrain + I(meanmax + meanmin), family = binomial,
  data = frogs, na.action = "na.fail", x = TRUE, y = TRUE
)

# global shrinkage shrinks all parameters with the same factor
shrink_res1 <- shrink(frog.glm3_s, type = "global")
shrink_res1
# reproduce results
coef(frog.glm3_s)[-1] * shrink_res1$ShrinkageFactors

# parameterwise shrinkage (combined with backward selection)
shrink_res2 <- shrink(frog.glm3_s, type = "parameterwise")
shrink_res2

############################
# Shrinkage with the LASSO #
############################
library(glmnet)
# fit model with lasso
# prepare data set with explanatory variables
data_env <- frogs[, names(frogs) %in% c(
  "logdistance", "logNoPools",
  "NoOfSites", "avrain"
) ]
data_env$meanmaxmin <- frogs$meanmin + frogs$meanmax

# now we can fit the model
lasso_mod <- glmnet(as.matrix(data_env), frogs$pres.abs, family = "binomial")

# plot result
par(cex = 1.5)
plot(lasso_mod, label = TRUE, xvar = "lambda")

# coefficients for different models
coef(lasso_mod)

# overall, variables 3, 4 and 5 (equivalent to logNoPools, logdistance, meanmaxmin)
# can be interpreted as most important,
# as their betas are unequal to zero for highest penalty

# we determine the optimal lambda based on cross-validation
set.seed(222)
cvfit <- cv.glmnet(as.matrix(data_env), frogs$pres.abs, family = "binomial")
plot(cvfit)
# both the minimum lambda and that within 1se have 4 variables

# values for both lambdas
cvfit$lambda.min
cvfit$lambda.1se

# we extract the betas for these lambdas
coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.1se")
# regression coefficients are stronger shrunk
# for the higher lambda (related to 1se)

###########################################
# GLM Model diagnostics and interpretation #
###########################################

library(statmod)
# 1. let us check the dispersion of the model
summary(frog.glm3)
# if we divide the Residual Deviance by the
# degrees of freedom we yield a dispersion parameter
# of approximately 1 which is very good
# (remember that >2 or <0.5 would indicate over- or underdispersion)

frog.glm3$deviance / frog.glm3$df.resid
# formal calculation

# let us look at the qq plot
qr <- qres.binom(frog.glm3)
# we compute Dunn-Smyth-Residuals that are most appropriate for GLMs
qqnorm(qr)
# also relatively good fit to the qqline

# in case of over- or underdispersion you would refit
# the model with another distribution e.g. negative binomial
# or quasi likelihood estimation (e.g. for the binomial model use glm(..., family="quasibinomial")
# see the notes to the slides for how to decide between negative binomial and
# the quasibinomial model

# another way to check the residuals are so-called rootogramms
# these are particularly relevant for count data (e.g. to decide between the poisson
# and the negative binomial model)
# but are very useful for the binomial data we have here
# uncomment to install package
# install.packages("countreg", repos="http://R-Forge.R-project.org")
# library(countreg)
# rootogram(<glm.model>)

## you can also look at the residuals per component of the model
residualPlots(frog.glm3, type = "pearson")
# the conditional mean function should be constant as we move across the plot.
# This is largely fulfilled in our case
# See Fox & Weisberger (2011) p. 317-320 for more details

# 2. checking for linearity
crPlots(frog.glm3)
# no sign of non-linearity
# red dashed line indicates linear fit, green line smoothed line

## 3. check for influential observations
influenceIndexPlot(frog.glm3, vars = c("Cook", "hat"), id.n = 3)
# 182 predictor outlier, 77 has highest Cook`s distance
# but nevertheless nothing to worry

plot(frog.glm3, which = 5)
# confirms evaluation, shows that points 77,76 and a few others are outliers

influencePlot(frog.glm3, id.n = 3)
# similar plot - cooks distance is plotted as increasing bobble

# we could also look at the values numerically:
influence.measures(frog.glm3)

# we can check whether removal of the points with the highest CookD and hat values results in model change:
compareCoefs(frog.glm3, frog.glm5 <- update(frog.glm3, subset = -c(71, 75, 76, 77, 86, 145, 182)))
# model fit changes slightly - new estimates within one standard error
# but note that you would need to shrink again (post selection shrinkage)
# or run the lasso with the reduced data
# however, assuming that new data would also contain some outliers
# predictive accuracy is higher if keeping these data points

# plot to show response to individual variables
par(mfrow = c(2, 2))
termplot(frog.glm3)

# How to interpret the values on the y axis?
library(faraway)
ilogit(-2)
# 0.12
ilogit(4)
# 0.98

# matches approximately with value range
min(frog.glm3$fitted.values)
max(frog.glm3$fitted.values)

# the termplot also helps us to understand what the regression coefficients
# in the GLM summary mean: they provide the slope for the lines
# i.e. they indicate how the log odds change for one unit of increase in the
# explanatory variable.
summary(frog.glm3)
# For example, an increase of logdistance by 1 unit (e.g. from 7 to 8)
# results in a change of logodds of -0.85 (compare to the plot)

# Finally, we employ cross validation to determine the prediction error
CVbinary(frog.glm3, nfolds = 10)
# 10-fold cross validation
# gives proportion that are correctly predicted (true present, true absent)

########################################
#    Further diagnostics and analysis  #
########################################
# Further examples and diagnostics are provided in Logan (2010, chapter 17),
# Maindonald (2010, chapter 8) and Fox & Weisberger (2011, chapter 6)

#########################################
#   		 Exercise  				  #
#########################################

# Evaluate which environmental variables are most important for the occurence of the Bradypus?

# read presence absence data for the Bradypus (column pa) and environmental variables
# see http://www.worldclim.org/bioclim for more information on the individual environmental variables
pa_data <- read.csv("http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/envtrain.csv/at_download/file")

## check data
head(pa_data)
str(pa_data)

# biome should be a factor!
pa_data$biome <- factor(pa_data$biome)

# ###create plot for slides
# library(faraway)
# png("GLM.png", bg="transparent")
# par(cex=1.5, las=1)
# # thinning of absence data with annual mean temperature > 20 in order to get a better looking figure
# nr_bio <- nrow(pa_data[pa_data$bio1 > 20 & pa_data$pa==0,])
# # 591 rows
# ind_bio <- which(pa_data$bio1 > 20 & pa_data$pa==0)
# # index of absence with annual mean temp > 20
# thinning <- sample(ind_bio, nr_bio-10)
# padat2 <- pa_data[-thinning, ]
# plot(pa ~ bio1, padat2, ylim = c(0,1), xlab="Annual mean temperature", ylab="Probability of Occurence")
# lmod <- lm(pa ~ bio1, padat2)
# abline(lmod, lwd=2, col="red")
# logitmod <- glm(pa ~ bio1, family=binomial, padat2)
# x <- seq(min(padat2$bio1), max(padat2$bio1),0.1)
# lines(x,ilogit(logitmod$coefficients[1]+logitmod$coefficients[2]*x), lwd=2, col="green")
# dev.off()
# # model diagnostics for linear regression model
# par(mfrow=c(2,1))
# plot(lmod, which=2:3)
