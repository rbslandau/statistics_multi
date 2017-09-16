#################################################################################
# Script for Session 2 of the course: Applied multivariate Statistics with R	#
#					by Ralf B. Schäfer											#
# 					WS 2017/18													#
# The script describes multiple linear regression analysis in R					#
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

# Research question: 
# Which patterns and environmental variables control the diversity of marine ostracods?
# Hence, we look for the most important variables, assuming a linear relationship of the
# explanatory variables with the response variable (species richness after rarefaction)
# 
# load data
data_oc <- read.csv("http://datadryad.org/bitstream/handle/10255/dryad.39576/OstraMRegS400JB.txt?sequence=1", 
                    sep = "\t")
# the data relates to the article:
# Yasuhara M, Hunt G, van Dijken G, Arrigo KR, Cronin TM, Wollenburg JE (2012)
# Patterns and controlling factors of species diversity in the Arctic Ocean. 
# Journal of Biogeography 39(11): 2081-2088. 
# http://dx.doi.org/10.1111/j.1365-2699.2012.02758.x
# and has been taken from the DRYAD repository
# Yasuhara M, Hunt G, van Dijken G, Arrigo KR, Cronin TM, Wollenburg JE (2012) 
# Data from: Patterns and controlling factors of species diversity in the Arctic Ocean. 
# Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.9gc21

# the variables are given as:
# E100: species richness rarefied to 100 individuals, 
# SR: species richness, 
# MDS 1: multidimensional scaling first axis, 
# MDS 2: multidimensional scaling second axis, 
# DP: water depth (m), 
# RC: region code, 
# BT: bottom water temperature, 
# SA: salinity, 
# SP: seasonality of productivity, 
# IC: number of ice-free days per year, 
# P: surface productivity, 
# DCA 1: detrended correspondence analysis first axis, 
# DCA 2: detrended correspondence analysis second axis, 
# LAT: latitude, 
# LON: longitude.

#############################################################
# Checking distribution of variables and transformation     #
#############################################################

# We first remove the community information that is captured in the DCA and MDS
data_oc2 <- data_oc[ , !names(data_oc) %in% c("MDS1", "MDS2", "DCA1","DCA2")]

# now we look whether any of the variables should be transformed
par(mfrow = c(3, 3))
hist(data_oc2$DP)
hist(data_oc2$RC)
hist(data_oc2$BT)
hist(data_oc2$SA)
hist(data_oc2$SP)
hist(data_oc2$IC)
hist(data_oc2$P)
hist(data_oc2$LAT)
hist(data_oc2$LON)
# particularly DP and BT appear skewed
# only DP has a relatively wide range
dev.off()
# we transform DP using the log
hist(log10(data_oc2$DP))
# less skewed - add to data set
data_oc2$DPlog <- log10(data_oc$DP)

##############################
# Checking for collinearity  #
##############################

# we continue by checking for collinearity between exploratory variables
# we remove the response variable and non-transformed depth before checking collinearity
data_check <- data_oc2[ , !names(data_oc2) %in% c("E100", "DP", "SR")]

pairs(data_check)
# displays scatterplot between variables,

# code for a more informative plot for exploratory analyses was posted on the R graph gallery
# but is now offline, a similar function can be found in the help of ?pairs [ see 
# examples at bottom of the help]
# the function is given as:
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    
    test <- cor.test(x,y)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.05, 0.1, 1),
                  symbols = c("*", ".", " "))
    
    text(0.5, 0.5, txt, cex = cex * r)
    text(.8, .8, Signif, cex = cex, col = 2)
}
pairs(data_check, 
      lower.panel = panel.smooth, upper.panel = panel.cor)
# The lower triangle shows the raw data and a smooth function (red line) indicates the main trend
# the upper panel shows the correlation (size varies with the value)
# and stars indicate the results of a correlation test
# note that the significance tests should be interpreted with caution!

# high correlation between SP and IC, IC and P,
# SP and P, IC and BT

cor(data_check, use = "complete.obs") 
# shows correlation matrix
# we will encounter another useful technique later in the course: variable clustering

# however, we check the VIF when considering all variables
# we fit a model with all variables
mod_vif <- lm(data_oc2$E100 ~ ., data = data_check)

library(car) 
# library that contains several functions for linear model diagnosis
vif(mod_vif) 
# variance inflation factor, above 5 for several variables:
# BT, SP, IC

# we can reproduce the result for IC VIF:
mod_vifb <- lm(data_check$IC ~ ., data = data_check[ , names(data_check) != "IC"])
summary.lm(mod_vifb)#
1 / (1 - 0.9363)
# gives similar VIF

# IC correlates with many variables, including temperature
# we assume that temperature, productivity are more important 
# than IC, thus we remove IC
data_env <- data_check[ , names(data_check) != "IC"]

# now we fit the full model
# according to the paper, no interactions need to be considered
# for the sake of demonstration we include an interaction term here
# between temperature and DPlog - assuming an interaction between depth
# and temperature (irrespective whether this makes much sense)
mod_1 <- lm(data_oc2$E100 ~ . + BT:DPlog, data = data_env, na.action = 'na.fail')
summary(mod_1)
# all variables except for BT, Latitude and longitude significant (or close to)

# Harrell 2015 suggests that for each variable 10-20 cases should be available
# for linear regression to yield reliable estimates (in the sense of predictive power) 
# In our case study, we are within this range and do not apply measures to reduce
# the number of variables

######################################
# Traditional modelling approaches   #
######################################

# ideally, we would construct a few alternative models a priori
# driven by previous studies. 
# However, the authors of the original studies did not follow this
# strategy and given my limited knowledge on marine ostracods, 
# we do not implement this strategy
# Notwithstanding, this would be fairly easy from the perspective of scripting
# - you would fit and compare several models using one of the goodness of fit measures described below

#######################################################
# Computation of all possible models  #
#######################################################

library(MuMIn)
# we use the mod_1 (see above) as most complex model. 
# The dredge function computes all simpler models and evaluates them in terms of goodness of fit
models <- dredge(mod_1)
print(models)
# best fit model contains all variables except LAT and BT:DPl

##########################################################
# Average model parameters for several best-fit models   #
##########################################################

# theoretical background for model averaging can be found here: 
#
# Kenneth P. Burnham and David R. Anderson 2002: Model Selection and Multimodel Inference (2nd ed.). 
# New York: Springer-Verlag, 2002
# See also a paper of Dormann for an ecological application of this topic: 
# Dormann, C.F. et al. 2008. Prediction uncertainty of environmental change effects 
# on temperate European biodiversity. Ecology Letters 11, 235-244.
# Note also the criticism on model averaging:
# Cade B.S. (2015) Model averaging and muddled multimodel inferences. Ecology 96, 2370–2382.

# get all models that are a maximum of 2 in terms of AICc higher 
mods <- get.models(models, subset = delta < 2)
mods
# and average model parameters across these models 
avg_model <- model.avg(mods) 
summary(avg_model)
# also variables that were not significant are included in the model now

# compute r2 for averaged model (full average)
# # first we compute the predicted values using predict
fit_y <-  predict(avg_model)

# finally we compute the residual and total variance
res_var_avg <- sum((data_oc2$E100 - fit_y)^2)
tot_var_avg <- sum((data_oc2$E100 - mean(data_oc2$E100))^2)
# and the Rsquared
1 - res_var_avg/tot_var_avg
# not necessarily highest r2 (when compared to models above)

##################################################################
# stepwise backward model selection using hypothesis testing on variables
##################################################################

# we start with the full model
summary(mod_1)
# all variables except for BT, Latitude and longitude significant (or close to)

### modifying the model- removing interaction where t-test 
### does not reject H0 that beta = 0
mod_2 <- update(mod_1, ~. -BT:DPlog) 
# this removes the specified variable from the model
summary.lm(mod_2)
# LAT not significant
mod_3 <- update(mod_2, ~. -LAT)
summary.lm(mod_3)
# now all variables except SA statistically significant
# we remove SA as well
mod_4 <- update(mod_3, ~. -SA)
summary.lm(mod_4)
# now all variables except RC statistically significant
# despite close to significance, we remove RC
# as we intend to identify only the most important variables
mod_5 <- update(mod_4, ~. -RC)
summary.lm(mod_5)
# now all variables statistically significant

par(mfrow = c(2, 2))
plot(mod_5)
# model looks ok
dev.off()

# again the qreference function is helpful to check whether the residuals 
# follow a normal distribution
library(DAAG)
qreference(residuals(mod_5), nrep = 8)
# here we cannot detect any difference to the reference distribution


##############################################
# Extraction of information from the model   #
##############################################

coef(mod_5)
# gives the coefficients for the regression equation

confint(mod_5)
# gives the confidence intervals for the estimated coefficients

fitted(mod_5)
# fitted values

str(mod_5)
# overview on all available information



#################################################################
# stepwise backward model selection relying on partial F-test 	#
#################################################################

# is the explained variance significantly different between the models?
anova(mod_1, mod_2)
# no significant difference in the model fit
# this works only for nested models i.e. when the model with more variables includes 
# all variables of the model with less variables
anova(mod_2, mod_3)
# no significant difference in fit between 2 and 3
anova(mod_1, mod_3)
# or 1 and 3
anova(mod_1, mod_4)
# or 1 and 4
anova(mod_1, mod_5)
# or 1 and 5

###################################################################
# Stepwise backward selection using information theoretic approach#
###################################################################

AIC(mod_1)
# calculation of AIC
AIC(mod_2) 
# AIC higher slightly lower
AIC(mod_3) 
AIC(mod_4) 
AIC(mod_5) 
# AIC lowest for model 4

BIC(mod_1)
# calculation of BIC
BIC(mod_2) 
# stronger reduction compared to AIC
BIC(mod_3) 
BIC(mod_4) 
BIC(mod_5) 
# BIC lowest for model 5

# automatic calculation of the AICc can be done with the following functions
library(MuMIn)
AICc(mod_1)
AICc(mod_2)
AICc(mod_3)
AICc(mod_4)
AICc(mod_5)
# result similar to that for AIC

# in addition, the corrected AIC can be calculated manually according to:
# number of parameters (+1 for the estimated variance)
K <- length(mod_5$coefficients) + 1 
# calculation of the sample size
n = nrow(data_check) 

# extract number of parameters in the model
AIC(mod_5, k = 2) + 2 * K * (K + 1) / (n - K - 1)
# same as AICc above

##################################################
# automatic model building - stepwise modelling
##################################################

# we have set the full model above

nullmodel <- lm(data_oc2$E100 ~ 1, data = data_env)
# nullmodel: no variables, only mean

step(mod_1, direction = "backward", trace = 100, 
     scope = list(upper = mod_1, lower = nullmodel), k = log(n))
# automatic backward model building with the BIC - yields to same model
# you could start with any other model that has less variables
# Tools for explore different models visually are contained in the package meifly
# see http://had.co.nz/classifly/ for more details
# coincides with manual model building based on hypothesis testing

step(nullmodel, direction = "forward", trace = 100, 
     scope = list(upper = mod_1, lower = nullmodel), k = log(n))
# automatic forward model building with the BIC
# yields to a less parsimonious model, also including SA
# BIC of backward selection is lower
# main reason is the high correlation between SA and DPlog
# if DPlog is in the model, SA is not required it seems
# however, when starting with the null model, SA has the highest individual 
# explanatory power

# as discussed in the lecture stepwise selection has several problems
# we first run post-selection shrinkage that may yield to more accurate
# regression coefficients
library(shrink)
# we need to fit the model again with x = TRUE and y = TRUE, which returns 
# information (model matrix and the response)
# required by the shrink function
mod_5_s <- lm(data_oc2$E100 ~ .- LAT -RC - SA, data = data_env, 
              x = TRUE, y = TRUE)

# global shrinkage shrinks all parameters with the same factor
# should be used 
shrink_res1 <- shrink(mod_5_s, type = "global")
shrink_res1
# reproduce results
coef(mod_5_s)[-1] * shrink_res1$ShrinkageFactors 
# note that the intercept has been removed here because result would be incorrect

# parameterwise shrinkage (combined with backward selection) resulted in 
# the most accurate results in a previous simulation study
shrink_res2 <- shrink(mod_5_s, type = "parameterwise")
shrink_res2
# for categorical variables represented by dummy variables use "joint" shrinkage
# see Dunkler D., Sauerbrei W. & Heinze G. (2016) Global, Parameterwise and 
# Joint Shrinkage Factor Estimation. Journal of Statistical Software 69, 1–19. 
# free to download: https://www.jstatsoft.org/index.php/jss/article/view/v069i08

# as an alternative to post-selection shrinkage, bootstrapping and cross validation 
# of the whole stepwise selection procedure have been suggested to reduce overfitting

###########################################
# bootstrapping of stepwise modeling      #
###########################################

library(bootStepAIC)
# we need the response variable and explanatory variables in one dataframe
data_env2 <- data_env
data_env2$resp_e100 <- data_oc2$E100
set.seed(1111) # to allow for reproducible analysis
boot.stepAIC(mod_1, data_env2, k = log(n), direction = "backward")
# 'Covariates selected' provides information how often a covariate was selected
# for the final model in stepwise selection
# interestingly, for the BIC, only few variables were selected

boot.stepAIC(mod_1, data_env2, k = 2, direction = "backward")
# less restrictive GOF measure (AIC), still no clear picture regarding
# which variables should be included/removed
# Note that Austin (2008) cautioned against the use of bootstrapping
# to fix the problems of stepwise selection and see Harrell 2015 for
# related problems

# we do not look into the details of cross-validation for stepwise 
# regression, but the implementation of this in R has been featured in 
# a blog post
# https://www.r-bloggers.com/variable-selection-using-cross-validation-and-other-techniques/

#####################################
# accounting for sequential testing #
# in hypothesis testing    		   #
#####################################

# a new procedure corrects p-values for the sequential testing and provides
# coefficients and confidence intervals for predictors as they enter the model
# see lecture slides for more information
library(selectiveInference)
# fit model
# function requires that predictors are of class matrix and response variable of class vector
fs_model <- fs(as.matrix(data_env), data_oc2$E100)
plot(fs_model)
# displays change in coefficients over stepwise process

# inference for model
fsInf(fs_model)
# provides stopping point for chosen alpha (default = 0.1)
# and coefficients with confidence intervals for the variable
# as it enters the model
# in our case, the result coincides with that for ordinary stepwise forward selection

# extract coefficients of the final model intercepts
coef(fs_model, s = 6)
# note that here the coefficients are computed accounting for all 
# variables in the model. This explains the change in coefficients
# of several variables

############################
# Shrinkage with the LASSO #
############################
library(glmnet)
# fit model with lasso
lasso_mod <- glmnet(as.matrix(data_env), data_oc2$E100)
dev.off()
par(cex = 1.5)
plot(lasso_mod, label = TRUE)
# plotting coefficients against L1-norm (sum of absolute betas) on x axis
# top numbers provide information on betas different from 0
# 
plot(lasso_mod, label = TRUE, xvar = "lambda")
# plotting against lambda (the penalty term)
# 
plot(lasso_mod, label = TRUE, xvar = "dev")
# plotting against explained variance
# 
coef(lasso_mod)
# coefficients for different models

# overall, variables BT, SA, DPlog and P can be interpreted
# as most important, as their betas are unequal to zero for highest penalty
# SA would be single most important

# we determine the optimal lambda based on cross-validation
set.seed(222)
cvfit <- cv.glmnet(as.matrix(data_env), data_oc2$E100)
plot(cvfit)
# given that the lowest lambda corresponds to the minimum lambda used,
# we run the function again, including lower lambdas
set.seed(111)
cvfit2 <- cv.glmnet(as.matrix(data_env), data_oc2$E100, lambda = 10^seq(-5,0, 0.1))
plot(cvfit2)

# provides information on number of variables in model (above plot)
# CV MSE, log lambda and two dashed lines indicating the lowest CV MSE
# and the model within one standard error
cvfit2$lambda.min
cvfit2$lambda.1se

# we extract the coefficients for these lambdas
coef(cvfit2, s = "lambda.min")
coef(cvfit2, s = "lambda.1se")
coef(cvfit2)
# contains all variables except LAT

# compare to coefficients from ordinary stepwise model building
# and from parameterwise shrinkage
coefs <- data.frame(coef(mod_3), coef(cvfit2)@x)
names(coefs) <- c("Ord_stepwise", "Lasso_shrink")
coefs
# generally lasso leads to strongest shrinkage except for the most important variable SA

# a hypothesis testing approach for the lasso is introduced in:
# Taylor J. & Tibshirani R. (2016) Post-selection inference for L1-penalized 
# likelihood models. ArXiv e-prints.
# freely avaible from: https://arxiv.org/abs/1602.07358

#################################################
# Model assessment with k-fold cross validation #
#################################################
# function provided by Kabacoff 2011:214
# prediction accuracy is validated using cross validation, execute all code below to set the function
# this makes only sense for models where CV has not been employed in model selection
# (i.e.not meaningful for lasso model above)
shrinkage <- function(fit, k = 10){
  library(bootstrap)
  theta.fit <- function(x, y){
    lsfit(x, y)
  }
  theta.predict <- function(fit, x){
    cbind(1, x) %*% fit$coef
  }
  x <- fit$model[ , 2:ncol(fit$model)] 
  y <- fit$model[ , 1]
  results <- crossval(x, y, theta.fit, theta.predict, ngroup = k) 
  r2 <- cor(y, fit$fitted.values)^2 
  r2cv <- cor(y, results$cv.fit)^2 
  cat("Original R-square =", r2, "\n")
  cat(k, "Fold Cross-Validated R-square =", r2cv, "\n") 
  cat("Change =", r2 - r2cv, "\n") 
}
#### end of function

#calculate validation for model 4 and 5
shrinkage(mod_4) 
# r-square is about ~5% reduced
shrinkage(mod_5) 
# r-square is about ~6% reduced

# CV does often underestimate the true test error, because data from the 
# same data set is used in validation
# nevertheless, it is useful for identifying the correct level of flexibility
# as implemented above

#####################################
# Relative importance of a variable #
#####################################

library(relaimpo)
pred_imp <- calc.relimp(mod_5, 
                        type = c("lmg", "first", "last", "betasq"), 
                        rela = TRUE)
# check out the help for this function for information on these relative importance measures
?calc.relimp
plot(pred_imp)
# BT, P and DPlog as most important variables

# pmvd is currently only available from the website 
# http://prof.beuth-hochschule.de/groemping/software/relaimpo/
# due to copyright restrictions

# lmg is a special case of hierarchical partitioning for linear models and R2
# hierarchical partitioning is implemented in another package that also allows for 
# the choice of different gof measures whereas relaimpo relies on R2
# note that Grömping (2015) cautions against using more than 9 variables in the hier.part function
# Grömping, U. (2015). Variable importance in regression models. WIREs Comput Stat 7, 137-152

# let us regard model 3 here
library(hier.part)
dev.off()
hier.part(data_oc2$E100, data_env[ , !(names(data_env) %in% "LAT" )], 
          gof = "Rsqu")

# one reason for the different results is the correlation between 
# DPlog and SA, if we knew that only one of these is important and
# the other one a spurious effect, we could remove one of them before analysis

################################################################################
# Exercise 4 : For an effective environmental protection scheme you need to know 
#  causes of pollution. 
# In this example, the aim is to find the variables that are most important to 
# predict the SO2 air concentrations.
# Load the data set US Air Pollution and model the SO2 concentration
# as response and the other variables as predictors. 
# Compare the results for manual model building using hypothesis testing,
# automatic backward model selection with post-selection shrinkage and the 
# Lasso shrinkage 
################################################################################

data("USairpollution", package = "HSAUR2") 
# see package information for details on data set
head(USairpollution)
