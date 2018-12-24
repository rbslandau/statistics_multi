###################################################################################
# Script for Session 5 of the course: Applied multivariate Statistics with R  	  #
# 						by Ralf B. Schäfer,	WS 2018/19							  #
# 						Redundancy Discrimnant Analysis					  		  #
###################################################################################

# we first set a working directory i.e. a directory where we store all files
setwd("~/Gitprojects/Teaching/Statistics_multi/Code")
# you have to set a path to a working directory on your local machine here
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

# parts of the code of this script
# have been adapted from Borcard et al. 2011
# and Zuur et al. 2007

library(vegan)

# load lichen data
data(varespec)
head(varespec)

# Here we are interested in species with a strong relationship with the environmental gradient,
# consequently we remove species that occur at 4 or less of the 24 sites.
# Compare the discussion on the removal of rare species in the literature:
# Cao, Y.; Larsen, D. P.; Thorne, R. S., Rare species in multivariate analysis for bioassessment:
# some considerations. J. N. Am. Benthol. Soc. 2001, 20, (1), 144-153.
# Marchant, R., Do rare species have any place in multivariate analysis for bioassessment?
# J. N. Am. Benthol. Soc. 2002, 21, (2), 311-313.
# Rare species may especially cause problems in Correspondence Analysis and Canonical Correspondence Analysis
# In this analysis, we remove them primarily to simplify the data set

# presence-absence transformation to calculate species number per site
varespec_pa <- decostand(varespec, "pa")
# calculate sum per species
vare_sum <- apply(varespec_pa, 2, sum)
sort(vare_sum)

# removal of rare species

# remove species that occur at less than 5 sites
varespec_fin <- varespec[, !vare_sum < 5]

sort(apply(varespec_fin, 2, max))
# strong differences in order of magnitude of species abundances
sort(apply(varespec_fin, 2, sd))
# species with highest abundances also dominate the variation
# without standardisation (which includes mean centering and standardisation to unit variance)
# the analysis would be dominated by the species with the highest variation
# We can standardise the data directly within the rda command (see below)

# now we look at the environmental variables
data(varechem)
head(varechem)
# environmental data from soils

# collinearity may hamper interpretation of the RDA with respect to the relevance of individual variables
# we therefore check for collinearity
# note that we could also conduct a PCA and work with orthogonal environmental variables,
# though this does not necessarily simplify interpretation

### function for an overview of environmental variables
# see script on multiple linear regression for more information
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 0.8 / strwidth(txt)

  test <- cor.test(x, y)
  # borrowed from printCoefmat
  Signif <- symnum(
    test$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.05, 0.1, 1),
    symbols = c("*", ".", " ")
  )

  text(0.5, 0.5, txt, cex = cex * r)
  text(.8, .8, Signif, cex = cex, col = 2)
}

# since there are too many variables for plotting in one go, we use subsets
pairs(varechem[, 1:7], lower.panel = panel.smooth, upper.panel = panel.cor)
# high correlations (> 0.75) between Mg and Ca as well as for K and S
# note that for a serious analysis, past research on lichens should be considered
# here we omit Ca and K based on ecological intuition (which may be elusive)

pairs(varechem[, 8:14], lower.panel = panel.smooth, upper.panel = panel.cor)
# no other problematic bivariate relationships

# check again after omission of Ca and K
varechem_red <- varechem[, -c(3:4)]

pairs(varechem_red, lower.panel = panel.smooth, upper.panel = panel.cor)
# highest correlation between Al and Fe. Given the potential toxicity of Al, we omit Fe now

varechem_fin <- varechem_red[, -6]

# we can also check the VIF
library(faraway)
sort(vif(varechem_fin))
# according to the VIF, more should be done about multicollinearity,
# but since we are no lichen specialists, we leave as is
# Borcard et al. 2011: 175 argue that in the case of RDA, VIFs > 10 should be avoided
# please refer to the part on multiple linear regression on methods to deal with multicollinearity

# checking gradient length with detrended correspondence analysis
decorana(varespec_fin)
# gives the average standard deviation
# a completely unimodal gradient has approximately 4 SD
# linear analyses are presumably ok for SD up to 2-3

# check whether Hellinger transformation decreases gradient length
# details on this transformation can be found in:
# Legendre, P.; Gallagher, E. D., Ecologically meaningful transformations for ordination of species data.
# Oecologia 2001, 129, (2), 271-280.
varespec_hel <- decostand(varespec_fin, "hellinger")

# check axis length again
decorana(varespec_hel)
# transformation decreases axis length
# nevertheless, we use original data because of easier interpretation

####################
#   conduct RDA    #
####################

# since we did not standardise before, we should scale the variables
# this implies that changes in species abundances obtain equal consideration
# if we had used the hellinger transformation, we would not need to scale
var_rda <- rda(varespec_fin ~., data = varechem_fin, scale = TRUE)
summary(var_rda)
# we obtain 11 RDA (constrained) and 12 PCA (unconstrained) axes.
# Further PCA axes are not displayed as their eigenvalues are too small
# and would exceed the number of observations

# 53% explained by environmental variables
# first two axes explain 22% of total variation
# and 22/53 (~42%) of constrained variation

# species scores represent the position of species on the axis in bi- or triplots
# site scores give the coordinates of sites in the space of the species (response)
# site constraints give the coordinates of sites in the space of the environmental variables (explanatory)
# biplot scores give the coordinates for the environmental variable arrows

# if including rare species via:
var_rda_test <- rda(varespec ~., data = varechem_fin, scale = TRUE)
summary(var_rda_test)
# the total constrained variance is higher, but explanatory
# power on first axes is lower, which suggests that
# species response data different from the main responses
# has been added

# Extraction of canonical coefficients from rda object
coef(var_rda)

# Global test of the RDA result
set.seed(111)
anova.cca(var_rda, step = 1000)
# full model not statistically significant

# Tests of all canonical axes
set.seed(111)
anova.cca(var_rda, by = "axis", step = 1000)
# only first two RDA axes statistically significant
# for more information see: Legendre, P.; Oksanen, J.; ter Braak, C. J. F.,
# Testing the significance of canonical axes in redundancy analysis. Methods in Ecology and Evolution 2011, 2, (3), 269-277.

###############################
# Triplots of the rda results #
###############################
# again there are different scalings available
# we compare two commonly used scalings

# Scaling 1: distance triplot
par(mfrow = c(1, 2))
plot(var_rda, scaling = 1, main = "Triplot RDA scaling 1 - wa scores")
var.sc <- scores(var_rda, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, var.sc[, 1], var.sc[, 2], length = 0, lty = 1, col = "red")
# lines for species can be omitted

# How to interpret this plot?
# Zuur et al 2007, page 211 give the following rules (content modified by RBS):
# 1. Sites can be projected perpendicularly on the species or explanatory lines
#    and approximate the fitted value of the site along the variable
# 2. Distances between sites represent a two-dimensional approximation
#    of their fitted Euclidean distances.
# 3. Angles between lines of species and explanatory variables reflect their correlation.
# 4. Angles between lines for species or explanatory variables do not represent correlations.
# an additional technical explanation is given in Legendre & Legendre 2012: 639-641

# Scaling 2 (default): correlation triplot

plot(var_rda, main = "Triplot RDA scaling 2 - wa scores")
var2.sc <- scores(var_rda, choices = 1:2, display = "sp")
arrows(0, 0, var2.sc[, 1], var2.sc[, 2], length = 0, lty = 1, col = "red")
# lines for species can be omitted

# 1. Sites can be projected perpendicularly on lines (species or explanatory)
#    indicating the values of the species or explanatory variables at those sites.
# 2. Distances between sites cannot be compared directly with each other.
# 3. All angles between lines reflect their correlation.
#    But these are not as accurate as the ones obtained by a PCA of the explanatory variables.
# 4. Length of lines is not important

##########################
# 		Extra Section    #
#     WA vs LC scores    #
##########################

# In the previous plots the positions of sites are based
# on the weighted averages of species positions
# there is an open discussion whether the position of the sites should
# better be expressed as linear combinations of the explanatory variables:
# see
vignette("decision-vegan", package = "vegan")
# and compare to Borcard et al. 163-164

# This would be the code for the plots in the lc scale

# Scaling 1
quartz()
# use x11() on windows to open new graphics device
plot(
  var_rda, scaling = 1, display = c("sp", "lc", "cn"),
  main = "Triplot RDA - scaling 1 - lc scores"
)
arrows(0, 0, var.sc[, 1], var.sc[, 2], length = 0, lty = 1, col = "red")
# the rules from above still apply for interpretation

# Scaling 2
quartz()
# use x11() on windows to open new graphics device
plot(
  var_rda, display = c("sp", "lc", "cn"),
  main = "Triplot RDA - scaling 2 - lc scores"
)
arrows(0, 0, var2.sc[, 1], var2.sc[, 2], length = 0, lty = 1, col = "red")
# the rules from above still apply for interpretation

##########################
#  End	Extra Section    #
# 						 #
##########################


#########################
# RDA and model building #
#########################

# the model building frameworks of the linear and generalized linear models in R
# can also be used in the case of RDA
# This makes sense if the purpose is to identify the most important variables
# and the most parsimonious RDA model

# Forward selection using the adj. R2 as goodness of fit measure
step_both_R2 <- ordiR2step(rda(varespec_fin ~ Al, data = varechem_fin), scope = formula(var_rda), direction = "both", pstep = 1000)
# this yields a best-fit model consisting of Al + P + N + Baresoil
# Note that the stepwise function has several stopping criteria:
# 1. adj r2 decreases
# 2. adj r2 would exceed that of the model defined in the scope
# this can be checked with RsquareAdj(var_rda)
# 3. exceeding selected p-value

#####################
# Model calculation #
# Parsimonious RDA  #
#####################

var_rda_pars <- rda(varespec_fin ~ Al + P + N + Baresoil, data = varechem_fin, scale = TRUE)
var_rda_pars
# 4 variables still explain 23% of total variation
# to put this into context, consider an average of ~5% explained variance in ecological studies
# see: Moller, A. P.; Jennions, M. D.,
# How much variance can be explained by ecologists and evolutionary biologists? Oecologia 2002, 132, (4), 492-500

set.seed(111)
anova.cca(var_rda_pars, step = 1000)
# model is statistically significant
set.seed(111)
anova.cca(var_rda_pars, step = 1000, by = "axis")
vif.cca(var_rda_pars)

# Scaling 1: distance triplot
quartz()
# use x11() on windows to open new graphics device
par(cex = 1.5)
plot(var_rda_pars, scaling = 1, main = "Triplot RDA scaling 1 - wa scores")
var_pars_sc <- scores(var_rda_pars, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, var_pars_sc[, 1], var_pars_sc[, 2], length = 0, lty = 1, col = "red")
summary(var_rda_pars)

###################
# Partial RDA     #
###################
# Evaluating the effect of a variable
# after accounting for the effect of
# another variable or set of variables

# Here we account for the effect of
# soil variables (via Condition(<variables>)) to evaluate the effects of
# habitat (baresoils)
var_rda_partial <- rda(varespec_fin ~ Baresoil + Condition(Al + P + N), data = varechem_fin, scale = TRUE)
summary(var_rda_partial)
# gives explained variation after accounting for soil chemistry variables
#
plot(var_rda_partial, scaling = 1, main = "Triplot RDA scaling 1 - wa scores")
var.sc <- scores(var_rda_partial, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, var.sc[, 1], var.sc[, 2], length = 0, lty = 1, col = "red")

# allows to interpret effect of Baresoil after accounting for other effects
# see Borcard et al 2011: 180ff on how to attribute variation to different
# variables in RDA

##################
#   Exercise     #
##################
# Conduct a complete RDA for the RIKZ data set (taken from Zuur et al. 2007).
# How many RDA axes are needed?
# Which variables are most important to explain the community pattern?

RIKZ <- read.table("https://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/teaching/RIKZ_data/at_download/file", header = TRUE)
# Species data
Species <- RIKZ[, 2:5]
ExplVar <- RIKZ[, 9:15]
# Background on the data - taken from Zuur et al. 2007:
# Sea level rise may affect the ecology of the Dutch coastal system
# The Dutch governmental institute RIKZ therefore started a research project
# on the relationship between some abiotic aspects
# (e.g., sediment composition, slope of the beach) and the benthic fauna.
# NAP is the height of the sampling station relative to the mean tidal level.
# Exposure is an index that is composed of the following elements:
# wave action, length of the surf zone, slope, grain size and the depth of the anaerobic layer.
# Humus constitutes the amount of organic material.
