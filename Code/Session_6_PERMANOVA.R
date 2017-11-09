###################################################################################
# Script for Session 6 of the course: Applied multivariate Statistics with R  	  #
#		  created by Eduard Szöcs,	 modified by RBS, 	WS 2017/18				  #
# 						        PERMANOVA					  		  	  		  #
###################################################################################

### ------------ Load data and package -----------------------------------------
### Load data and package
require(vegan)
werra_sp <- read.table('https://raw.githubusercontent.com/rbslandau/statistics_multi/master/Data/River_sp.csv', sep = ';', 
                       header = TRUE, row.names = 1)
werra_env <- read.table('https://raw.githubusercontent.com/rbslandau/statistics_multi/master/Data/River_env.csv', sep = ';', header = TRUE)


### ------------ Distance matrix ----------------------------------------------
### Compute distance matrix using Bray-Curtis on double-root transformed abundances
range(werra_sp)
range(werra_sp^0.5)
range(werra_sp^0.25)
# We use a double root transformation to reduce the range of the data
dist_werra <- vegdist(werra_sp^0.25, method = 'bray')
# otherwise highly abundant taxa would dominate the distance measures

### ------------ NMDS ----------------------------------------------------------
### Run NMDS
nmds <- metaMDS(dist_werra)
### Plot NMDS
op <- ordiplot(nmds, type = 'n')

# points
cols = c('red', 'green')
points(nmds, cex = 2, pch = 16, col = cols[werra_env$Position])

# decoration
ordispider(nmds, groups = werra_env$Position, label = TRUE)
ordihull(nmds, groups = werra_env$Position, lty = 'dotted')
legend("bottomleft", pch = 16, col = cols, legend = levels(werra_env$Position))
# upstream and downstream communities are well separated (indicates an effect!)
# spread may be different (lower for upstream section)

### ------------ PERMANOVA -----------------------------------------------------
pmv <- adonis(werra_sp^0.25 ~ Position, data = werra_env, 
              permutations = 999, 
              method = 'bray')
pmv
# communities differ statsitically significant between upstream/downstream
# the position explains 30.7% of variance

# plot permuted F-values
densityplot(permustats(pmv))
plot(density(pmv))

## But are the assumptions met?
# check assumptions visually
nmds <- metaMDS(dist_werra)
op <- ordiplot(nmds, type = 'n')
# points
cols = c('red', 'green')
points(nmds, cex = 2, pch = 16, col = cols[werra_env$Position])
# decoration
ordispider(nmds, groups = werra_env$Position, label = TRUE)
ordihull(nmds, groups = werra_env$Position, lty = 'dotted')
legend("bottomleft", pch = 16, col = cols, legend = levels(werra_env$position))
# graphically we would say that upstream has a slightly lower spread then downstream


### ------------ Distance based dispersion test -------------------------------
bd <- betadisper(dist_werra, werra_env$Position)
bd
# also an eigenvalue based method

boxplot(bd)
# boxplot of Average distance to median shows also that there might be lower dispersion 
# for upstream sites

# F-Test
anova(bd)

# permutaion test
permutest(bd)
# We cannot find a statistically significantly different dispersion 
# -> assumption of homogeneity is met



### ------------ SIMPER --------------------------------------------------------
sim <- simper(werra_sp^0.25, group=werra_env$Position)
summary(sim)
# contr :   contribution to dissimilarity between upstream and downstream
# sd    :   standard deviation of contribution (is the species response consitent?)
# ratio :   ratio between contr and sd (high ratio = high, consisten contribution)
# av.   : average abundance per groups
# cumsum:   cumulative contribution (rule of thumb : species till 70% are investigated)

# many species contribute to the difference between upstream and downstream sections
# Lype sp. decreases (from 6.29 to 2.49) and has the highest contribution
# Lasiocephala.basalis increases (from 3.24 to 7.09)
# Gammarus sp. descreases to zero (from 2.78 to 0) and also shows this consitently (high contribution to sd ratio)
