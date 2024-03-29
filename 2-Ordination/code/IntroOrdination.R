### ----------------------------------------------------------------------------
### R Script for the workshop: "Data analysis in freshwater ecology using R" @SEFS9
### Part 2: A (brief) introduction to ordination and the vegan package
### Written by Eduard Szöcs, 03.07.2015

setwd("/home/edisz/Documents/Uni/Projects/PHD/CONFERENCES/SEFS9_Geneva/2-Ordination/data/")
### --------------
### Just a small function to reset the graphics window
### Written by Gavin Simpson
### Source: http://stackoverflow.com/questions/5789982/reset-par-to-the-default-values-at-startup
### Usage: par(resetPar()) 
resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}
### ----------------------------------------------------------------------------
### Datasets

### ---------------
### Doubs river data
# Load Doubs data sets
Dabu <- read.table('doubsAbu.csv', sep = ',', header = TRUE)
Denv <- read.table('doubsEnv.csv', sep = ',', header = TRUE)
Dspa <- read.table('doubsSpa.csv', sep = ',', header = TRUE)
# Dimensions of tables (30 sites,  27 taxa)
dim(Dabu)
head(Dabu[ , 1:18])
# Dimension and first rows of Environmental data
dim(Denv)
head(Denv)

### ---------------
### Melbourne data
# setwd('3-Ordination/data/')
abu <- read.table('melbourneAbu.csv', sep = ';', header = TRUE)
env <- read.table('melbourneEnv.csv', sep = ';', header = TRUE)
# dimensions of data.frame
dim(env)
dim(abu)
### first rows of environmental data (only first 10 variables)
head(env[ , 1:10])


### ---------------
### Artificial data
# Load dummy data
dummy <- read.table('dummydata.csv', header = TRUE, sep = ';')
# plot dummy data
matplot(dummy[ , -1], type = 'l', xlab = 'Site', ylab = 'Abundance', 
        lty = 'solid', lwd = 2, col = 1:9)
legend('topright', legend = colnames(dummy)[-1],
        col = 1:9, lty = 'solid', lwd = 2)



### ----------------------------------------------------------------------------
### Unconstrained Ordination

### ----------
# Slide 12 - plot environmental variables along doubs river
par(mfrow = c(2, 2))
plot(Dspa, asp = 1, pch = 21, cex = 5*Denv$alt / max(Denv$alt), bg = 'darkred',
     main = 'Altitude')
lines(Dspa, col = 'steelblue', lwd = 2)
plot(Dspa, asp = 1, pch = 21, cex = 5*Denv$deb / max(Denv$deb), bg = 'darkred',
     main = 'Discharge')
lines(Dspa, col = 'steelblue', lwd = 2)
plot(Dspa, asp = 1, pch = 21, cex = 5*Denv$nit / max(Denv$nit), bg = 'darkred',
     main = 'Nitrate')
lines(Dspa, col = 'steelblue', lwd = 2)
plot(Dspa, asp = 1, pch = 21, cex = 5*Denv$pho / max(Denv$pho), bg = 'darkred',
     main = 'Phosphate')
lines(Dspa, col = 'steelblue', lwd = 2)
par(resetPar()) 
### 3-D Plot
require(rgl)
plot3d(Denv[ , c('das', 'alt', 'pho')], pch = 16, size = 5)

### ----------
### PCA
require(vegan)
PCA <- rda(Denv, scale = TRUE)
plot(PCA, scaling = 3)

## Create biplot with arrows. 
## Scaling 3 is symmetric scaling, all interpretations are approximates.
biplot(PCA, cex = 5, scaling = 1)



### -----------
## Create biplot with arrows and color ellipses added
# biplot
bp <- biplot(PCA, cex = 5, scaling = 3)
# colors to use
cols <- rainbow(5, start = 2/6)
# group sites
g <-  c(rep('g1', 10), rep('g2', 6), 
                        rep('g3', 6), rep('g4', 4),
                        rep('g5', 4))
# for every group add one ellipse
for (i in 1:5) {
  ordiellipse(bp, groups = g, 
              show.groups = unique(g)[i], 
              col = cols[i], 
              lwd = 4)
}
# See ?ordiellipse for other possible and usefull decorations

### ----------
### Numerical summary output
summary(PCA)
summary(PCA, display = NULL, scaling = 3)

### ----------
### Exercise 1 - Solution

### -----------------------------------------------------------------------------
### Excursus | principal component regression (PCR)

### PCA from Exercise 1
take <- env[ , !names(env) %in% c('ID', 'logCond', 'logmaxTU')]
PCA <- rda(take, scale = TRUE)
# calculate shannon diversity index
div <- diversity(abu[ , -1], index = 'shannon')
pc <- scores(PCA, choices = c(1, 2), scaling = 1, display = 'sites')
model_data <- data.frame(div, pc, logCond = env$logCond, logmaxTU = env$logmaxTU)
model <- lm(div ~ PC1 + PC2 + logCond + logmaxTU, data = model_data)
summary(model)



### -----------------------------------------------------------------------------
### Ordination of abundance data

### -----------
### Plot abundance of four species along the Doubs river

par(mfrow = c(2, 2))
plot(Dspa, asp = 1, pch = 21, cex = 3*Dabu$TRU / max(Dabu$TRU), bg = 'darkred',
     main = 'Brown trout')
lines(Dspa, col = 'steelblue', lwd = 2)
plot(Dspa, asp = 1, pch = 21, cex = 3*Dabu$BCO / max(Dabu$BCO), bg = 'darkred',
     main = 'Bream')
lines(Dspa, col = 'steelblue', lwd = 2)
plot(Dspa, asp = 1, pch = 21, cex = 3*Dabu$CHA / max(Dabu$CHA), bg = 'darkred',
     main = 'Bullhead')
lines(Dspa, col = 'steelblue', lwd = 2)
arrows(25, 133, 94, 45, code = 3, lwd = 7, col = 'darkorange')
text(30, 70, '?', cex = 2, col = 'darkorange')
plot(Dspa, asp = 1, pch = 21, cex = 3*Dabu$OMB / max(Dabu$OMB), bg = 'darkred',
     main = 'Grayling')
lines(Dspa, col = 'steelblue', lwd = 2)
arrows(25, 133, 94, 45, code = 3, lwd = 7, col = 'darkorange')
text(30, 70, '?', cex = 2, col = 'darkorange')
par(resetPar()) 

### -----------
### (Dis-) Similarity measures

# Example with dummy data (3 sites, 3 species)
mat <- matrix(c(0, 4, 8, 0, 1, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
colnames(mat) <- c('Spe1', 'Spe2', 'Spe3')
rownames(mat) <- c('sit1', 'sit2', 'sit3')
mat
# Calculate euclidean distance between sites
vegdist(mat, method = 'euclidean')

# Calculate Bray-Curtis distance between sites
vegdist(mat, method = 'bray')


### -----------------------------------------------------------------------------
### Principal Coordiante Analysis (PCoA)

# Distance matrix
Dabu_dist <- vegdist(Dabu, method = 'bray')

# PCoA
PCOA <- cmdscale(Dabu_dist, eig = TRUE)

# Create plot
plot(PCOA$points, type = 'n', 
     xlab = 'PCOA1', ylab = 'PCOA2')
text(PCOA$points, 
     labels = rownames(Dabu), cex = 0.9)
abline(h = 0 , lty = 'dotted')
abline(v = 0 , lty = 'dotted')
# Add species as weighted averages
wa <- wascores(PCOA$points, Dabu)
text(wa, labels = colnames(Dabu), 
     col = 'red', cex = 0.7)
# explained variance
(PCOA$eig / sum(PCOA$eig))[1:2] * 100


### -----------------------------------------------------------------------------
### Nonmetric Multidimensional Scaling (NMDS)

# Distance matrix
Dabu_0 <- Dabu[!rowSums(Dabu) == 0, ]
Dabu_dist <- vegdist(Dabu_0, method = 'bray')

# NMDS
NMDS <- metaMDS(Dabu_dist, k = 2)

# Plot
plot(NMDS, type = 't')

# Add species as weighted averages
wa <- wascores(NMDS$points, Dabu_0)
text(wa, labels = colnames(Dabu), 
     col = 'red', cex = 0.7)

# Stress value
NMDS$stress

### ----------
### Exercise 2 - Solution



### -----------------------------------------------------------------------------
### Indirect Gradient Analysis

# PCoA of fish community data
Dabu_dist <- vegdist(Dabu, method = 'bray')
PCOA <- cmdscale(Dabu_dist, eig = TRUE)

# plot PCoA
plot(PCOA$points, type = 'n', 
     xlab = 'PCOA1', ylab = 'PCOA2', main = 'PCoA of fish community data')
text(PCOA$points, 
     labels = rownames(Dabu), cex = 0.9)
wa <- wascores(PCOA$points, Dabu)
text(wa, labels = colnames(Dabu), 
     col = 'red', cex = 0.7)
abline(h = 0 , lty = 'dotted')
abline(v = 0 , lty = 'dotted')

# Plot PCoA with Superimposed Altitude (displayed as point size)
plot(PCOA$points,
     xlab = 'PCOA1', ylab = 'PCOA2', 
     cex = 5*Denv$alt / max(Denv$alt),
     main = 'PCoA with Altitude superimposed', 
     bg = 'grey75', pch = 21)
wa <- wascores(PCOA$points, Dabu)
text(wa, labels = colnames(Dabu), 
     col = 'red', cex = 0.7)
abline(h = 0 , lty = 'dotted')
abline(v = 0 , lty = 'dotted')

# Plot PCoA with Superimposed Altitude (displayed as point size)
plot(PCOA$points,
     xlab = 'PCOA1', ylab = 'PCOA2', 
     cex = 5*Denv$alt / max(Denv$alt),
     main = 'PCoA with Altitude superimposed', 
     bg = 'grey75', pch = 21)
abline(h = 0 , lty = 'dotted')
abline(v = 0 , lty = 'dotted')

# Fit Altitude to site-scores
ef <- envfit(PCOA, Denv)
plot(ef)
# summary with R^2 and permutation tests
ef

# Fit Generalized Additive Model to ordination
ordisurf(PCOA, Denv$alt, add = TRUE)



### ----------------------------------------------------------------------------
### Direct Gradient Analysis
### ----------
### RDA
RDA <- rda(Dabu ~ ., data = Denv, scale = TRUE)
plot(RDA, scaling = 3)

# summary
summary(RDA)
# * Partitioning of variance (constrained / unconstrained)
# * Explained variance
# * RDA axes: Variance explained by the model 
# * PCA axes: variance in the residuals (=PCA on residuals)
# * fairly high amount of residual variance

### ----------
### Transformations
mat
decostand(mat, 'hellinger')

### ----------------------------------------------------------------------------
### transformation-based RDA (tb-RDA)
# Hellinger transformation
Dabu_h <- decostand(Dabu, 'hellinger')
# RDA on Hellinger transformed abundances
tbRDA <- rda(Dabu_h ~ ., data = Denv)
# no additional scaling!

# Plot
plot(tbRDA, type = 't')


### ----------------------------------------------------------------------------
### distance-based RDA (db-RDA)
dbRDA <- capscale(Dabu ~ ., data = Denv, 
                  distance = 'bray')
plot(dbRDA, type = 't')
summary(dbRDA)

### ----------
### Exercise 3 - Solution



### ----------------------------------------------------------------------------
### Model building / diagnostics / testing 

# visual inspetion via PCA
biplot(rda(Denv, scale = TRUE), 
       display = 'species', scaling = 2)
# RDA on Hellinger transformed abundances
tbRDA <- rda(decostand(Dabu, 'hellinger') ~ ., data = Denv)
# Variance inflation factors
vif.cca(tbRDA)
# GOF for each species on 
goodness(tbRDA)[ , 1:3]
# total var explained by constraints
goodness(tbRDA, summarize = TRUE)[1:6]
# variance explained by constrained (CCA) and unconstrained (CA) axes
inertcomp(tbRDA, proportional = TRUE)
densityplot(permustats(permutest(tbRDA)))
# RDA on Hellinger transformed abundances
tbRDA <- rda(decostand(Dabu, 'hellinger') ~ alt + oxy + pH + nit + pho,
             data = Denv)
vif.cca(tbRDA)
# Tests if the overall model is significant
anova(tbRDA)
# Tests RDA axes
anova(tbRDA, by = 'axis')
# Tests RDA axes
anova(tbRDA, by = 'terms')
# Tests RDA axes
anova(tbRDA, by = 'margin')
