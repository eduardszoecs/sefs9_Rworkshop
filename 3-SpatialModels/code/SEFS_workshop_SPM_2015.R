#########################################################
#####Spatial autocorrelation in ecological modelling#####
#########################################################

## Set the working directory

setwd("/Users/avitbhowmik/PhD/Conferences/SEF2015/Workshop/Data")

## Import data

snouter.df <- read.csv("snouter.csv", header=T, sep=",")

# Alternatively, directly from GitHub
library(RCurl)
snouter <- getURL("https://raw.githubusercontent.com/EDiLD/sefs9_Rworkshop/master/3-SpatialModels/data/snouter.csv")
snouter.df <- read.csv(text=snouter)

#snouter.df <- read.csv("https://raw.githubusercontent.com/EDiLD/sefs9_Rworkshop/master/3-SpatialModels/data/snouter.csv",
	#header=T, sep=",")

## Look through the data.
str(snouter.df)
head(snouter.df)
nrow(snouter.df)
summary(snouter.df)

## Important: is there spatial information?

## Map snouter abundance
# Let's convert the snouter data into a spatial data based on spatial information

library(sp)
snouter.sp <- snouter.df
coordinates(snouter.sp) <- ~x.coord+y.coord
class(snouter.sp)
head(snouter.sp@coords)
head(snouter.sp@data)

spplot(snouter.sp["snouter.ab"], scales=list(draw=T), colorkey=T, main="Snouter abundance")

# With ggplot

library(ggplot2)
ggplot(data=snouter.df, aes(x=x.coord, y=y.coord, color=snouter.ab))+
geom_point() + coord_equal()

## Try to map snouter presence-absence

snouter.sp@data$snouter.pa.fac <- factor(snouter.sp@data$snouter.pres.abs)
spplot(snouter.sp["snouter.pa.fac"], scales=list(draw=T), main="Snouter presence-absence")

## Visualize and quantify spatial autocorrelation
# Let's extract the residuals from a non-spatial model

nonsp.snout.ols <- lm(snouter.ab ~ prec + dist.jung, data=snouter.df)
nonsp.snout.ols.res <- residuals(nonsp.snout.ols)
nonsp.snout.ols.stres <- rstandard(nonsp.snout.ols)

# Then map the residuals
snouter.sp@data$nonsp.ab.stres <- nonsp.snout.ols.stres
spplot(snouter.sp["nonsp.ab.stres"], scales=list(draw=T), colorkey=T, main="Non-spatial model residual for abundance")

# Let's make a correlogram
library(ncf)
correlog.snout.ab <- correlog(snouter.df$x.coord, snouter.df$y.coord, nonsp.snout.ols.res , na.rm=T, increment=1, resamp=0)

# and plot for the first 20 distance classes

plot(correlog.snout.ab$correlation[1:20], type="b", xlab="distance", ylab="Moran's I",  col="blue")
abline(h=0, lty=2, col="red")

## Now let's make a variogram

library(gstat)
var.snouter.ab <- variogram(nonsp.ab.stres~1, snouter.sp)
plot(var.snouter.ab)

## Less elegant way

library(nlme)
gls.nsp <- gls(snouter.ab ~ prec + dist.jung, data=snouter.df)
var.snouter.ab.lme <- Variogram(gls.nsp, form=~x.coord+y.coord, resType = "pearson")
plot(var.snouter.ab.lme, xlim=c(0,20))

## Now let's calculate global Moran's I

library(spdep)
snouter.ab.nb <- dnearneigh(snouter.sp@coords, 0, 20)
snouter.ab.listw <- nb2listw(snouter.ab.nb)
Glob.Moran.I.snout.ab <- moran.test(nonsp.snout.ols.stres , listw=snouter.ab.listw)
Glob.Moran.I.snout.ab$estimate
Glob.Moran.I.snout.ab$statistic

?corClasses
gls.snouter.ab.exp <- gls(snouter.ab ~ prec + dist.jung, data=snouter.df, correlation=corExp(form=~x.coord+y.coord))
summary(gls.snouter.ab.exp)

gls.snouter.ab.spher <- gls(snouter.ab ~ prec + dist.jung, data=snouter.df, correlation=corSpher(form=~x.coord+y.coord))
summary(gls.snouter.ab.spher)

AIC(gls.snouter.ab.exp, gls.snouter.ab.spher)
## Try with Gaussian structure and let me know the result

summary(gls.nsp)
summary(nonsp.snout.ols)

## Compute explained variance

library(MuMIn)
r.squaredLR(gls.nsp)
r.squaredLR(gls.snouter.ab.exp)

######################################################################################
### GLMM ###

nonsp.snout.glm <- glm(snouter.pres.abs ~ prec + dist.jung, data=snouter.df, family=binomial)
nonsp.snout.glm.res <- residuals(nonsp.snout.glm)
nonsp.snout.glm.stres <- rstandard(nonsp.snout.glm)

# Then map the residuals
snouter.sp@data$nonsp.pa.stres <- nonsp.snout.glm.stres
spplot(snouter.sp["nonsp.pa.stres"], scales=list(draw=T), colorkey=T, main="Non-spatial model residual for abundance")

# Correlogram
correlog.snout.pa <- correlog(snouter.df$x.coord, snouter.df$y.coord, nonsp.snout.glm.stres , na.rm=T, increment=1, resamp=0)

# and plot for the first 20 distance classes

plot(correlog.snout.pa$correlation[1:18], type="b", xlab="distance", ylab="Moran's I",  col="blue")
abline(h=0, lty=2, col="red")

## Variogram

var.snouter.pa <- variogram(nonsp.pa.stres~1, snouter.sp)
plot(var.snouter.pa)

## Now let's calculate global Moran's I

snouter.pa.nb <- dnearneigh(snouter.sp@coords, 0, 18)
snouter.pa.listw <- nb2listw(snouter.pa.nb)
Glob.Moran.I.snout.pa <- moran.test(nonsp.snout.glm.stres , listw=snouter.pa.listw)
Glob.Moran.I.snout.pa$estimate
Glob.Moran.I.snout.pa$statistic

library(MASS)

group <- factor(rep("a",nrow(snouter.df)))
snouter.df <- data.frame(snouter.df, group)

#exponential correlation structure
glmm.snouter.pa.exp<- glmmPQL(snouter.pres.abs ~ prec + dist.jung, random=~1|group, data=snouter.df, 
	correlation=corExp(form=~x.coord+y.coord), family=binomial)

#spherical correlation structure
glmm.snouter.pa.exp<- glmmPQL(snouter.pres.abs ~ prec + dist.jung, random=~1|group, data=snouter.df, 
	correlation=corSpher(form=~x.coord+y.coord), family=binomial)

AIC(glmm.snouter.pa.exp, glmm.snouter.pa.exp)
## Try with Gaussian structure and let me know the result

summary(nonsp.snout.glm)

## Compute explained variance
r.squaredLR(nonsp.snout.glm)
r.squaredLR(gls.snouter.ab.exp)