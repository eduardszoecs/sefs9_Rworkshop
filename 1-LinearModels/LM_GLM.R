#################################################################################
# Script for linear model and generalised linear model for the SEFS				#
# 				Workshop "Analysing Ecological data with R" 					#
#							by Ralf B. Schäfer									#
# 																				#
#################################################################################

# let us first set a working directory i.e. a directory where we store all files
setwd("~/Arbeit/Vortraege/2015/SEFS/Workshop")
# you have to set a path to a working directory on your local machine here
# to simplify the identification of your path, you can use the following function
file.choose()
# and select a file in your desired working directory. Subsequently, copy the path
# without (!!!) the file reference into the setwd function

# we load the data into R
data <- read.csv("Mibi.csv")

# 
str(data)
# check that everything has been loaded successfully
head(data)
# look at data
# details can be found in Voss et al. 2015 (see slides)
# this is a simplified example from the data analysis conducted in the paper
# Research question: Can we explain the microbial decomposition (k_micro) from other variables?

# first we look whether any of the explanatory variables should be transformed
# check for range
summary(data)
# most data range less than a factor of 10 from min to max
# we look closer at EC, NO2, NO3 and P

par(mfrow=c(2, 2))
plot(density(data$EC))
plot(density(data$NO2))
plot(density(data$NO3))
plot(density(data$P))

# check log transformation for NO2, NO3 and P
plot(density(data$EC))
plot(density(log10(data$NO2)))
plot(density(log10(data$NO3)))
plot(density(log10(data$P)))
# looks better

# add log transformed variables to dataframe
data$NO2log <-  log10(data$NO2)
data$NO3log <- log10(data$NO3)
data$Plog <- log10(data$P)

# we continue by looking at the collinearity between exploratory variables
pairs(data[ , c(4:8, 12,14:16)])
# displays scatterplot between variables,
# some variables exhibit high intercorrelation (EC with O2)

# code for a more informative plot for exploratory analyses was posted on the R graph gallery
# but is now offline
# the function is given as:
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    
    test <- cor.test(x,y)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.05, 0.1, 1),
                  symbols = c("*", ".", " "))
    
    text(0.5, 0.5, txt, cex = cex * r)
    text(.8, .8, Signif, cex=cex, col=2)
}
pairs(data[ , c(4:8, 12,14:16)], 
      lower.panel = panel.smooth, upper.panel = panel.cor)
# note that the significance tests should be interpreted with caution!

# a similar function is
library(car)
spm(data[ , c(4:8, 12,14:16)])

# based on ecological knowledge we rather remove EC
# however, we check the VIF when considering all variables
mod_1 <- lm(k_micro ~ Temp + stream.width + pH + O2 + EC + velo + NO2log + NO3log +Plog, data = data)
summary.lm(mod_1) 

library(car) 
# library that contains several functions for linear model diagnosis
vif(mod_1) 
# variance inflation factor, above 4 (and 5) for EC

mod_1b <- lm(k_micro ~ Temp + stream.width + pH + O2 + velo + NO2log + NO3log +Plog, data = data)
vif(mod_1b)
# now ok

# we reproduce the VIF for EC:
mod_1_v <- lm(EC ~ Temp + stream.width + pH + O2 + velo + NO2log + NO3log +Plog, data = data)
summary.lm(mod_1_v)#
1/(1-0.8388)
# gives same VIF

##############################################
# Extraction of information from the model   #
##############################################

coefficients(mod_1)
# gives the coefficients for the regression equation

confint(mod_1)
# gives the confidence intervals for the betas

fitted(mod_1)
# fitted values

str(mod_1)
# overview on all available information

##########################################
#              Modelling				 #
##########################################

##########################################
#	Comparison of all possible models    #
##########################################

# prepare data, function requires matrix of explanatory variables
library(leaps)
exp_mat <- as.matrix(data[ , c(4:7, 12,14:16)])
mod_all <- leaps(exp_mat, data$k_micro, method="adjr2", names = names(data[ , c(4:7, 12,14:16)]))
# you could also use Mellows Cp or R2
print(cbind(mod_all$adjr2, mod_all$which))
# we plot the results
par(mai = c(1.4,1.4,0.4,0.4), cex=1.5, las=1)
plot(mod_all$size-1,mod_all$adjr2, xlab = "Number of variables p-1", ylab = expression("adj. R"^2))
# highest R2 for pH, P and NO3

#######################################
#### Model Checking
######################################
mod_all_fin <- lm(k_micro ~ pH + NO3log +Plog, data = data)
par(mfrow = c(2, 2))
plot(mod_all_fin)
# model checking looks ok

# the qreference function is helpful to check whether the residuals follow a normal distribution
library(DAAG)
qreference(residuals(mod_all_fin), nrep = 8)

####################################
#	Stepwise variable selection    #
####################################

#####################################################
# Hypthesis-based I: Removing variables  				#
# where t-test does not reject H0 that beta = 0	    #
#####################################################

mod_1b <- lm(k_micro ~ Temp + stream.width + pH + O2 + velo + NO2log + NO3log +Plog, data = data)
summary(mod_1b) 
# remove O2 because it has the highest p-values

### modifying the model- removing variables where t-test does not reject H0 that beta = 0
mod_2 <- update(mod_1b, ~. -O2) 
# this removes the specified variable from the mod
summary(mod_2)
# in the next step, you would remove stream width

mod_3 <- update(mod_2, ~. -stream.width)
summary(mod_3)
## and so on

#####################################################
# Hypthesis-based II: Comparing models based on   	#
# on explained variance							    #
#####################################################

# For example, you might have constructed a priori the following
# nested models:
mod_4 <- lm(k_micro ~ Landuse_type + Temp + NO3log + Plog, data = data)
mod_5 <- lm(k_micro ~ Temp + NO3log+Plog, data = data)
mod_6 <- lm(k_micro ~ Landuse_type, data = data)
summary(mod_6)
# the estimated coefficient for factors
# are given as Intercept + estimate for factor level

# use F test to check for significant reduction in explained variance
anova(mod_4, mod_5)
# no significance
anova(mod_5, mod_6)
# does not work because same degrees of freedom
# and should not be done because models are not nested
anova(mod_4, mod_6)
# not significantly different

####################################
#Use information theoretic approach#
####################################

# 1) for stepwise model selection
# we start with model 4
drop1(mod_4)
# computes AIC for single term deletions
# add1 computes the addition of single terms
# removing Landuse_type would yield lowest AIC

# to use BIC
n = nrow(data) 
drop1(mod_4, k = log(n))
# k gives the penalty term
# although the output says AIC, it is the BIC

#2) to compare models

AIC(mod_4)
# calculation of AIC
AIC(mod_5) 
AIC(mod_6)
# AIC higher -> would suggest model 5 as best fit

# for BIC use
AIC(mod_4)
AIC(mod_5) 
AIC(mod_6)
# BIC lower -> would suggest model 4 as best fit

# automatic calculation of the corrected AIC for smaller sample sizes (AICc) can be done with the following functions
library(MuMIn)
AICc(mod_4)
AICc(mod_5)


##################################################
# automatic model building - stepwise modelling  #
##################################################

# note that automatic model selection is no panacea to 
# the problem of finding the best model

# in case of missng values, we first remove missing values 
# to avoid comparison of models with different cases
complete.cases(data)
# in case we needed to remove
data_comp <- data[complete.cases(data), ]
# would remove cases that give "FALSE"

# we first specify the full model
# omit some variables
names(data)
data_fin <- data[ ,-c(1,8:11,13)]

# including all variables is simple by using:'.'
fullmodel <- lm(k_micro ~ ., data = data_fin) 
# all variables that could be in the model

nullmodel <- lm(k_micro ~ 1, data = data_fin)
# nullmodel: no variables, only mean (intercept)

step(fullmodel, direction = "both", trace = 100, 
     scope = list(upper = fullmodel, lower = nullmodel), k = log(n))
# automatic forward and backward model building with the BIC - 
# you could start with any other model that has less variables
# Tools for explore different models visually are contained in the package meifly
# see http://had.co.nz/model-vis/ for more details

#####################################
# Relative importance of a variable #
#####################################

library(relaimpo)
pred_imp <- calc.relimp(mod_5, type = c("lmg", "first", "last", "betasq"), rela = TRUE)
# check out the help for this function for information on these relative importance measures
plot(pred_imp)

# 'lmg' is equivalent to Hierarchical partioning, which is recommended for use, together with pmvd
# Chevan, A., and Sutherland, M. (1991). Hierarchical Partitioning. Am. Stat. 45, 90–96.
# Grömping, U. (2007). Estimators of Relative Importance in Linear Regression Based on Variance Decomposition. The American Statistician 61, 139–147.

# pmvd is currently only available from the website 
# http://prof.beuth-hochschule.de/groemping/software/relaimpo-relative-importance-of-regressors/download-package-non-us-version/ 
# due to copyright restrictions

# hierarchical partitioning is also implemented in another package that allows for 
# the choice of different gof measures whereas relaimpo relies on R2
library(hier.part)
hier.part(data$k_micro, data[ , c("Temp","NO3log","Plog")], gof = "Rsqu")
# goodness of fit (gof) can be any of "RMSPE", Root-mean-square 'prediction' error
# "logLik", Log-Likelihood or "Rsqu", R-squared

#################################################
# Exercise: 									#
# Conduct a multiple regression analysis using  #
# "inv_decomp" (invertebrate decomposition) as  #
# response variable								#
#################################################


#######################
# Additional material #
#######################

####################
# Cross validation #
####################
# function provided by Kabacoff 2011:214
# prediction accuracy is validated using cross validation, execute all code below to set the function
shrinkage <- function(fit, k=10){
  library(bootstrap)
  theta.fit <- function(x, y){
    lsfit(x, y)
  }
  theta.predict <- function(fit,x){
    cbind(1, x) %*% fit$coef
  }
  x <- fit$model[ , 2:ncol(fit$model)] 
  y <- fit$model[ , 1]
  results <- crossval(x, y, theta.fit, theta.predict, ngroup = k) 
  r2 <- cor(y, fit$fitted.values)^2 
  r2cv <- cor(y, results$cv.fit)^2 
  cat("Original R-square =", r2, "\n")
  cat(k, "Fold Cross-Validated R-square =", r2cv, "\n") 
  cat("Change =", r2-r2cv, "\n") 
}
#### end of function

#calculate cross-validation for best-fit model from automatic model building
mod_bf <- lm(k_micro ~ pH, data = data)
shrinkage(mod_bf, k=2) 
# note that the sample size becomes rather small for 2 groups
# if you repeat this code, you obtain a wide range of changes in R2

# you can also use leave-one-out (prediction sum of squares) PRESS statistic
library(MPV)
PRESS(mod_bf)
PRESS(mod_bf)/deviance(mod_bf)
# deviance provides the mean sum of square error
# model is valid when both are similar 

######################################################
#Average model parameters for several best-fit models#
######################################################
library(MuMIn)
# theoretical background for model averaging can be found here: 
#
# Kenneth P. Burnham and David R. Anderson 2002: Model Selection and Multimodel Inference (2nd ed.). 
# New York: Springer-Verlag, 2002
# See also a paper of Dormann for an ecological application of this topic: 
# Dormann, C.F. et al. 2008. Prediction uncertainty of environmental change effects 
# on temperate European biodiversity. Ecology Letters 11, 235-244.
# Note also the criticism:
# Cade 2015: Model Averaging and Muddled Multimodel Inferences, Ecological monographs, in press

# we use the fullmodel (see above) as most complex model. 
# The dredge function computes all simpler models and evaluates them in terms of goodness of fit
# however, we have to set the na.fail argument
fullmodel <- lm(k_micro ~ ., data = data_fin, na.action = "na.fail") 

models <- dredge(fullmodel)
print(models)
# gives all possible models

avg_model <- model.avg(get.models(models, subset = delta < 2)) 
# get all models that are a maximum of 2 in terms of AICc higher and average model parameters across these models 
summary(avg_model)
# also variables that are not significant are included in the model now

# compute r2 for averaged model
fit_y <- avg_model$coef.shrinkage[1]+ avg_model$coef.shrinkage[2]*avg_model$x[ ,2] + avg_model$coef.shrinkage[3]*avg_model$x[ ,3] + avg_model$coef.shrinkage[4]*avg_model$x[ ,4] + avg_model$coef.shrinkage[5]*avg_model$x[ ,5]

res_var_avg <- sum((data$k_micro - fit_y)^2)
tot_var_avg <- sum((data$k_micro - mean(data$k_micro))^2)
1- res_var_avg/tot_var_avg
# not necessarily highest r2 (when compared to models above)

#################################################################################
# 						Generalised linear models in R							#
#################################################################################

## Demonstration - code modified from Maindonald 2010, chapter 8
## multiple regression for GLM
library(DAAG)
# information on data set that will be used
# Southern Corroboree frog, a critically endangered species
?frogs

# many of the procedures for the lm also apply to the glm
# they are not repeated here but subject to an exercise

head(frogs)

# construction of a glm is very similar to that of a lm
frog.glm <- glm(pres.abs ~ distance + NoOfPools+ avrain+ altitude,  family = binomial, data = frogs, na.action= "na.fail")
# na.action is only provided if later used with the dredge function
# family displays all options beside binomial including link functions
?family()

# updating and VIF work as for the lm
vif(frog.glm)
frog.glm2 <- update(frog.glm, ~ . -NoOfPools)

### so let us check the initial model fit
summary(frog.glm)
# Wald z test results can be used for single term deletions
# all other modelling functions work as for the lm outlined above
drop1(frog.glm)

# check with anova - chisquare test for deviance of both models
anova(frog.glm, frog.glm2, test = "Chisq")
# significant difference, simplification questionable
# for quasi-likelihood models, you should use test = "F" for quasi-likelihood based models
anova(frog.glm, test = "Chisq")
# gives reduction in deviance, when terms are added sequentially (Type I ANOVA)
# makes not much sense to draw inference from this, as it is sensitive to the order of the variables

# rather look at Type II ANOVA
Anova(frog.glm)
# is the same as checking whether one term can be dropped without reducing the model fit significantly:
drop1(frog.glm, test = "Chisq")
# all deletions result in a significant reduction in explained deviance

# also for the information theoretic approaches, 
# the functions are the same as for the linear model
AIC(frog.glm)
# calculation of AIC
AIC(frog.glm2) 
BIC(frog.glm)
# calculation of BIC
BIC(frog.glm2) 

### finally we look at automatic model building using the BIC
full <- glm(pres.abs ~ distance + NoOfPools + NoOfSites + avrain, family = binomial, data = frogs, na.action= "na.fail")
null <- glm(pres.abs ~ 1, family = binomial, data = frogs, na.action= "na.fail")

# calculate sample size
n <- length(frogs$pres.abs)
step(full, direction = "both", trace = 100, scope= list(upper = full, lower = null), k = log(n))
# automatic forward and backward model building with the BIC 

#########################
#GLM Model diagnostics #
#########################
# 1. let us check the dispersion of the model
summary(frog.glm)
# if we divide the Residual Deviance by the degrees of freedom we yield a dispersion parameter 
# of approximately 1 which is very good (remember that >2 or <0.5 would indicate over- or underdispersion)

frog.glm$deviance / frog.glm$df.resid
# formal calculation

# let us look at the qq plot
plot(frog.glm, which = 2)
# also relatively good fit to the qqline


## you can also look at the residuals
residualPlots(frog.glm, type = "pearson")
# the conditional mean function should be constant as we move across the plot.
# See Fox & Weisberger (2011) p. 317-320 for more details
# in case of over- or underdispersion you would refit the model with another distribution
# e.g. negative binomial distribution or if this does not help use quasi likelihood 
# estimation (e.g. for the binomial model use glm(..., family="quasibinomial")

# 2. checking for linearity
crPlots(frog.glm)
# no sign of non-linearity

## 3. check for influential observations
influenceIndexPlot(frog.glm, vars = c("Cook", "hat"), id.n = 3)
# 193 predictor outlier, 77 has highest Cook`s distance, but nevertheless nothing to worry

plot(frog.glm, which=5)
# confirms evaluation, shows that points 77 is outlier

influencePlot(frog.glm, id.n = 3)
# similar plot - cooks distance is plotted as increasing bobble 

# we can check whether removal of the points with the highest CookD results in model change:
compareCoefs(frog.glm, frogb.glm <- update(frog.glm, subset = -c(77)))
# model fit changes slightly - new estimates almost identical

# plot regarding contribution of variables
par(mfrow = c(2,2))
termplot(frog.glm)

# How to interpret the values on the y axis?
library(faraway)
ilogit(-10)
# close to 0
ilogit(5)
# 0.99

# matches approximately with value range
min(frog.glm$fitted.values)
max(frog.glm$fitted.values)

######################################################
#### cross validation and model averaging       ######
######################################################
CVbinary(frog.glm)
# gives proportion that are correctly predicted (true present, true absent)

# further techniques such as model averaging or contribution of variable are also applicable for the GLM
### model averaging
model.avg(get.models(dredge(frog.glm), subset = delta < 10))
# usually delta is set to 2, but in this case would yield only 1 model

#############################################################
# Exercise: 												#
# Conduct a Generalized multiple regression analysis using  #
# for the frog data set										#
#############################################################
