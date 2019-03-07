# Plot mean fitness fitness with standard errors for G1Y13, G2Y14, and G1Y14,
# using PhenMod_G values.
# TK and AN
# 2018-09-16

# Requires: PhenMod_G.rda
# Input: PhenMod_G model objects
# Output: Plot of mean fitness with standard error bars

#setwd("~/Desktop/Dropbox/GitHubRepositories/Nashoba_Kono_Phenology/")
setwd("/Users/tomkono/Dropbox/GitHub/Nashoba_Kono_Phenology")
library(aster)
load("Results/RDA/PhenMod_G.rda")

# in order to get predicted mean fitness, we must create a data frame with hypothetical individuals
#newdata_2013 <- data.frame(
#   Initial=1,
#   Block=factor(1, levels=as.character(1:7)),
#   CrossType=factor("Nes", levels=c("Nes", "Fac")),
#   Germ=1,
#   Pods=1,
#   NumPods=1,
#   NumSeeds=1,
#   EarliestGerm=131,
#   EarliestFlowPod=210)
##### PLAY AROUND ############################
# This is the "hypothetical" data that we will use for predict(). Developed with DE, using actually observed numbers for Block, CrossType, EarliestGerm, and EarliestFlowPod.
newdata_g1y13 <- data.frame(
    Initial=1,
    #Block=as.factor(just2013$Block),
    Block=as.factor(4),
    CrossType=factor("Nes", levels=c("Nes", "Fac")),
    Germ=1,
    #Germ=just2013$Germ,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g1y13_final$EarliestGerm,
    #EarliestGerm=153,
    EarliestFlowPod=g1y13_final$EarliestFlowPod)
    #EarliestFlowPod=219)

##### END PLAY AROUND ############################

newdata_g2y14 <- data.frame(
    Initial=1,
    #Block=as.factor(just2013p$Block),
    Block=as.factor(4),
    Germ=1,
    #Germ=median(just2013p$Germ),
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g2y14_final$EarliestGerm,
    #EarliestGerm=140,
    EarliestFlowPod=g2y14_final$EarliestFlowPod)
    #EarliestFlowPod=223)
    #EarliestFlowPod=209)


newdata_g1y14 <- data.frame(
    Initial=1,
    #Block=as.factor(just2014$Block),
    Block=as.factor(4),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g1y14_final$EarliestGerm,
    #EarliestGerm=167,
    EarliestFlowPod=g1y14_final$EarliestFlowPod)
    #EarliestFlowPod=213)

# reshape newdata
renewdata_g1y13 <- reshape(
    newdata_g1y13,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
renewdata_g2y14 <- reshape(
    newdata_g2y14,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
renewdata_g1y14 <- reshape(
    newdata_g1y14,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
# Put the 'fit' variable onto renewdata - the row in which varb==NumSeeds will be 1, others will be 0
fit_g1y13 <- as.numeric(grepl("NumSeed", renewdata_g1y13$varb))
renewdata_g1y13$fit <- fit_g1y13
fit_g2y14 <- as.numeric(grepl("NumSeed", renewdata_g2y14$varb))
renewdata_g2y14$fit <- fit_g2y14
fit_g1y14 <- as.numeric(grepl("NumSeed", renewdata_g1y14$varb))
renewdata_g1y14$fit <- fit_g1y14

#   Since each cohort has a different "hypothetical" data frame to use for prediction, each cohort should have its own amat.
nind_g1y13 <- nrow(newdata_g1y13)
nnode_g1y13 <- length(vars)
amat_g1y13 <- array(0, c(nind_g1y13, nnode_g1y13, nind_g1y13))
### This is our version of the "for (i in 1:nind) amat[i, grep("hdct", vars), i] <- 1" line TR644, page 53.
### Instead of using a for loop over the number of individuals, we just use 1, since we are making only hypothetical individual.
### Instead of using grep("hdct", vars), we are using just 4, since we only have one year of our fitness proxy (hdct was over four years)
for(i in 1:nrow(newdata_g1y13)) amat_g1y13[i, 4, i] <- 1

nind_g2y14 <- nrow(newdata_g2y14)
nnode_g2y14 <- length(vars)
amat_g2y14 <- array(0, c(nind_g2y14, nnode_g2y14, nind_g2y14))
for(i in 1:nrow(newdata_g2y14)) amat_g2y14[i, 4, i] <- 1

nind_g1y14 <- nrow(newdata_g1y14)
nnode_g1y14 <- length(vars)
amat_g1y14 <- array(0, c(nind_g1y14, nnode_g1y14, nind_g1y14))
for(i in 1:nrow(newdata_g1y14)) amat_g1y14[i, 4, i] <- 1

# Now we predict. We will do separate predictions for 2013, 2013p, and 2014, using their separate model outputs
pout_g1y13 <- predict(
    g1y13_aoutSeedG,
    varvar=varb,
    idvar=id,
    root=Initial,
    newdata=renewdata_g1y13,
    se.fit=TRUE,
    amat=amat_g1y13,
    parm.type="mean.value",
    info.tol=1e-15)
    
pout_g2y14 <- predict(
    g2y14_aoutSeedG,
    varvar=varb,
    idvar=id,
    root=Initial,
    newdata=renewdata_g2y14,
    se.fit=TRUE,
    amat=amat_g2y14,
    parm.type="mean.value",
    info.tol=1e-15)

pout_g1y14 <- predict(
    g1y14_aoutSeedG,
    varvar=varb,
    idvar=id,
    root=Initial,
    newdata=renewdata_g1y14,
    se.fit=TRUE,
    amat=amat_g1y14,
    parm.type="mean.value",
    info.tol=1e-10)

# Now, we get means and se for each cohort
## mean predicted value (D.E.)
wbar_g1y13sampEst <- mean(pout_g1y13$fit)
wbar_g1y13 <- ((wbar_g1y13sampEst)*2)
#[1] 5.544531
## variance of the mean predicted value (D.E.)
grad_g1y13 <- pout_g1y13$gradient
var_g1y13 <- g1y13_aoutSeedG$fisher  
predvar_g1y13 <- grad_g1y13 %*% solve(var_g1y13) %*% t(grad_g1y13)
se_g1y13 <- t(rep(1/nind_g1y13,nind_g1y13)) %*% predvar_g1y13 %*% rep(1/nind_g1y13,nind_g1y13) 
se_g1y13 <- sqrt(se_g1y13)
#           [,1]
#[1,] 0.03917693

wbar_g2y14 <- mean(pout_g2y14$fit)
grad_g2y14 <- pout_g2y14$gradient
var_g2y14 <- g2y14_aoutSeedG$fisher
predvar_g2y14 <- grad_g2y14 %*% solve(var_g2y14) %*% t(grad_g2y14)
se_g2y14 <- t(rep(1/nind_g2y14,nind_g2y14)) %*% predvar_g2y14 %*% rep(1/nind_g2y14,nind_g2y14)
se_g2y14 <- sqrt(se_g2y14)

wbar_g1y14 <- mean(pout_g1y14$fit)
grad_g1y14 <- pout_g1y14$gradient
var_g1y14 <- g1y14_aoutSeedG$fisher
predvar_g1y14 <- grad_g1y14 %*% solve(var_g1y14) %*% t(grad_g1y14)
se_g1y14 <- t(rep(1/nind_g1y14,nind_g1y14)) %*% predvar_g1y14 %*% rep(1/nind_g1y14,nind_g1y14)
se_g1y14 <- sqrt(se_g1y14)


print(wbar_g1y13) # 3.386746
print(wbar_g2y14) # 2.934508
print(wbar_g1y14) # 7.55011

# put them on a figure
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
group_names <- c(expression(paste("G"[1], "Y"[13])), expression(paste("G"[2], "Y"[14])), expression(paste("G"[1], "Y"[14])))
i <- seq(along=group_names)
i[1] <- i[1] + 0.25
i[3] <- i[3] - 0.25
fit_values <- c(wbar_g1y13, wbar_g2y14, wbar_g1y14)
fit_se <- c(se_g1y13, se_g2y14, se_g1y14)
fit_upper_err <- fit_values + crit*fit_se
fit_lower_err <- fit_values - crit*fit_se
pdf(file="Results/Figures/MeanFitnessWithSEBars.pdf", height=3, width=3)
par(mar=c(4, 4, 0.1, 0.8))
plot(c(i, i), c(fit_upper_err, fit_lower_err), type="n", axes=FALSE, xlab="", ylab="", xlim=c(1, 3))
segments(i, fit_upper_err, i, fit_lower_err)
err_bar_width <- 0.075
segments(i-err_bar_width, fit_lower_err, i+err_bar_width, fit_lower_err)
segments(i-err_bar_width, fit_upper_err, i+err_bar_width, fit_upper_err)
segments(i-err_bar_width, fit_values, i+err_bar_width, fit_values)
axis(side=2)
title(ylab="Estimated Mean Fitness")
axis(side=1, at=i, labels=group_names)
title(xlab="Cohort")
dev.off()
