# Examine the Pearson residuals from PhenMod_G. We will be following tech
# report 661.

# Depends: PhenMod_G.R
# Input: PhenMod_G.rda
# Output: Figures with residual plots for each cohort

# 2018-07-29
# AN and TK

library(aster)

# Set the working directory to be the root of the GitHub repository
setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

# Read in PhenMod_G.rda
load("Results/RDA/PhenMod_G.rda")

# The aster model objects that we will use to check residuals are named
# g1y13_aoutSeedG, g1y14_aoutSeedG, and g2y14_aoutSeedG. We will use these
# objects in the code that starts on page 23 of the tech report.
g1y13.xi.hat <- predict(
    g1y13_aoutSeedG,
    model.type="cond",
    parm.type="mean")
g1y13.xi.hat <- matrix(
    g1y13.xi.hat,
    nrow=nrow(g1y13_aoutSeedG$x),
    ncol=ncol(g1y13_aoutSeedG$x))

# Next, we need the estimates of theta, which is the conditional canoncial
# parameter. This is from page 17 of the tech report. We also need to define a
# newdata, like we are building a fitness landscape. Note that "vars" and
# g1y13_final are being loaded from the PhenMod_G.rda file.
g1y13_newdata <- data.frame(
    Initial=1,
    Block=as.factor(4),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g1y13_final$EarliestGerm,
    EarliestFlowPod=g1y13_final$EarliestFlowPod)
g1y13_renewdata <- reshape(
    g1y13_newdata,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
g1y13_fit <- as.numeric(grepl("NumSeeds", g1y13_renewdata$varb))
g1y13_renewdata$fit <- g1y13_fit
g1y13_theta <- predict(
    g1y13_aoutSeedG,
    newdata=g1y13_renewdata,
    varvar=varb,
    idvar=id,
    root=Initial,
    model.type="conditional",
    parm.type="canonical")
g1y13_theta <- matrix(
    g1y13_theta,
    nrow=nrow(g1y13_renewdata),
    ncol=ncol(g1y13_aoutSeedG$x))

# Charlie then refers to a variable called woof. This is the terminal nodes in
# the graphical model, for all individuals that don't have a 0 in any of the
# internal nodes.
g1y13_woof <- g1y13_final$NumSeeds[g1y13_final$NumSeeds > 0]
g1y13_nwoof <- length(g1y13_woof)
# We use 4 for the column index because NumSeeds is the 4th variable of the
# graphical model.
g1y13_woof_theta <- g1y13_theta[g1y13_final$NumSeeds > 0, 4]
g1y13_woof_xi <- g1y13.xi.hat[g1y13_final$NumSeeds > 0, 4]
g1y13_wgrad <- double(g1y13_nwoof)
g1y13_winfo <- double(g1y13_nwoof)
for(i in 1:g1y13_nwoof) {
    g1y13_wgrad[i] <- famfun(fam.poisson(), deriv=1, g1y13_woof_theta[i])
    g1y13_winfo[i] <- famfun(fam.poisson(), deriv=2, g1y13_woof_theta[i])
}
all.equal(g1y13_woof_xi, g1y13_wgrad)
g1y13_pearson <- (g1y13_woof - g1y13_woof_xi)/sqrt(g1y13_winfo)

# Make the plot and save it to the proper place
pdf(file="Results/PhenFigures/G1Y13_Pearson_Residuals.pdf", height=6, width=6)
plot(
    g1y13_pearson~g1y13_woof_xi,
    xlab="Fitted Values",
    ylab="Pearson Residuals",
    main="G1Y13 Residuals")
abline(h=0, lty=2, lwd=2)
dev.off()
pdf(file="Results/PhenFigures/G1Y13_Pearson_Residuals_Hist.pdf", height=6, width=6)
hist(g1y13_pearson,
    xlab="Pearson Residuals",
    main="G1Y13 Residuals",
    breaks=20,
    freq=FALSE)
x <- seq(-200, 200, length.out=1000)
y <- dnorm(x, mean=mean(g1y13_pearson), sd=sd(g1y13_pearson))
lines(x, y, col="blue", lwd=2)
dev.off()

# Do the same for G1Y14
g1y14.xi.hat <- predict(
    g1y14_aoutSeedG,
    model.type="cond",
    parm.type="mean")
g1y14.xi.hat <- matrix(
    g1y14.xi.hat,
    nrow=nrow(g1y14_aoutSeedG$x),
    ncol=ncol(g1y14_aoutSeedG$x))
# Next, we need the estimates of theta, which is the conditional canoncial
# parameter. This is from page 17 of the tech report. We also need to define a
# newdata, like we are building a fitness landscape.
g1y14_newdata <- data.frame(
    Initial=1,
    Block=as.factor(4),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g1y14_final$EarliestGerm,
    EarliestFlowPod=g1y14_final$EarliestFlowPod)
g1y14_renewdata <- reshape(
    g1y14_newdata,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
g1y14_fit <- as.numeric(grepl("NumSeeds", g1y14_renewdata$varb))
g1y14_renewdata$fit <- g1y14_fit
g1y14_theta <- predict(
    g1y14_aoutSeedG,
    newdata=g1y14_renewdata,
    varvar=varb,
    idvar=id,
    root=Initial,
    model.type="conditional",
    parm.type="canonical")
g1y14_theta <- matrix(
    g1y14_theta,
    nrow=nrow(g1y14_renewdata),
    ncol=ncol(g1y14_aoutSeedG$x))

# Charlie then refers to a variable called woof. This is the terminal nodes in
# the graphical model, for all individuals that don't have a 0 in any of the
# internal nodes.
g1y14_woof <- g1y14_final$NumSeeds[g1y14_final$NumSeeds > 0]
g1y14_nwoof <- length(g1y14_woof)
# We use 4 for the column index because NumSeeds is the 4th variable of the
# graphical model.
g1y14_woof_theta <- g1y14_theta[g1y14_final$NumSeeds > 0, 4]
g1y14_woof_xi <- g1y14.xi.hat[g1y14_final$NumSeeds > 0, 4]
g1y14_wgrad <- double(g1y14_nwoof)
g1y14_winfo <- double(g1y14_nwoof)
for(i in 1:g1y14_nwoof) {
    g1y14_wgrad[i] <- famfun(fam.poisson(), deriv=1, g1y14_woof_theta[i])
    g1y14_winfo[i] <- famfun(fam.poisson(), deriv=2, g1y14_woof_theta[i])
}
all.equal(g1y14_woof_xi, g1y14_wgrad)
g1y14_pearson <- (g1y14_woof - g1y14_woof_xi)/sqrt(g1y14_winfo)
pdf(file="Results/PhenFigures/G1Y14_Pearson_Residuals.pdf", height=6, width=6)
plot(
    g1y14_pearson~g1y14_woof_xi,
    xlab="Fitted Values",
    ylab="Pearson Residuals",
    main="G1Y14 Residuals")
abline(h=0, lty=2, lwd=2)
dev.off()
pdf(file="Results/PhenFigures/G1Y14_Pearson_Residuals_Hist.pdf", height=6, width=6)
hist(g1y14_pearson,
    xlab="Pearson Residuals",
    main="G1Y14 Residuals",
    breaks=20,,
    freq=FALSE)
x <- seq(-200, 200, length.out=1000)
y <- dnorm(x, mean=mean(g1y14_pearson), sd=sd(g1y14_pearson))
lines(x, y, col="blue", lwd=2)
dev.off()

# And G2Y14
g2y14.xi.hat <- predict(
    g2y14_aoutSeedG,
    model.type="cond",
    parm.type="mean")
g2y14.xi.hat <- matrix(
    g2y14.xi.hat,
    nrow=nrow(g2y14_aoutSeedG$x),
    ncol=ncol(g2y14_aoutSeedG$x))
# Next, we need the estimates of theta, which is the conditional canoncial
# parameter. This is from page 17 of the tech report. We also need to define a
# newdata, like we are building a fitness landscape.
g2y14_newdata <- data.frame(
    Initial=1,
    Block=as.factor(4),
    CrossType=factor("Nes", levels=c("Nes", "Fac")),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g2y14_final$EarliestGerm,
    EarliestFlowPod=g2y14_final$EarliestFlowPod)
g2y14_renewdata <- reshape(
    g2y14_newdata,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
g2y14_fit <- as.numeric(grepl("NumSeeds", g2y14_renewdata$varb))
g2y14_renewdata$fit <- g2y14_fit
g2y14_theta <- predict(
    g2y14_aoutSeedG,
    newdata=g2y14_renewdata,
    varvar=varb,
    idvar=id,
    root=Initial,
    model.type="conditional",
    parm.type="canonical")
g2y14_theta <- matrix(
    g2y14_theta,
    nrow=nrow(g2y14_renewdata),
    ncol=ncol(g2y14_aoutSeedG$x))

# Charlie then refers to a variable called woof. This is the terminal nodes in
# the graphical model, for all individuals that don't have a 0 in any of the
# internal nodes.
g2y14_woof <- g2y14_final$NumSeeds[g2y14_final$NumSeeds > 0]
g2y14_nwoof <- length(g2y14_woof)
# We use 4 for the column index because NumSeeds is the 4th variable of the
# graphical model.
g2y14_woof_theta <- g2y14_theta[g2y14_final$NumSeeds > 0, 4]
g2y14_woof_xi <- g2y14.xi.hat[g2y14_final$NumSeeds > 0, 4]
g2y14_wgrad <- double(g2y14_nwoof)
g2y14_winfo <- double(g2y14_nwoof)
for(i in 1:g2y14_nwoof) {
    g2y14_wgrad[i] <- famfun(fam.poisson(), deriv=1, g2y14_woof_theta[i])
    g2y14_winfo[i] <- famfun(fam.poisson(), deriv=2, g2y14_woof_theta[i])
}
all.equal(g2y14_woof_xi, g2y14_wgrad)
g2y14_pearson <- (g2y14_woof - g2y14_woof_xi)/sqrt(g2y14_winfo)
pdf(file="Results/PhenFigures/G2Y14_Pearson_Residuals.pdf", height=6, width=6)
plot(
    g2y14_pearson~g2y14_woof_xi,
    xlab="Fitted Values",
    ylab="Pearson Residuals",
    main="G2Y14 Residuals")
abline(h=0, lty=2, lwd=2)
dev.off()
pdf(file="Results/PhenFigures/G2Y14_Pearson_Residuals_Hist.pdf", height=6, width=6)
hist(g2y14_pearson,
    xlab="Pearson Residuals",
    main="G2Y14 Residuals",
    breaks=20,
    freq=FALSE)
x <- seq(-1000, 1000, length.out=1000)
y <- dnorm(x, mean=mean(g2y14_pearson), sd=sd(g2y14_pearson))
lines(x, y, col="blue", lwd=2)
dev.off()
