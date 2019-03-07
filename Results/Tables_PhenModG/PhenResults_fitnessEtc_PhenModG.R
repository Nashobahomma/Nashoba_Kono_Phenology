# Code to calculate numbers for Table 1:
# Number of seeds planted, proportion germinated, number of seeds produced,
# and mean fitness of each cohort. We will use the dataset that has not yet
# been filtered for EG and EF missingness.

# Depends: PhenMod_G
# Input: PhenMod_G.rda
# Output: Numbers for table 1 (PhenResults_fitnessEtc)

# Set working directory
setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

library(aster)
load("Results/RDA/PhenMod_G.rda")

# The total number of seeds planted for G1Y13
sum(g1y13_trimmed$Initial)
# 2370

# The total number of seeds germinated
sum(g1y13_trimmed$Germ)
# 694

# The proportion of seeds germinated
sum(g1y13_trimmed$Germ) / sum(g1y13_trimmed$Initial)
# 0.292827

# The tota number of seeds produced (*2 for sampling)
sum(g1y13_trimmed$NumSeeds) * 2
# 7088

# Standard deviation of pre-sampling NumSeeds
sd(g1y13_trimmed$NumSeeds)
# 25.30496

# We have to use predict() to get mean fitness for the cohort
#   First, make the "hypothetical data"
newdata_g1y13 <- data.frame(
    Initial=1,
    Block=as.factor(4),
    CrossType=factor("Nes", levels=c("Nes", "Fac")),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g1y13_final$EarliestGerm,
    EarliestFlowPod=g1y13_final$EarliestFlowPod)
#   Then reshape it to "long" format
renewdata_g1y13 <- reshape(
    newdata_g1y13,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
#   Make the "fitness" variable
fit_g1y13 <- as.numeric(grepl("NumSeed", renewdata_g1y13$varb))
renewdata_g1y13$fit <- fit_g1y13
#   Make the amat for prediction
nind_g1y13 <- nrow(newdata_g1y13)
nnode_g1y13 <- length(vars)
amat_g1y13 <- array(0, c(nind_g1y13, nnode_g1y13, nind_g1y13))
for(i in 1:nrow(newdata_g1y13)) amat_g1y13[i, 4, i] <- 1
#   Use predict() to get fitness
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
#   Calculate the mean fitness (*2 for sampling). This is Wbar!
mean(pout_g1y13$fit) * 2
# 3.386746


# The total number of seeds planted for G1Y14
sum(g1y14_trimmed$Initial)
# 5193

# The total number of seeds germinated
sum(g1y14_trimmed$Germ)
# 378

# The proportion of seeds germinated
sum(g1y14_trimmed$Germ) / sum(g1y14_trimmed$Initial)
# 0.07279029

# The tota number of seeds produced (*2 for sampling)
sum(g1y14_trimmed$NumSeeds)
# 4207

# Standard deviation of G1Y14 NumSeeds
sd(g1y14_trimmed$NumSeeds)
# 18.47093

# We have to use predict() to get mean fitness for the cohort
#   First, make the "hypothetical data"
newdata_g1y14 <- data.frame(
    Initial=1,
    Block=as.factor(4),
    CrossType=factor("Nes", levels=c("Nes", "Fac")),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g1y14_final$EarliestGerm,
    EarliestFlowPod=g1y14_final$EarliestFlowPod)
#   Then reshape it to "long" format
renewdata_g1y14 <- reshape(
    newdata_g1y14,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
#   Make the "fitness" variable
fit_g1y14 <- as.numeric(grepl("NumSeed", renewdata_g1y14$varb))
renewdata_g1y14$fit <- fit_g1y14
#   Make the amat for prediction
nind_g1y14 <- nrow(newdata_g1y14)
nnode_g1y14 <- length(vars)
amat_g1y14 <- array(0, c(nind_g1y14, nnode_g1y14, nind_g1y14))
for(i in 1:nrow(newdata_g1y14)) amat_g1y14[i, 4, i] <- 1
#   Use predict() to get fitness
pout_g1y14 <- predict(
    g1y14_aoutSeedG,
    varvar=varb,
    idvar=id,
    root=Initial,
    newdata=renewdata_g1y14,
    se.fit=TRUE,
    amat=amat_g1y14,
    parm.type="mean.value",
    info.tol=1e-15)
#   Calculate the mean fitness. This is Wbar!
mean(pout_g1y14$fit)
# 7.55011


# The total number of seeds planted for G2Y14
sum(g2y14_trimmed$Initial)
# 6842

# The total number of seeds germinated
sum(g2y14_trimmed$Germ)
# 950

# The proportion of seeds germinated
sum(g2y14_trimmed$Germ) / sum(g2y14_trimmed$Initial)
# 0.1388483

# The total number of seeds produced
sum(g2y14_trimmed$NumSeeds)
# 6352

# The standard deviation of G2Y14 NumSeeds
sd(g2y14_trimmed$NumSeeds)
# 157.2211

# We have to use predict() to get mean fitness for the cohort
#   First, make the "hypothetical data"
newdata_g2y14 <- data.frame(
    Initial=1,
    Block=as.factor(4),
    CrossType=factor("Nes", levels=c("Nes", "Fac")),
    Germ=1,
    Pods=1,
    NumPods=1,
    NumSeeds=1,
    EarliestGerm=g2y14_final$EarliestGerm,
    EarliestFlowPod=g2y14_final$EarliestFlowPod)
#   Then reshape it to "long" format
renewdata_g2y14 <- reshape(
    newdata_g2y14,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
#   Make the "fitness" variable
fit_g2y14 <- as.numeric(grepl("NumSeed", renewdata_g2y14$varb))
renewdata_g2y14$fit <- fit_g2y14
#   Make the amat for prediction
nind_g2y14 <- nrow(newdata_g2y14)
nnode_g2y14 <- length(vars)
amat_g2y14 <- array(0, c(nind_g2y14, nnode_g2y14, nind_g2y14))
for(i in 1:nrow(newdata_g2y14)) amat_g2y14[i, 4, i] <- 1
#   Use predict() to get fitness
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
#   Calculate the mean fitness. This is Wbar!
mean(pout_g2y14$fit)
# 2.934508
