# Calculate beta (selection gradient) as in Lande and Arnold analysis using
# fitness values from aster models. We will use this to predict the response to
# selection, then compare the prediction to the observed response.
# TK and AN
# 2018-08-12

# Depends: PhenMod_G.R
# Input: PhenMod_G.rda
# Output: Estimates of beta for EG and EF, for input into LA_pred_obs.tex

setwd("/Users/tomkono/Dropbox/GitHub/Nashoba_Kono_Phenology")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Nashoba_Kono_Phenology")
library(aster)

# Load the RDA file, which has the aster model objects.
load("Results/RDA/PhenMod_G.rda")

# We will estimate Lande and Arnold betas for G1Y13 first. The aster object
# is called g1y13_aoutSeedG. Use predict() to get maximum-likelihood
# estimates of lifetime fitness for each individual.
g1y13_pout <- predict(g1y13_aoutSeedG)
g1y13_pout <- matrix(
    g1y13_pout,
    nrow=nrow(g1y13_aoutSeedG$x),
    ncol=ncol(g1y13_aoutSeedG$x)
    )
# Set the column names of the pout matrix using the names of the life history
# stages.
colnames(g1y13_pout) <- colnames(g1y13_aoutSeedG$x)
# We will use "NumSeeds" as the surrogate for fitness.
g1y13_mufit <- g1y13_pout[, "NumSeeds"]
# For G1Y13, there was ~50% sampling imposed, so multiply by 2
g1y13_mufit <- g1y13_mufit * 2
# Translate the absolute fitness estimates into relative fitness estimates by
# dividing by the mean
g1y13_relfit <- g1y13_mufit/mean(g1y13_mufit)

# Center the phenotypic predictors to 0, too
g1y13_eg_cent <- g1y13_final$EarliestGerm - mean(g1y13_final$EarliestGerm)
g1y13_ef_cent <- g1y13_final$EarliestFlowPod - mean(g1y13_final$EarliestFlowPod)

# This implements the Lande-Arnold OLS method to get estimates of beta. They are
# the coefficients of the centered predictors in the regression
g1y13_ols <- lm(g1y13_relfit ~ 0 + g1y13_eg_cent + g1y13_ef_cent)
# We can plot the diagnostics of the model to examine the residuals. Perhaps
# this may explain why the betas are so different from the Aster model
pdf(file="Results/Figures/G1Y13_LA_Residuals.pdf", height=6, width=6)
plot(g1y13_ols)
dev.off()

# The betas are the coefficients of the OLS
g1y13_betas <- g1y13_ols$coefficients

# Then, make a matrix from the phenotypic predictor variables. We use it to
# calculated beta with:
#   beta = inverse(P) * cov(w, z)
# where P is the phenotypic variance-covariance matrix, w is relative fitness,
# and z is the phenotypic predictors. The phenotypic predictors will come from
# 'g1y13_final'
phen_mat <- cbind(
    g1y13_eg_cent,
    g1y13_ef_cent,
    deparse.level=0)
# Then, use the matrix of phenotypic values to make the variance-covariance
# matrix. The  %*% operator is the matrix multiplication operator
phen_var_covar <- solve(t(phen_mat) %*% phen_mat/nrow(phen_mat))
# Calculate the covariance between relative fitness and the phenotypic predictor
# variables
# To get standard errors for beta, we will have to re-predict with a different
# model matrix. We will build that next. Note: TK doesn't really know exactly
# what is being done here
mu <- predict(g1y13_aoutSeedG)
amat <- matrix(0, nrow=length(g1y13_mufit), ncol=length(mu))
blank <- matrix(0, nrow=nrow(g1y13_pout), ncol=ncol(g1y13_pout))
blank.idx <- grep("NumSeeds", colnames(g1y13_pout))
for(i in 1:nrow(amat)) {
    boom <- blank
    boom[i, blank.idx] <- 1
    amat[i,] <- as.vector(boom)
}
bmat <- phen_var_covar %*% t(phen_mat) %*% amat/nrow(phen_mat)
cmat <- apply(amat, 2, sum)/nrow(phen_mat)
cmat <- rbind(cmat)
dmat <- rbind(bmat, cmat, deparse.level=0)
d3way <- array(as.vector(t(dmat)), dim=c(dim(g1y13_aoutSeedG$modmat)[1:2], nrow(dmat)))
# We need to add info.tol=1e-9 to get standard errors here.
dout <- predict(g1y13_aoutSeedG, amat=d3way, se.fit=TRUE, info.tol=1e-9)
zeta1 <- dout$fit[1]
zeta2 <- dout$fit[2]
zeta3 <- dout$fit[3]
jacobian <- rbind(
    c(1/zeta3, 0, -zeta1/zeta3^2),
    c(0, 1/zeta3, -zeta2/zeta3^2)
    )
dvar <- dout$gradient %*% solve(g1y13_aoutSeedG$fisher) %*% t(dout$gradient)
# This is the variance-covariance matrix for the beta estimates. The diagonal
# elements give the variances; we take the square root to get the standard error
# of the beta estimates.
beta_var <- jacobian %*% dvar %*% t(jacobian)
beta_SEs <- sqrt(diag(beta_var))
# Print out the beta estimates with the standard errors
names(g1y13_betas) <- c("EG", "EF")
names(beta_SEs) <- c("EG", "EF")
print(rbind(g1y13_betas, beta_SEs))


# Do the same for G1Y14
g1y14_pout <- predict(g1y14_aoutSeedG)
g1y14_pout <- matrix(
    g1y14_pout,
    nrow=nrow(g1y14_aoutSeedG$x),
    ncol=ncol(g1y14_aoutSeedG$x)
    )
# Set the column names of the pout matrix using the names of the life history
# stages.
colnames(g1y14_pout) <- colnames(g1y14_aoutSeedG$x)
# We will use "NumSeeds" as the surrogate for fitness.
g1y14_mufit <- g1y14_pout[, "NumSeeds"]
# Translate the absolute fitness estimates into relative fitness estimates by
# dividing by the mean
g1y14_relfit <- g1y14_mufit/mean(g1y14_mufit)

# Center the phenotypic predictors to 0, too
g1y14_eg_cent <- g1y14_final$EarliestGerm - mean(g1y14_final$EarliestGerm)
g1y14_ef_cent <- g1y14_final$EarliestFlowPod - mean(g1y14_final$EarliestFlowPod)

# This implements the Lande-Arnold OLS method to get estimates of beta. They are
# the coefficients of 
g1y14_ols <- lm(g1y14_relfit ~ 0 + g1y14_eg_cent + g1y14_ef_cent)
# We can plot the diagnostics of the model to examine the residuals. Perhaps
# this may explain why the betas are so different from the Aster model
pdf(file="Results/Figures/G1Y14_LA_Residuals.pdf", height=6, width=6)
plot(g1y14_ols)
dev.off()

# The betas are the coefficients of the OLS
g1y14_betas <- g1y14_ols$coefficients

# Then, make a matrix from the phenotypic predictor variables. We use it to
# calculated beta with:
#   beta = inverse(P) * cov(w, z)
# where P is the phenotypic variance-covariance matrix, w is relative fitness,
# and z is the phenotypic predictors. The phenotypic predictors will come from
# 'g1y14_final'
phen_mat <- cbind(
    g1y14_final$EarliestGerm,
    g1y14_final$EarliestFlowPod,
    deparse.level=0)
# Then, use the matrix of phenotypic values to make the variance-covariance
# matrix. The  %*% operator is the matrix multiplication operator
phen_var_covar <- solve(t(phen_mat) %*% phen_mat/nrow(phen_mat))
# To get standard errors for beta, we will have to re-predict with a different
# model matrix. We will build that next. Note: TK doesn't really know exactly
# what is being done here
mu <- predict(g1y14_aoutSeedG)
amat <- matrix(0, nrow=length(g1y14_mufit), ncol=length(mu))
blank <- matrix(0, nrow=nrow(g1y14_pout), ncol=ncol(g1y14_pout))
blank.idx <- grep("NumSeeds", colnames(g1y14_pout))
for(i in 1:nrow(amat)) {
    boom <- blank
    boom[i, blank.idx] <- 1
    amat[i,] <- as.vector(boom)
}
bmat <- phen_var_covar %*% t(phen_mat) %*% amat/nrow(phen_mat)
cmat <- apply(amat, 2, sum)/nrow(phen_mat)
cmat <- rbind(cmat)
dmat <- rbind(bmat, cmat, deparse.level=0)
d3way <- array(as.vector(t(dmat)), dim=c(dim(g1y14_aoutSeedG$modmat)[1:2], nrow(dmat)))
# We need to add info.tol=1e-9 to get standard errors here.
dout <- predict(g1y14_aoutSeedG, amat=d3way, se.fit=TRUE, info.tol=1e-15)
zeta1 <- dout$fit[1]
zeta2 <- dout$fit[2]
zeta3 <- dout$fit[3]
jacobian <- rbind(
    c(1/zeta3, 0, -zeta1/zeta3^2),
    c(0, 1/zeta3, -zeta2/zeta3^2)
    )
dvar <- dout$gradient %*% solve(g1y14_aoutSeedG$fisher) %*% t(dout$gradient)
# This is the variance-covariance matrix for the beta estimates. The diagonal
# elements give the variances; we take the square root to get the standard error
# of the beta estimates.
beta_var <- jacobian %*% dvar %*% t(jacobian)
beta_SEs <- sqrt(diag(beta_var))
# Print out the beta estimates with the standard errors
names(g1y14_betas) <- c("EG", "EF")
names(beta_SEs) <- c("EG", "EF")
print(rbind(g1y14_betas, beta_SEs))

# Do it all again for G2Y14.
g2y14_pout <- predict(g2y14_aoutSeedG)
g2y14_pout <- matrix(
    g2y14_pout,
    nrow=nrow(g2y14_aoutSeedG$x),
    ncol=ncol(g2y14_aoutSeedG$x)
    )
# Set the column names of the pout matrix using the names of the life history
# stages.
colnames(g2y14_pout) <- colnames(g2y14_aoutSeedG$x)
# We will use "NumSeeds" as the surrogate for fitness.
g2y14_mufit <- g2y14_pout[, "NumSeeds"]
# For g2y14, there was ~50% sampling imposed, so multiply by 2
g2y14_mufit <- g2y14_mufit * 2
# Translate the absolute fitness estimates into relative fitness estimates by
# dividing by the mean
g2y14_relfit <- g2y14_mufit/mean(g2y14_mufit)

# Center the phenotypic predictors to 0, too
g2y14_eg_cent <- g2y14_final$EarliestGerm - mean(g2y14_final$EarliestGerm)
g2y14_ef_cent <- g2y14_final$EarliestFlowPod - mean(g2y14_final$EarliestFlowPod)

# This implements the Lande-Arnold OLS method to get estimates of beta. They are
# the coefficients of the centered predictors in the regression
g2y14_ols <- lm(g2y14_relfit ~ 0 + g2y14_eg_cent + g2y14_ef_cent)
# We can plot the diagnostics of the model to examine the residuals. Perhaps
# this may explain why the betas are so different from the Aster model
pdf(file="Results/Figures/G2Y14_LA_Residuals.pdf", height=6, width=6)
plot(g2y14_ols)
dev.off()

# The betas are the coefficients of the OLS
g2y14_betas <- g2y14_ols$coefficients
phen_mat <- cbind(
    g2y14_eg_cent,
    g2y14_ef_cent,
    deparse.level=0)
# Then, use the matrix of phenotypic values to make the variance-covariance
# matrix. The  %*% operator is the matrix multiplication operator
phen_var_covar <- solve(t(phen_mat) %*% phen_mat/nrow(phen_mat))
# Calculate the covariance between relative fitness and the phenotypic predictor
# variables
# To get standard errors for beta, we will have to re-predict with a different
# model matrix. We will build that next. Note: TK doesn't really know exactly
# what is being done here
mu <- predict(g2y14_aoutSeedG)
amat <- matrix(0, nrow=length(g2y14_mufit), ncol=length(mu))
blank <- matrix(0, nrow=nrow(g2y14_pout), ncol=ncol(g2y14_pout))
blank.idx <- grep("NumSeeds", colnames(g2y14_pout))
for(i in 1:nrow(amat)) {
    boom <- blank
    boom[i, blank.idx] <- 1
    amat[i,] <- as.vector(boom)
}
bmat <- phen_var_covar %*% t(phen_mat) %*% amat/nrow(phen_mat)
cmat <- apply(amat, 2, sum)/nrow(phen_mat)
cmat <- rbind(cmat)
dmat <- rbind(bmat, cmat, deparse.level=0)
d3way <- array(as.vector(t(dmat)), dim=c(dim(g2y14_aoutSeedG$modmat)[1:2], nrow(dmat)))
# We need to add info.tol=1e-9 to get standard errors here.
dout <- predict(g2y14_aoutSeedG, amat=d3way, se.fit=TRUE, info.tol=1e-9)
zeta1 <- dout$fit[1]
zeta2 <- dout$fit[2]
zeta3 <- dout$fit[3]
jacobian <- rbind(
    c(1/zeta3, 0, -zeta1/zeta3^2),
    c(0, 1/zeta3, -zeta2/zeta3^2)
    )
dvar <- dout$gradient %*% solve(g2y14_aoutSeedG$fisher) %*% t(dout$gradient)
# This is the variance-covariance matrix for the beta estimates. The diagonal
# elements give the variances; we take the square root to get the standard error
# of the beta estimates.
beta_var <- jacobian %*% dvar %*% t(jacobian)
beta_SEs <- sqrt(diag(beta_var))
# Print out the beta estimates with the standard errors
names(g2y14_betas) <- c("EG", "EF")
names(beta_SEs) <- c("EG", "EF")
print(rbind(g2y14_betas, beta_SEs))

# Compare to the OLS standard errors
summary(g1y13_ols)
summary(g1y14_ols)
summary(g2y14_ols)

# Predict the response to selection by G*beta, where G is the additive genetic
# variance-covariance matrix (G matrix) for the traits, from Quercus.
# Cross product: row-sum of the G matrix multiplied by the column-sum of the beta vector
G_matrix <- matrix(
    c(3.834964, 10.577751,
    10.577751, 15.706558),
    byrow=T,
    nrow=2
    )
g1y13_resp <- crossprod(G_matrix, g1y13_betas)
g1y14_resp <- crossprod(G_matrix, g1y14_betas)

# Print the values for the table so that it's easy to enter them into the table
tab <- matrix(
    c(
        mean(g1y13_trimmed$EarliestGerm, na.rm=TRUE),
        sd(g1y13_trimmed$EarliestGerm, na.rm=TRUE),
        mean(g1y13_trimmed$EarliestGerm, na.rm=TRUE)+g1y13_resp[1,1],
        mean(g2y14_trimmed$EarliestGerm, na.rm=TRUE),
        sd(g2y14_trimmed$EarliestGerm, na.rm=TRUE),
        mean(g1y14_trimmed$EarliestGerm, na.rm=TRUE),
        sd(g1y14_trimmed$EarliestGerm, na.rm=TRUE),
        mean(g1y14_trimmed$EarliestGerm, na.rm=TRUE)+g1y14_resp[1,1],

        mean(g1y13_trimmed$EarliestFlowPod, na.rm=TRUE),
        sd(g1y13_trimmed$EarliestFlowPod, na.rm=TRUE),
        mean(g1y13_trimmed$EarliestFlowPod, na.rm=TRUE)+g1y13_resp[2,1],
        mean(g2y14_trimmed$EarliestFlowPod, na.rm=TRUE),
        sd(g2y14_trimmed$EarliestFlowPod, na.rm=TRUE),
        mean(g1y14_trimmed$EarliestFlowPod, na.rm=TRUE),
        sd(g1y14_trimmed$EarliestFlowPod, na.rm=TRUE),
        mean(g1y14_trimmed$EarliestFlowPod, na.rm=TRUE)+g1y14_resp[2,1]),
    nrow=2,
    byrow=T)
print(round(tab, 1))
# Add column and row names to the response and print them out
colnames(g1y13_resp) <- c("G1Y13")
rownames(g1y13_resp) <- c("Resp.EG", "Resp.EF")
colnames(g1y14_resp) <- c("G1Y14")
rownames(g1y14_resp) <- c("Resp.EG", "Resp.EF")
print(g1y13_resp)
print(g1y14_resp)
