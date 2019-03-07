# Make fitness landscapes using PhenMod_G. 

# Depends: PhenMod_G.R
# Input: PhenMod_G.rda
# Output: Fitness landscape figures

# 2018-07-29
# AN and TK

library(aster)
library(fields)
library(graphics)

# Set the working directory to be the root of the GitHub repo
setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

load(file="Results/RDA/PhenMod_G.rda")

# Make the G1Y13 landscape
pdf(file="Results/PhenFigures/G1Y13_FitnessLandscape_PhenMod_G.pdf", height=4, width=4)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2.5, 1, 0))
 # (4) Create a grid
plotlimits <- c(125, 190, 190, 250)
num_x_points <- 150
num_y_points <- 150
x_grid_points <- seq(plotlimits[1], plotlimits[2], length=num_x_points)
y_grid_points <- seq(plotlimits[3], plotlimits[4], length=num_y_points)
# (5) calculate the outer product of the grid points. These will serve as "dummies" for EarliestGerm and EarliestFlowPod
dummy_earliestgerm <- as.vector(outer(x_grid_points, y_grid_points^0))
dummy_earliestflowpod <- as.vector(outer(x_grid_points^0, y_grid_points))
# (6) Start building the "dummy" data frame for calculating the "height" of the fitness landscape
#  First, put in Initial, EarliestGerm (dummy), and EarliestFlowPod (dummy)
dummy_widedata <- data.frame(
 Initial=rep(1, length(dummy_earliestgerm)),
 Block=factor(rep(1, length(dummy_earliestgerm)), levels=levels(as.factor(g1y13_final$Block))),
 CrossType=factor("Nes", levels=c("Nes", "Fac")),
 Germ=rep(1, length(dummy_earliestgerm)),
 Pods=rep(1, length(dummy_earliestgerm)),
 NumPods=rep(1, length(dummy_earliestgerm)),
 NumSeeds=rep(1, length(dummy_earliestgerm)),
 EarliestGerm=dummy_earliestgerm,
 EarliestFlowPod=dummy_earliestflowpod)
# (7) Reshape the dummy data to feed to aster
dummy_longdata <- reshape(dummy_widedata, varying=list(vars), direction="long", timevar="varb", times=as.factor(vars), v.names="resp")
#(8) Add the 'fit' variable to the dummy long data
#This is 0 for every row that is not 'NumSeeds' in the dummy long data, and 1 for every row that is 'NumSeeds'
#  Every individual will be treated as if they are all in the same block
dummy_fit <- as.numeric(dummy_longdata$varb=="NumSeeds")
# (9) We will not do the transform() thing because the manual says not to use it if you don't want "unanticipated consequences"
#  We instead built the entire widedata for prediction in its complete form. Still have to add fit though
dummy_longdata <- data.frame(dummy_longdata, fit=dummy_fit) 
# (10) "Predict" the fitness at each point on the grid
fit_landscape <- predict(g1y13_aoutSeedG, newdata=dummy_longdata, varvar=varb, idvar=id, root=Initial)
fit_landscape <- matrix(fit_landscape, nrow = nrow(dummy_widedata))
colnames(fit_landscape) <- vars
fit_landscape <- fit_landscape[ , "NumSeeds"]
#pbaz <- pbaz
zz <- matrix(fit_landscape, num_x_points, num_y_points)
# (11)
heatmap_colors <- gray(c(0:7 / 8))

# (8) Add the 'fit' variable to the dummy long data
#  This is 0 for every row that is not 'NumSeeds' in the dummy long data, and 1 for every row that is 'NumSeeds'
#  Every individual will be treated as if they are all in the same block
dummy_fit <- as.numeric(dummy_longdata$varb=="NumSeeds")
# (9) We will not do the transform() thing because the manual says not to use it if you don't want "unanticipated consequences"
#  We instead built the entire widedata for prediction in its complete form. Still have to add fit though
dummy_longdata <- data.frame(dummy_longdata, fit=dummy_fit) 
# (10) "Predict" the fitness at each point on the grid
fit_landscape <- predict(g1y13_aoutSeedG, newdata=dummy_longdata, varvar=varb, idvar=id, root=Initial)
fit_landscape <- matrix(fit_landscape, nrow = nrow(dummy_widedata))
colnames(fit_landscape) <- vars
fit_landscape <- fit_landscape[ , "NumSeeds"]
#pbaz <- pbaz
zz <- matrix(fit_landscape, num_x_points, num_y_points)
# (11)
heatmap_colors <- gray(c(0:7 / 8))
# Need to change this from a filled.contour to 
#    # Hopefully we only x2 for this G1Y13 group.... (for sampling etc...)
contour(
    x_grid_points,
    y_grid_points,
    (zz*2),  
    levels=seq(0, 25, 5),
    labcex=1.5,
    col=heatmap_colors,
    vfont=c("sans serif","bold"),
    method="edge",
    xlab="Earliest Germinant (Julian Date)",
    ylab="Earliest Flower/Pod (Julian Date)",
    xlim=plotlimits[1:2],
    ylim=plotlimits[3:4])

points(jitter(g1y13_final$EarliestFlowPod, 0.75) ~ jitter(g1y13_final$EarliestGerm, 0.75), cex=1, col="white", bg="black", pch=21, lwd=0.25)
           
# the direction of selection arrow. These values come from the model summary
# output. They are coefficients of fit:EG and fit:EF
g1y13_eg_sel <- -0.006626
g1y13_ef_sel <- -0.014417
# The slope of the arrow is fit:EF / fit:EG.
g1y13_slope <- g1y13_ef_sel / g1y13_eg_sel
# Define the endpoints of the arrow. We set the increase in X to be 10, so the
# increase in Y is 10*Slope
x0 <- 165
x1 <- x0+10
y0 <- 210
y1 <- y0 + (10*g1y13_slope)
arrows(x0=x0, y0=y0, x1=x1, y1=y1, angle=25, code=1, lwd=2, length=0.1)
dev.off()


# Make the G1Y14 landscape
pdf(file="Results/PhenFigures/G1Y14_FitnessLandscape_PhenMod_G.pdf", height=4, width=4)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2.5, 1, 0))
 # (4) Create a grid
plotlimits <- c(125, 190, 190, 250)
num_x_points <- 150
num_y_points <- 150
x_grid_points <- seq(plotlimits[1], plotlimits[2], length=num_x_points)
y_grid_points <- seq(plotlimits[3], plotlimits[4], length=num_y_points)
# (5) calculate the outer product of the grid points. These will serve as "dummies" for EarliestGerm and EarliestFlowPod
dummy_earliestgerm <- as.vector(outer(x_grid_points, y_grid_points^0))
dummy_earliestflowpod <- as.vector(outer(x_grid_points^0, y_grid_points))
# (6) Start building the "dummy" data frame for calculating the "height" of the fitness landscape
#  First, put in Initial, EarliestGerm (dummy), and EarliestFlowPod (dummy)
dummy_widedata <- data.frame(
 Initial=rep(1, length(dummy_earliestgerm)),
 Block=factor(rep(1, length(dummy_earliestgerm)), levels=levels(as.factor(g1y14_final$Block))),
 Germ=rep(1, length(dummy_earliestgerm)),
 Pods=rep(1, length(dummy_earliestgerm)),
 NumPods=rep(1, length(dummy_earliestgerm)),
 NumSeeds=rep(1, length(dummy_earliestgerm)),
 EarliestGerm=dummy_earliestgerm,
 EarliestFlowPod=dummy_earliestflowpod)
# (7) Reshape the dummy data to feed to aster
dummy_longdata <- reshape(dummy_widedata, varying=list(vars), direction="long", timevar="varb", times=as.factor(vars), v.names="resp")
#(8) Add the 'fit' variable to the dummy long data
#This is 0 for every row that is not 'NumSeeds' in the dummy long data, and 1 for every row that is 'NumSeeds'
#  Every individual will be treated as if they are all in the same block
dummy_fit <- as.numeric(dummy_longdata$varb=="NumSeeds")
# (9) We will not do the transform() thing because the manual says not to use it if you don't want "unanticipated consequences"
#  We instead built the entire widedata for prediction in its complete form. Still have to add fit though
dummy_longdata <- data.frame(dummy_longdata, fit=dummy_fit) 
# (10) "Predict" the fitness at each point on the grid
##
fit_landscape <- predict(g1y14_aoutSeedG, newdata=dummy_longdata, varvar=varb, idvar=id, root=Initial)
fit_landscape <- matrix(fit_landscape, nrow = nrow(dummy_widedata))
colnames(fit_landscape) <- vars
fit_landscape <- fit_landscape[ , "NumSeeds"]
#pbaz <- pbaz
zz <- matrix(fit_landscape, num_x_points, num_y_points)

# (11)
heatmap_colors <- gray(c(0:7 / 8))
# Need to change this from a filled.contour to 
#    # Don't need to multiply zz times two for the G1Y14 and G2Y14 cohorts
contour(
    x_grid_points,
    y_grid_points,
    (zz),  
    levels=seq(0, 25, 5),
    labcex=1.5,
    col=heatmap_colors,
    vfont=c("sans serif","bold"),
    method="edge",
    xlab="Earliest Germinant (Julian Date)",
    ylab="Earliest Flower/Pod (Julian Date)",
    xlim=plotlimits[1:2],
    ylim=plotlimits[3:4])

points(jitter(g1y14_final$EarliestFlowPod, 0.75) ~ jitter(g1y14_final$EarliestGerm, 0.75), cex=1, col="white", bg="black", pch=21, lwd=0.25)

# the direction of selection arrow. These values come from the model summary
# output. They are coefficients of fit:EG and fit:EF
g1y14_eg_sel <- -3.081e-04
g1y14_ef_sel <- -1.845e-03
# The slope of the arrow is fit:EF / fit:EG.
g1y14_slope <- g1y14_ef_sel / g1y14_eg_sel
# Define the endpoints of the arrow. We set the increase in X to be 3, so the
# increase in Y is 2*Slope. Change in X is arbitrary - we set it so that the
# arrow looks reasonable on the plot.
x0 <- 165
x1 <- x0+3
y0 <- 210
y1 <- y0 + (3*g1y14_slope)
arrows(x0=x0, y0=y0, x1=x1, y1=y1, angle=25, code=1, lwd=2, length=0.1)
dev.off()


# Make the G2Y14 landscape
pdf(file="Results/PhenFigures/G2Y14_FitnessLandscape_PhenMod_G.pdf", height=4, width=4)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2.5, 1, 0))
 # (4) Create a grid
plotlimits <- c(125, 190, 190, 250)
num_x_points <- 150
num_y_points <- 150
x_grid_points <- seq(plotlimits[1], plotlimits[2], length=num_x_points)
y_grid_points <- seq(plotlimits[3], plotlimits[4], length=num_y_points)
# (5) calculate the outer product of the grid points. These will serve as "dummies" for EarliestGerm and EarliestFlowPod
dummy_earliestgerm <- as.vector(outer(x_grid_points, y_grid_points^0))
dummy_earliestflowpod <- as.vector(outer(x_grid_points^0, y_grid_points))

# (6) Start building the "dummy" data frame for calculating the "height" of the fitness landscape
#  First, put in Initial, EarliestGerm (dummy), and EarliestFlowPod (dummy)
dummy_widedata <- data.frame(
 Initial=rep(1, length(dummy_earliestgerm)),
 Block=factor(rep(1, length(dummy_earliestgerm)), levels=levels(as.factor(g2y14_final$Block))),
 Germ=rep(1, length(dummy_earliestgerm)),
 Pods=rep(1, length(dummy_earliestgerm)),
 NumPods=rep(1, length(dummy_earliestgerm)),
 NumSeeds=rep(1, length(dummy_earliestgerm)),
 EarliestGerm=dummy_earliestgerm,
 EarliestFlowPod=dummy_earliestflowpod)

# (7) Reshape the dummy data to feed to aster
dummy_longdata <- reshape(dummy_widedata, varying=list(vars), direction="long", timevar="varb", times=as.factor(vars), v.names="resp")

#(8) Add the 'fit' variable to the dummy long data
#This is 0 for every row that is not 'NumSeeds' in the dummy long data, and 1 for every row that is 'NumSeeds'
#  Every individual will be treated as if they are all in the same block
dummy_fit <- as.numeric(dummy_longdata$varb=="NumSeeds")

# (9) We will not do the transform() thing because the manual says not to use it if you don't want "unanticipated consequences"
#  We instead built the entire widedata for prediction in its complete form. Still have to add fit though
dummy_longdata <- data.frame(dummy_longdata, fit=dummy_fit) 

# (10) "Predict" the fitness at each point on the grid
##
fit_landscape <- predict(g2y14_aoutSeedG, newdata=dummy_longdata, varvar=varb, idvar=id, root=Initial)
fit_landscape <- matrix(fit_landscape, nrow = nrow(dummy_widedata))
colnames(fit_landscape) <- vars
fit_landscape <- fit_landscape[ , "NumSeeds"]
#pbaz <- pbaz

zz <- matrix(fit_landscape, num_x_points, num_y_points)

# (11)
heatmap_colors <- gray(c(0:7 / 8))
#    # Don't need to multiply zz times two for the G1Y14 and G2Y14 cohorts
contour(
    x_grid_points,
    y_grid_points,
    (zz),  
    levels=seq(0, 5, 1),
    labcex=1.5,
    col=heatmap_colors,
    vfont=c("sans serif","bold"),
    method="edge",
    xlab="Earliest Germinant (Julian Date)",
    ylab="Earliest Flower/Pod (Julian Date)",
    xlim=plotlimits[1:2],
    ylim=plotlimits[3:4])

points(jitter(g2y14_final$EarliestFlowPod, 0.75) ~ jitter(g2y14_final$EarliestGerm, 0.75), cex=1, col="white", bg="black", pch=21, lwd=0.25)

# the direction of selection arrow. These values come from the model summary
# output. They are coefficients of fit:EG and fit:EF
g2y14_eg_sel <- 0.0003195
g2y14_ef_sel <- -0.0048522
# The slope of the arrow is fit:EF / fit:EG.
g2y14_slope <- g2y14_ef_sel / g2y14_eg_sel
# Define the endpoints of the arrow. We set the increase in X to be 1, so that
# the arrow is comparable to the other landscapes. Note that we also have to
# flip x0-x1 and y0-y1 because the slope is negative.
x1 <- 165
x0 <- x1+1
y1 <- 210
y0 <- y1 + (1*g2y14_slope)
arrows(x0=x0, y0=y0, x1=x1, y1=y1, angle=25, code=1, lwd=2, length=0.1)

dev.off()
