# Run Aster analyses using EarliestGerminant (EG) and EarliestFlowPod (EF)
# as predictor variables. The graphical model in use is:
#   Initial -> Germ -> Pods -> NumPods -> NumSeeds
# We will fit aster models to each cohort of G1Y13, G1Y14, and G2Y14.

# Depends: None
# Input: Model-ready data CSVs
# Output: PhenMod_G.rda, with aster objects and cleaned data frames

# Changes from PhenMod_F:
# - Fix an issue where G2Y14 numbers planted may be higher than G1Y13 numbers of
#   seed produced. To fix this, we will increase the G1Y13 NumSeeds to be at
#   least as high as G2Y14 number planted.
# - Be much more clear about where figures and data for tables get written.
# - Include statements that save model summaries to files, so they stay
#   consistent for lookup when writing the manuscript text.
# - Do not test for LateGerm and CrossType, because they were not significant
#   in PhenMod_F

# 2018-07-26
# AN and TK

library(aster)

# Set the working directory to be the base of the GitHub repository
setwd("/Users/tomkono/Dropbox/GitHub/Nashoba_Kono_Phenology")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Nashoba_Kono_Phenology")

# Set up the graphical model. We will uset he same graphical model for each
# cohort, so it does not need to change.
vars <- c("Germ", "Pods", "NumPods", "NumSeeds")
pred <- c(0, 1, 2, 3)
# The families for the functions are labeled as follows:
# Bernoulli (1), Bernoulli (1), 0-Poisson (3), Poisson (2)
fam <- c(1, 1, 3, 2)

##### G1Y13 #####
# Read in the G1Y13 data.
g1y13_wide_data <- read.csv("Data/G1Y13_Data.csv", header=TRUE)

# Now, trim down the datasheet to just the columns needed for the Aster model.
g1y13_trimmed <- g1y13_wide_data[,c("Loc", "ID", "Pat", "Mat", "Block", "Initial",
                                "Germ", "Pods", "NumPods", "NumSeeds",
                                "EarliestGerm", "EarliestFlowPod", "CrossType")]

# Note here: 2.1.11 has a problem where 5 seeds were planted in 2014, but there
# were no seeds collected in 2013. A note in the 2014 MasterA sheet says that
# the seeds probably belong to 2.2.1
g1y13_trimmed$NumSeeds[g1y13_trimmed$Loc == "2.1.11"] <- g1y13_trimmed$NumSeeds[g1y13_trimmed$Loc == "2.1.11"] - 5
g1y13_trimmed$NumSeeds[g1y13_trimmed$Loc == "2.2.1"] <- g1y13_trimmed$NumSeeds[g1y13_trimmed$Loc == "2.2.1"] + 5

# 26 August 2018
# AN and TK: 2.2.3 for G1Y13 should be set to 0 NumSeeds. 2.2.3 in G2Y14
# does not have any data beyond seeds planted, so it is likely a mis-entry from
# a different cohort.
g1y13_trimmed$NumSeeds[g1y13_trimmed$Loc == "2.2.3"] <- 0
# One last QC step: keep only rows that have no NAs in EG and EF so that we
# can fit an Aster model with predictors
keep <- !is.na(g1y13_trimmed$EarliestGerm) & !is.na(g1y13_trimmed$EarliestFlowPod)
g1y13_final <- g1y13_trimmed[keep,]

# Reshape the wide data to long format.
g1y13_long_data <- reshape(
    g1y13_final,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")

# Set "NumSeeds" as the fitness variable, and add the "fit" column
g1y13_fit <- as.numeric(g1y13_long_data$varb == "NumSeeds")
g1y13_long_data$fit <- g1y13_fit

# And fit the model
g1y13_aoutSeedG <- aster(
    resp~varb+fit:(as.factor(Block)+EarliestGerm+EarliestFlowPod),
    pred, fam, varb, id, Initial, data=g1y13_long_data)

# Write the summary to a text file
sink("Results/Model_Outputs/G1Y13_PhenMod_G_Summary.txt")
summary(g1y13_aoutSeedG, info.tol=1e-09)
sink()

##### G1Y14 #####
# Read in the G1Y14 data
g1y14_wide_data <- read.csv("Data/G1Y14_Data.csv",header=TRUE)
# Trim it down, just like the G1Y13 data
g1y14_trimmed <- g1y14_wide_data[, c("Loc", "Block", "Initial", "Germ", "Pods",
                                  "Mat", "MatL", "MatDam", "MatSire", "Pat",
                                  "PatL", "PatDam", "PatSire", "NumPods",
                                  "NumSeeds", "EarliestGerm", "EarliestFlowPod")]
# Keep only rows with no NA in EG and EF
keep <- !is.na(g1y14_trimmed$EarliestGerm) & !is.na(g1y14_trimmed$EarliestFlowPod)
g1y14_final <- g1y14_trimmed[keep,]

# Reshape and add the fit column
g1y14_long_data <- reshape(
    g1y14_final,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
g1y14_fit <- as.numeric(g1y14_long_data$varb == "NumSeeds")
g1y14_long_data$fit <- g1y14_fit

# Fit the model
g1y14_aoutSeedG <- aster(
    resp~varb+fit:(as.factor(Block)+EarliestGerm+EarliestFlowPod),
    pred, fam, varb, id, Initial, data=g1y14_long_data)

# Write the summary to a text file
sink("Results/Model_Outputs/G1Y14_PhenMod_G_Summary.txt")
summary(g1y14_aoutSeedG, info.tol=1e-10)
sink()

##### G2Y14 #####
# Read the G2Y14 data
g2y14_wide_data <- read.csv("Data/G2Y14_Data.csv",header=TRUE)

# Trim columns for the model
g2y14_trimmed <- g2y14_wide_data[, c("Loc2013p", "Block", "Initial", "Germ",
                                     "Pods", "NumPods", "NumSeeds",
                                     "EarliestGerm", "EarliestFlowPod")]

# We rename the columns of the data, since we want to change the Loc2013p to
# just Loc.
names(g2y14_trimmed) <- c("Loc", "Block", "Initial", "Germ", "Pods", "NumPods",
                          "NumSeeds", "EarliestGerm", "EarliestFlowPod")

# Remove rows with NA in EarliestGerm or EarliestFlowPod
keep <- !is.na(g2y14_trimmed$EarliestGerm) & !is.na(g2y14_trimmed$EarliestFlowPod)
g2y14_final <- g2y14_trimmed[keep, ]

# Reshape and add the "fit" variable
g2y14_long_data <- reshape(
    g2y14_final,
    varying=list(vars),
    direction="long",
    timevar="varb",
    times=as.factor(vars),
    v.names="resp")
g2y14_fit <- as.numeric(g2y14_long_data$varb == "NumSeeds")
g2y14_long_data$fit <- g2y14_fit

# Fit the model
g2y14_aoutSeedG <- aster(
    resp~varb+fit:(as.factor(Block)+EarliestGerm+EarliestFlowPod),
    pred, fam, varb, id, Initial, data=g2y14_long_data)

# Write the summary to a text file
sink("Results/Model_Outputs/G2Y14_PhenMod_G_Summary.txt")
summary(g2y14_aoutSeedG, info.tol=1e-10)
sink()

# Save the RDA file for later use
save(
    g1y13_aoutSeedG,
    g1y14_aoutSeedG,
    g2y14_aoutSeedG,
    g1y13_trimmed,
    g1y14_trimmed,
    g2y14_trimmed,
    g1y13_final,
    g1y14_final,
    g2y14_final,
    vars,
    pred,
    fam,
    file="Results/RDA/PhenMod_G.rda")
