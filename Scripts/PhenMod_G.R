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
setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

# Set up the graphical model. We will uset he same graphical model for each
# cohort, so it does not need to change.
vars <- c("Germ", "Pods", "NumPods", "NumSeeds")
pred <- c(0, 1, 2, 3)
# The families for the functions are labeled as follows:
# Bernoulli (1), Bernoulli (1), 0-Poisson (3), Poisson (2)
fam <- c(1, 1, 3, 2)

##### G1Y13 #####
# Read in the G1Y13 data.
g1y13_wide_data <- read.csv("Data_Files/2013/Model_Ready/2013ModelData.csv", header=TRUE)

# Next we want to adjust the G1Y13 NumSeeds values to match the G2Y14 number
# planted. These values come from the "MasterA" data sheet for G2Y14
g2y14_masterA <- read.csv("Data_Files/2013'/Raw/scg2013'MasterA_March5.csv", header=TRUE)

# The strategy we will take is to merge the G1Y13 and G2Y14 data frames, keeping
# all Locs. Then, we will check each Loc and update the NumSeeds value as
# necessary.
g1y13_g2y14 <- merge(
    x=g1y13_wide_data, y=g2y14_masterA,
    by.x="Loc", by.y="Loc2013p",
    all.x=TRUE)

# Define a function that returns the corrected value of NumSeeds for G1Y13
correctNumSeeds <- function(data_row) {
    raw_NumSeeds <- data_row["NumSeeds"]
    num_planted <- data_row["seedsplanted"]
    # If the number planted is NA, then we cannot do any more checking, and we
    # just return the number of seeds produced by G1Y13
    if(is.na(num_planted)) {
        return(raw_NumSeeds)
    } else if(raw_NumSeeds >= num_planted) {
        return(raw_NumSeeds)
    } else if(raw_NumSeeds < num_planted) {
        return(num_planted)
    } else {
        return(NA)
    }
}

# Apply the correction function across the merged data frames, save it into
# a new column called NumSeeds_Corrected
g1y13_g2y14$NumSeeds_Corrected <- as.numeric(apply(g1y13_g2y14, 1, correctNumSeeds))

# Now, trim down the datasheet to just the columns needed for the Aster model.
g1y13_trimmed <- g1y13_g2y14[,c("Loc", "ID", "Pat", "Mat", "Block", "Initial",
                                "Germ", "Pods", "NumPods", "NumSeeds_Corrected",
                                "EarliestGerm", "EarliestFlowPod", "CrossType")]
# And then replace "NumSeeds_Corrected" with just "NumSeeds"
names(g1y13_trimmed) <- c("Loc", "ID", "Pat", "Mat", "Block", "Initial",
                                "Germ", "Pods", "NumPods", "NumSeeds",
                                "EarliestGerm", "EarliestFlowPod", "CrossType")

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

# Save the trimmed G1Y13 data frame for other analyses (FTNS, etc).
write.csv(
    g1y13_trimmed,
    file="Data_Files/2013/Model_Ready/Correct_G1Y13_NumSeeds.csv",
    quote=FALSE,
    row.names=FALSE)

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
g1y14_wide_data <- read.csv("Data_Files/2014/Model_Ready/2014_ModelData_Ver3.csv",header=TRUE)
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
g2y14_wide_data <- read.csv("Data_Files/2013'/Model_Ready/2013PrimeModelData_NonBinned_Phen.csv",header=TRUE)

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
