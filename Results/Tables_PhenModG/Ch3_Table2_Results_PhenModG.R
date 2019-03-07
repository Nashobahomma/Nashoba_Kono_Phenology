# R code to make chapter 3 table 2, the first table in the Results section

# Depends: PhenMod_G
# Input: PhenMod_G.rda
# Output: Means of EG, EF for each cohort. Counts of Locs with EG, EF. Total NumSeeds for G2Y14.

# Set working directory
setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

# Load the RDA
load("Results/RDA/PhenMod_G.rda")

# The mean of EarliestGerm in G1Y13
mean(g1y13_trimmed$EarliestGerm, na.rm=TRUE)
# 135.6531

# The number of Locs with an observation of EG
sum(!is.na(g1y13_trimmed$EarliestGerm))
# 320

# The mean of EarliestFlowPod in G1Y13
mean(g1y13_trimmed$EarliestFlowPod, na.rm=TRUE)
# 213.9382

# The number of Locs with an observation of EF
sum(!is.na(g1y13_trimmed$EarliestFlowPod))
# 275


# The mean of EarliestGerm in G1Y14
mean(g1y14_trimmed$EarliestGerm, na.rm=TRUE)
# 161.987

# The number of Locs with an observation of EG
sum(!is.na(g1y14_trimmed$EarliestGerm))
# 230

# The mean of EarliestFlowPod in G1Y14
mean(g1y14_trimmed$EarliestFlowPod, na.rm=TRUE)
# 216.9583

# The number of Locs with an observation of EF
sum(!is.na(g1y14_trimmed$EarliestFlowPod))
# 120


# The mean of EarliestGerm in G2Y14
mean(g2y14_trimmed$EarliestGerm, na.rm=TRUE)
# 153.8333

# The number of Locs with an observation of EG
sum(!is.na(g2y14_trimmed$EarliestGerm))
# 114

# The mean of EarliestFlowPod in G2Y14
mean(g2y14_trimmed$EarliestFlowPod, na.rm=TRUE)
# 214.8361

# The number of Locs with an observation of EF
sum(!is.na(g2y14_trimmed$EarliestFlowPod))
# 61

# The tota number of seeds from G2Y14
sum(g2y14_trimmed$NumSeeds)
# 6352
