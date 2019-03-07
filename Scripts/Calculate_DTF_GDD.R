# Calculate the accumulated growing degree days (GDD) and days to flowering
# (DTF) for the C. f. phenology project. This script uses Julian dates of
# earliest germinant (EG) and earliest flower/pod (DF) and daily temperatures
# from PRISM for 2013 and 2014 at Shakopee, MN to calculate GDD. The base
# temperature is 10C, which is a value used to model soybean growth in MN.
# TK and AN

# Input: PhenMod_G.rda and PRISM daily min/max temps
# Output: Per-loc GDD and DTF values and distributions of GDD and DTF
# Requires: PhenMod_G

setwd("/Users/tomkono/Dropbox/GitHub/Nashoba_Kono_Phenology")
#setwd("/Users/amber-nashoba/Dropbox/GitHubRepositories/Nashoba_Kono_Phenology")

load("Results/RDA/PhenMod_G.rda")

################################################################################
################################################################################
################################################################################
#
# Section 1: Calculation of GDD and DTF for each Loc in each cohort
#
################################################################################
################################################################################
################################################################################

# Read the daily temperatures at SCG
temp_daily <- read.csv("Data/SCG_PRISM_Temp_Daily_MinMax.csv", header=TRUE)

#   Then, define a function that will claculate growing degree days. We use a
#   base temperature of 10C (50F), as is done for soybean.
gdd <- function(Loc, Year) {
    t_base <- 10
    #   Turn the Julian date into a month/day/year date for lookup in the
    #   PRISM daily temperatures data.
    #   First, define an "origin" - this is Julian day 0. It is Dec 31 of the
    #   previous year.
    origin <- paste(Year-1, '12-31', sep="-")
    #   Then, convert the Julian day numbers into m/d/y date format for lookup
    #   This is SUPER gross and ugly - here's a breakdown, from inside-out:
    #       - as.numeric, since Date objects just have numbers
    #       - as.Date to convert from Julian date to Y-M-D format
    #       - as.character so that we can split it and get month and day info
    #       - strsplit to separate month and day
    #       - unlist, since strsplit returns a list
    #       - as.numeric to remove leading zeroes from month and day
    germ <- as.numeric(unlist(strsplit(as.character(as.Date(as.numeric(Loc["EarliestGerm"]),origin=origin)), split="-")))
    flowpod <- as.numeric(unlist(strsplit(as.character(as.Date(as.numeric(Loc["EarliestFlowPod"]),origin=origin)), split="-")))
    #   The year is the first element in the converted dates. We need to take
    #   it from four digits down to two.
    y <- substr(germ[1], 3, 4)
    #   then make the appropriate m/d/y string to lookup
    germ <- paste(germ[2], germ[3], y, sep="/")
    flowpod <- paste(flowpod[2], flowpod[3], y, sep="/")
    #   Then, get the date range from germ to flowering. We assume the data are
    #   in order from early to late
    germ_day <- which(temp_daily$Date == germ)
    fp_day <- which(temp_daily$Date == flowpod)
    #   Then, returnt he accumulated GDDs over the period from germination to
    #   flowering.
    acc_gdds <- 0
    for(i in germ_day:fp_day) {
        mean_temp <- mean(temp_daily$T_Max[i], temp_daily$T_Min[i])
        gdd_day <- mean_temp - t_base
        acc_gdds <- acc_gdds + gdd_day
    }
    return(acc_gdds)
}

gdd_g1y13 <- apply(g1y13_final, 1, gdd, Year=2013)
gdd_g2y14 <- apply(g2y14_final, 1, gdd, Year=2014)
gdd_g1y14 <- apply(g1y14_final, 1, gdd, Year=2014)

# Append the GDDs accumulated from germ to flowering to the cohort data frames
g1y13_final <- data.frame(g1y13_final, GDD=gdd_g1y13)
g2y14_final <- data.frame(g2y14_final, GDD=gdd_g2y14)
g1y14_final <- data.frame(g1y14_final, GDD=gdd_g1y14)


# Make expressions for the G_Y_ names
name2013 <- c(expression(paste("G"[1], "Y"[13])))
name2013p <- c(expression(paste("G"[2], "Y"[14])))
name2014 <- c(expression(paste("G"[1], "Y"[14])))

days_to_maturity_g1y13 <- apply(g1y13_final, 1, function(x) { return(as.numeric(x["EarliestFlowPod"]) - as.numeric(x["EarliestGerm"]))})
days_to_maturity_g2y14 <- apply(g2y14_final, 1, function(x) { return(as.numeric(x["EarliestFlowPod"]) - as.numeric(x["EarliestGerm"]))})
days_to_maturity_g1y14 <- apply(g1y14_final, 1, function(x) { return(as.numeric(x["EarliestFlowPod"]) - as.numeric(x["EarliestGerm"]))})

# Save the data frames with accumulated GDD to flowering and DTF
g1y13_final <- data.frame(g1y13_final, DTF=days_to_maturity_g1y13)
g2y14_final <- data.frame(g2y14_final, DTF=days_to_maturity_g2y14)
g1y14_final <- data.frame(g1y14_final, DTF=days_to_maturity_g1y14)
save(g1y13_final, g2y14_final, g1y14_final, file="Results/RDA/Cohorts_with_GDDtoFlower.rda")


pdf(file="Results/PhenFigures/GDD_Distributions.pdf", height=6, width=6)
plot(
    density(gdd_g1y13),
    main="Accumulated GDD to Flowering",
    xlab="Total Accumulated GDD",
    ylab="Density",
    xlim=c(500, 1700),
    lwd=2,
    lty=1)
lines(density(gdd_g2y14), col="red", lwd=2, lty=2)
lines(density(gdd_g1y14), col="red", lwd=2, lty=1)
legend(
    "topright",
    c(name2013, name2013p, name2014),
    col=c("black", "red", "red"),
    lwd=2,
    lty=c(1, 2, 1))
dev.off()

pdf(file="Results/PhenFigures/DTF_Distributions.pdf", height=6, width=6)
plot(
    density(days_to_maturity_g1y13),
    main="Days From Germination to Flowering",
    xlab="Days to Flowering",
    ylab="Density",
    xlim=c(20, 110),
    lwd=2,
    lty=1)
lines(density(days_to_maturity_g2y14), col="red", lwd=2, lty=2)
lines(density(days_to_maturity_g1y14), col="red", lwd=2, lty=1)
legend(
    "topright",
    c(name2013, name2013p, name2014),
    col=c("black", "red", "red"),
    lwd=2,
    lty=c(1, 2, 1))
dev.off()

################################################################################
################################################################################
################################################################################
#
# Section 2: Mean GDD and DTF by maternal family
#
################################################################################
################################################################################
################################################################################

g1y13_mat <- list()
for(mat in unique(g1y13_final$Mat)) {
    #   Find where each Mat occurs in the data
    matloc <- g1y13_final$Mat == mat
    #   Get the corresponding Locs, as characters
    matloc <- as.character(g1y13_final$Loc[matloc])
    #   Remove NAs, append to the list
    g1y13_mat[[mat]] <- matloc[!is.na(matloc)]
}
# Then, for each maternal family, get the mean GDD of the Locs
mean_gdd_g1y13 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    mat_GDD_locs <- g1y13_final$Loc %in% Locs
    if(all(!mat_GDD_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    GDD <- mean(g1y13_final$GDD[mat_GDD_locs])
    return(GDD)
}
mean_dtf_g1y13 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    mat_dtf_locs <- g1y13_final$Loc %in% Locs
    if(all(!mat_dtf_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    dtf <- mean(g1y13_final$DTF[mat_dtf_locs])
    return(dtf)
}
mat_GDD_g1y13 <- lapply(g1y13_mat, mean_gdd_g1y13)
mat_DTF_g1y13 <- lapply(g1y13_mat, mean_dtf_g1y13)
# Stick 2013 into a data frame
matfam_GDD_DTF_g1y13 <- data.frame(
    Mat=names(mat_GDD_g1y13),
    GDD_G1Y13=as.numeric(mat_GDD_g1y13),
    DTF_G1Y13=as.numeric(mat_DTF_g1y13)
    )
#   Remove NA rows
matfam_GDD_DTF_g1y13 <- matfam_GDD_DTF_g1y13[!is.na(matfam_GDD_DTF_g1y13$Mat) & !is.na(matfam_GDD_DTF_g1y13$GDD) & !is.na(matfam_GDD_DTF_g1y13$DTF), ]

# Calculate G2Y14 maternal GDD means in the same way
mean_gdd_g2y14 <- function(Locs) {
    mat_GDD_locs <- g2y14_final$Loc %in% Locs
    if(all(!mat_GDD_locs)) {
        return(NA)
    }
    GDD <- mean(g2y14_final$GDD[mat_GDD_locs])
    return(GDD)
}
mean_dtf_g2y14 <- function(Locs) {
    mat_dtf_locs <- g2y14_final$Loc %in% Locs
    if(all(!mat_dtf_locs)) {
        return(NA)
    }
    dtf <- mean(g2y14_final$DTF[mat_dtf_locs])
    return(dtf)
}
mat_GDD_g2y14 <- lapply(g1y13_mat, mean_gdd_g2y14)
mat_DTF_g2y14 <- lapply(g1y13_mat, mean_dtf_g2y14)
matfam_GDD_DTF_g2y14 <- data.frame(
    Mat=names(mat_GDD_g2y14),
    GDD_G2Y14=as.numeric(mat_GDD_g2y14),
    DTF_G2Y14=as.numeric(mat_DTF_g2y14)
    )
#   Remove NA rows
matfam_GDD_DTF_g2y14 <- matfam_GDD_DTF_g2y14[!is.na(matfam_GDD_DTF_g2y14$Mat) & !is.na(matfam_GDD_DTF_g2y14$GDD) & !is.na(matfam_GDD_DTF_g2y14$DTF), ]

# Calculate the same for G1Y14
g1y14_mat <- list()
for(mat in unique(g1y14_final$Mat)) {
    #   Find where each Mat occurs in the data
    matloc <- g1y14_final$Mat == mat
    #   Get the corresponding Locs, as characters
    matloc <- as.character(g1y14_final$Loc[matloc])
    #   Remove NAs, append to the list
    g1y14_mat[[mat]] <- matloc[!is.na(matloc)]
}
mean_gdd_g1y14 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    mat_GDD_locs <- g1y14_final$Loc %in% Locs
    if(all(!mat_GDD_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    GDD <- mean(g1y14_final$GDD[mat_GDD_locs])
    return(GDD)
}
mean_dtf_g1y14 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    mat_dtf_locs <- g1y14_final$Loc %in% Locs
    if(all(!mat_dtf_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    dtf <- mean(g1y14_final$DTF[mat_dtf_locs])
    return(dtf)
}
mat_GDD_g1y14 <- lapply(g1y13_mat, mean_gdd_g1y14)
mat_DTF_g1y14 <- lapply(g1y13_mat, mean_dtf_g1y14)
matfam_GDD_DTF_g1y14 <- data.frame(
    Mat=names(mat_GDD_g1y14),
    GDD_G1Y14=as.numeric(mat_GDD_g1y14),
    DTF_G1Y14=as.numeric(mat_DTF_g1y14)
    )
#   Remove NA rows
matfam_GDD_DTF_g1y14 <- matfam_GDD_DTF_g1y14[!is.na(matfam_GDD_DTF_g1y14$Mat) & !is.na(matfam_GDD_DTF_g1y14$GDD) & !is.na(matfam_GDD_DTF_g1y14$DTF), ]


################################################################################
################################################################################
################################################################################
#
# Section 3: Mean GDD and DTF by paternal family. Note that there is no
#            paternal family information for G2Y14
#
################################################################################
################################################################################
################################################################################

g1y13_pat <- list()
for(pat in unique(g1y13_final$Pat)) {
    #   Find where each pat occurs in the data
    patloc <- g1y13_final$Pat == pat
    #   Get the corresponding Locs, as characters
    patloc <- as.character(g1y13_final$Loc[patloc])
    #   Remove NAs, append to the list
    g1y13_pat[[pat]] <- patloc[!is.na(patloc)]
}
# Then, for each paternal family, get the mean GDD of the Locs
mean_gdd_g1y13 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    pat_GDD_locs <- g1y13_final$Loc %in% Locs
    if(all(!pat_GDD_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    GDD <- mean(g1y13_final$GDD[pat_GDD_locs])
    return(GDD)
}
mean_dtf_g1y13 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    pat_dtf_locs <- g1y13_final$Loc %in% Locs
    if(all(!pat_dtf_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    dtf <- mean(g1y13_final$DTF[pat_dtf_locs])
    return(dtf)
}
pat_GDD_g1y13 <- lapply(g1y13_pat, mean_gdd_g1y13)
pat_DTF_g1y13 <- lapply(g1y13_pat, mean_dtf_g1y13)
# Stick 2013 into a data frame
patfam_GDD_DTF_g1y13 <- data.frame(
    Pat=names(pat_GDD_g1y13),
    GDD_G1Y13=as.numeric(pat_GDD_g1y13),
    DTF_G1Y13=as.numeric(pat_DTF_g1y13)
    )
#   Remove NA rows
patfam_GDD_DTF_g1y13 <- patfam_GDD_DTF_g1y13[!is.na(patfam_GDD_DTF_g1y13$Pat) & !is.na(patfam_GDD_DTF_g1y13$GDD) & !is.na(patfam_GDD_DTF_g1y13$DTF), ]

# Calculate the same for G1Y14
g1y14_pat <- list()
for(pat in unique(g1y14_final$Pat)) {
    #   Find where each pat occurs in the data
    patloc <- g1y14_final$Pat == pat
    #   Get the corresponding Locs, as characters
    patloc <- as.character(g1y14_final$Loc[patloc])
    #   Remove NAs, append to the list
    g1y14_pat[[pat]] <- patloc[!is.na(patloc)]
}
mean_gdd_g1y14 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    pat_GDD_locs <- g1y14_final$Loc %in% Locs
    if(all(!pat_GDD_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    GDD <- mean(g1y14_final$GDD[pat_GDD_locs])
    return(GDD)
}
mean_dtf_g1y14 <- function(Locs) {
    #   Get which Locs from the 2013 data are in the supplied vector of Locs
    pat_dtf_locs <- g1y14_final$Loc %in% Locs
    if(all(!pat_dtf_locs)) {
        return(NA)
    }
    #   Get the sum of the corresponding NumSeeds values
    dtf <- mean(g1y14_final$DTF[pat_dtf_locs])
    return(dtf)
}
pat_GDD_g1y14 <- lapply(g1y13_pat, mean_gdd_g1y14)
pat_DTF_g1y14 <- lapply(g1y13_pat, mean_dtf_g1y14)
patfam_GDD_DTF_g1y14 <- data.frame(
    Pat=names(pat_GDD_g1y14),
    GDD_G1Y14=as.numeric(pat_GDD_g1y14),
    DTF_G1Y14=as.numeric(pat_DTF_g1y14)
    )
#   Remove NA rows
patfam_GDD_DTF_g1y14 <- patfam_GDD_DTF_g1y14[!is.na(patfam_GDD_DTF_g1y14$Pat) & !is.na(patfam_GDD_DTF_g1y14$GDD) & !is.na(patfam_GDD_DTF_g1y14$DTF), ]


################################################################################
################################################################################
################################################################################
#
# Section 4: Print the numbers for the table
#
################################################################################
################################################################################
################################################################################

# Make a data frame so that the table prints nicely
tabledat <- data.frame(
    G1Y13=c(
        round(mean(g1y13_final$GDD, na.rm=TRUE), 1),
        round(sd(g1y13_final$GDD, na.rm=TRUE), 1),
        round(mean(matfam_GDD_DTF_g1y13$GDD_G1Y13, na.rm=TRUE), 1),
        round(sd(matfam_GDD_DTF_g1y13$GDD_G1Y13, na.rm=TRUE), 1),
        round(mean(patfam_GDD_DTF_g1y13$GDD_G1Y13, na.rm=TRUE), 1),
        round(sd(patfam_GDD_DTF_g1y13$GDD_G1Y13, na.rm=TRUE), 1),
        round(mean(g1y13_final$DTF, na.rm=TRUE), 1),
        round(sd(g1y13_final$DTF, na.rm=TRUE), 1),
        round(mean(matfam_GDD_DTF_g1y13$DTF_G1Y13, na.rm=TRUE), 1),
        round(sd(matfam_GDD_DTF_g1y13$DTF_G1Y13, na.rm=TRUE), 1),
        round(mean(patfam_GDD_DTF_g1y13$DTF_G1Y13, na.rm=TRUE), 1),
        round(sd(patfam_GDD_DTF_g1y13$DTF_G1Y13, na.rm=TRUE), 1)
        ),
    G2Y14=c(
        round(mean(g2y14_final$GDD, na.rm=TRUE), 1),
        round(sd(g2y14_final$GDD, na.rm=TRUE), 1),
        round(mean(matfam_GDD_DTF_g2y14$GDD_G2Y14, na.rm=TRUE), 1),
        round(sd(matfam_GDD_DTF_g2y14$GDD_G2Y14, na.rm=TRUE), 1),
        "NA",
        "NA",
        round(mean(g2y14_final$DTF, na.rm=TRUE), 1),
        round(sd(g2y14_final$DTF, na.rm=TRUE), 1),
        round(mean(matfam_GDD_DTF_g2y14$DTF_G2Y14, na.rm=TRUE), 1),
        round(sd(matfam_GDD_DTF_g2y14$DTF_G2Y14, na.rm=TRUE), 1),
        "NA",
        "NA"
        ),
    G1Y14=c(
        round(mean(g1y14_final$GDD, na.rm=TRUE), 1),
        round(sd(g1y14_final$GDD, na.rm=TRUE), 1),
        round(mean(matfam_GDD_DTF_g1y14$GDD_G1Y14, na.rm=TRUE), 1),
        round(sd(matfam_GDD_DTF_g1y14$GDD_G1Y14, na.rm=TRUE), 1),
        round(mean(patfam_GDD_DTF_g1y14$GDD_G1Y14, na.rm=TRUE), 1),
        round(sd(patfam_GDD_DTF_g1y14$GDD_G1Y14, na.rm=TRUE), 1),
        round(mean(g1y14_final$DTF, na.rm=TRUE), 1),
        round(sd(g1y14_final$DTF, na.rm=TRUE), 1),
        round(mean(matfam_GDD_DTF_g1y14$DTF_G1Y14, na.rm=TRUE), 1),
        round(sd(matfam_GDD_DTF_g1y14$DTF_G1Y14, na.rm=TRUE), 1),
        round(mean(patfam_GDD_DTF_g1y14$DTF_G1Y14, na.rm=TRUE), 1),
        round(sd(patfam_GDD_DTF_g1y14$DTF_G1Y14, na.rm=TRUE), 1)
        )
    )

rownames(tabledat) <- c(
    "Pop Mean GDD",
    "Pop SD GDD",
    "Mean Mat GDD",
    "SD Mat GDD",
    "Mean Pat GDD",
    "SD Pat GDD",
    "Pop Mean DTF",
    "Pop SD DTF",
    "Mean Mat DTF",
    "SD Mat DTF",
    "Mean Pat DTF",
    "SD Pat DTF"
    )

print(tabledat)
