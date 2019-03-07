# Generate the numbers needed for the PhenoFamTable.tex table

# Depends: None
# Input: G1Y13, G1Y14, G2Y14 model-ready data frames
# Output: Numbers for PhenoFamTable.tex: number of maternal families, paternal
#         families, mean dams per sire, etc.

# Note: this file is derived from "MethodsTableRevisions.r"


setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

#   Read in the model-ready files. We eventually want to create a data frame that
#   looks like this:
#       Mat  Initial Seeds_2013
#   keeping in mind that a single Mat may have been part of multiple families
dat_2013 <- read.csv("Data_Files/2013/Model_Ready/2013ModelData.csv", header=TRUE)
dat_2014 <- read.csv("Data_Files/2014/Model_Ready/2014_ModelData_Ver3.csv", header=TRUE)

#############################################################
# There are ~5 bad ids for misc=XZ(15 locs), misc2=ZZ(3 locs), not planted=SS(3 locs), mixed_families=CWBW(1 locs), no_family_info=GK(9 locs)..... length(which(dat_2013$ID == "GK"))
# ! # GK is pedigreed... need to update////  3 Locs not-planted : 18 misc-misc2 locs, 1 misc damily..... 19 locations of non/unresolved-pedigreed seeds 
#?????# [1] 2370 Loc # w/bums (dat_2013), w/o bums (mat3M: sum(dat_2013$Block=="#") B1=167/163, B2=127/122, B3=80/76, B4=55/50, B5=44/40, B6=35/29, B7=28/25
sum(dat_2013$Initial) 	# [1] 2370
sum(dat_2014$Initial)  	# [1] 5193


#############################################################
#2013
###### Portion A ######
#remove rows with Mat=misc, misc2, NA - also works in this case for Pat with those same things
mat1M <- dat_2013[!dat_2013$ID=="SS",]
dim(mat1M) # [1] 533  15:  3 locs (rows) had ID="SS" (total before was 536)
mat2M <- dat_2013[!dat_2013$Mat=="misc",]
dim(mat2M) # 521  15:  12 locs had mat="misc"
mat3M <- mat2M[!mat2M$Mat=="misc2",]
dim(mat3M) # 518  15: 3 locs had mat="misc2"
#get rid of Na mat values in data.frame
mat4M <- mat3M[!is.na(mat3M$Mat),]
dim(mat4M) # 514  15: 4 locs had na in mat some how... 
sum(mat4M$Initial) # 2280

length(unique(mat4M$Mat)) # 111 NOT 110 ... add GK mat fam 32A
length(unique(mat4M$Pat)) # 77 
sum(mat4M$Initial) #2280

#LOC for Each Unique Maternal value ...mat_list still has NA for some reason
mat_list <- list()
for(mat in unique(mat4M$Mat)) {
    #   Find where each Mat occurs in the data
    matloc <- mat4M$Mat == mat
    #   Get the corresponding Locs, as characters
    matloc <- as.character(mat4M$Loc[matloc])
    #   Remove NAs, append to the list
    mat_list[[mat]] <- matloc[!is.na(matloc)] 
} #length 110
mat_list_sum <- data.frame(summary(mat_list))
summary(as.numeric(mat_list_sum$Freq))


pat_list <- list()
for(pat in unique(mat4M$Pat)) {
	# Find where each Pat occurs in the data
	patloc <- mat4M$Pat == pat
	matloc <- mat4M$Mat == mat# Get the corresponding Locs, as characters
	patloc <- as.character(mat4M$Loc[patloc])
	# Remove NAs, append to the list
	pat_list[[pat]] <- patloc[!is.na(patloc)]
} #length 77
pat_list_sum <- data.frame(summary(pat_list))
summary(as.numeric(pat_list_sum$Freq))

#############################################################
###### Portion B ######
# unique Pats per block ; > bk1 <- mat4M[which(mat4M$Block=='1'),]; length(unique(bk1$Pat))
# Bk1=76, Bk2=64, Bk3=48, Bk4=35, Bk5=25, Bk6=20, Bk7=15

# unique Mats per block ; > bk7 <- mat4M[which(mat4M$Block=='7'),]; length(unique(bk7$Mat))
# NEW	update with GK pedigree: Bk1=107, Bk2=93, Bk3=63, Bk4=42, Bk5=33, Bk6=24, Bk7=19

# unique IDs (maternal FS) per block : # bk# <- mat4M[which(mat4M$Block=='#'),]; length(unique(bk#$ID))
# Bk1=154; Bk2=118; Bk3=71; Bk4=48; Bk5=34; Bk6= 24; Bk7=19

# Bk1= 154 of 163 kept Locs; Bk2= 118 of 122 kept Locs; Bk3= 71 of 76 kept Locs; Bk4=47 of 40 kept Locs; Bk5=33 of 40 kept Locs ; Bk6= 23 of 29 kept Locs; Bk7=18 of 25 kept Locs

#? Number of Loc with bad/no pedigree.. 
#easy for Nes - each Mat is its own FS family # separate out "fac"
#2013_dat.....NO FAM: blk1: 2 XZ,1 ZZ =3, bk2: 2+1=3, bk3: 2+1=3, bk4:2+0+CWBW=3, bk4:2+0+0, bk5: 3, bk6: 3, bk7: 2

NesPheno <- mat4M[(mat4M$CrossType %in% "Nes"),] # The Nes portion of this data is 326 records long
length(unique(NesPheno$Mat))  # 76 #aka unique FS families due to pat-nested-HS design
length(unique(NesPheno$Pat)) # 34
length(unique(NesPheno$ID))  # 76

FacPheno <- mat4M[(mat4M$CrossType %in% "Fac"),] # The Fac portion of this data is 179 records long
# number of unique Mat values		
# note: in this group ID with same number but different letters are at least Maternal Half-Sibs
# could code in Quercus - 140A, 140B, 140C with mother ID of 140
length(unique(FacPheno$Mat)) # 35 ...new info...GK pedigreed... old: 34
length(unique(FacPheno$Pat)) # 43
length(unique(FacPheno$ID))  # 82 bum ID removed... in mat4M 
# total IDs Fac plus nes = 76+82=152 .... below calculations = 144 FS families (with same indiv as Mat and same indiv as Pat)

FacPheno$PatMat <- paste(FacPheno$Pat, FacPheno$Mat)
FacPheno$MatPat <- paste(FacPheno$Mat, FacPheno$Pat)

Uni_FacPheno_PatMat <- unique(FacPheno$PatMat)
##################################################################
### Written by TK# Make a character vector that gives paternal and maternal IDs. They are
###### Portion C ######
#####################################AE NEED TO ADD results NUMBERS TO THIS
# separated by a single space
#parental_ids <- c("129B 330A", "330A 129B", "999A 888B", "888B 666D", "666D 888B", "666D 888B")

#Feed IN: Uni_FacPheno_PatMat <- unique(FacPheno$PatMat)
# Get all unique parental IDs.

# Get all unique parental IDs.
uniq_Uni_FacPheno_PatMat <- unique(unlist(strsplit(Uni_FacPheno_PatMat, " ")))

# Then, generate all pairwise combinations of unique parental IDs
pairwise_parental_comb <- combn(uniq_Uni_FacPheno_PatMat, 2)

# pairwise_parental_comb is a matrix - with columns as parental ID pairs. We
# want to apply() over the columns and compare those pairs to the ones we
# observe. This will give us a numeric vector of how many times each possble
# pairwise combination of parental IDs was observed.
parent_combination_counts <- apply(pairwise_parental_comb, 2, function(x) {
    # Generate vectors that give the two directions the cross could occur
    dir1 <- paste(x[1], x[2], sep=" ")
    dir2 <- paste(x[2], x[1], sep=" ")
    # Count how many were crossed in each direction
    dir1_crosses <- sum(Uni_FacPheno_PatMat == dir1)
    dir2_crosses <- sum(Uni_FacPheno_PatMat == dir2)
    # dir1_crosses + dir2_crosses is the total number of families that share
    # these two parents.
    total_shared_parents <- dir1_crosses + dir2_crosses
    return(c(dir1_crosses, dir2_crosses, total_shared_parents))
    })


# We will associate the parent ID combinations with the number of times that
# combination was observed. Do this in a data frame, with Parent_1 being
# the first row of pairwise_parental_comb and Parent_2 being the second row of
# pairwise_parental_comb.
# The number of crosses that is P1 x P2 is the first row of
# parent_combination_counts, the number that is P2 x P1 is the second row,
# and the total number is the third row.
parent_combinations <- data.frame(
    Parent_1=as.character(pairwise_parental_comb[1,]),
    Parent_2=as.character(pairwise_parental_comb[2,]),
    P1xP2=parent_combination_counts[1,],
    P2xP1=parent_combination_counts[2,],
    Num_Fams=parent_combination_counts[3,]
    )

# Then remove any combinations that are observed 0 times.
parent_combinations <- parent_combinations[parent_combinations$Num_Fams > 0,]
print(parent_combinations)


# AE maybe issue.... get 68 for fac, 76 for nest... =144 - but 
#############################################################
#############################################################
###### Portion D ######
########   MEAN DAMS PER SIRE
uniquePat <- unique(mat4M$Pat)
MatNumF <- rep(0,length(uniquePat))
for (i in 1:length(uniquePat))
	{
		placeholder <- subset(mat4M, Pat==uniquePat[i])
		droplevels(placeholder)
		MatNumF[i] <- length(unique(placeholder$Mat))
	}
summary(MatNumF)
#   Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
#   1.000   1.000    2.000    2.052    3.000      5.000 

MatNumF
# [1] 5 3 2 4 4 2 3 3 3 5 1 2 2 3 1 4 2 3 5 1 1 3 2 1 4 1 2 2 5 1 1 1 3 1 3 1 1 1
#[39] 2 3 2 2 2 1 3 3 1 2 1 1 4 1 3 1 2 1 5 2 1 1 2 1 3 2 1 1 1 2 1 1 1 1 1 1 1 1
#[77] 1  #### New ... includes GK pedigree
#############################################
########   MEAN INITIAL SEEDS PER SIRE (PER PAT)

star <- sapply(uniquePat, function(x) {
		uniquePatStar <- mat4M[mat4M$Pat==x,]
	return(sum(uniquePatStar$Initial, na.rm=TRUE))
		})
InitialPat <- data.frame(Pat=uniquePat, star=as.numeric(star))
row.names(InitialPat) <- NULL
# 
summary(InitialPat$star)
# Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
# 2.00    9.00     22.00     29.61   41.00      118.00 ### NEW with GK pedigree inclusion
length(InitialPat$star) # 77
PatMeanInitial <- mean(InitialPat$star) 						#[1] 29.61039  NEW
sdPatInitial <- sd(InitialPat$star) 							#[1] 27.31597  NEW
sePatInitial <- ((sdPatInitial)/sqrt(length(InitialPat$star))) 	#[1] 3.112943  NEW
#############################################
########   MEAN INITIAL SEEDS PER DAM (PER MAT)
uniqueMat <- unique(mat4M$Mat)

star <- sapply(uniqueMat, function(x) {
		uniqueMatStar <- mat4M[mat4M$Mat==x,]
	return(sum(uniqueMatStar$Initial, na.rm=TRUE))
		})
InitialMat <- data.frame(Mat=uniqueMat, star=as.numeric(star))
row.names(InitialMat) <- NULL
# 
summary(InitialMat$star)
#   Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
#   1.00    7.50     15.00     20.54   28.00      77.00 		NEW
length(InitialMat$star) 										# [1] 111 	 NEW
#############################################
########  MEAN SEEDS PLANTED PER SIRE
MatMeanInitial <- mean(InitialMat$star) 						#[1] 20.54054 NEW  
sdMatInitial <- sd(InitialMat$star) 							#[1] 16.47741 NEW
seMatInitial <- ((sdMatInitial)/sqrt(length(InitialMat$star))) 	#[1] 1.563967 NEW
#.
##########################################################################################
##########################################################################################
##########################################################################################
#		NOW FOR 2014 #
rm(list = ls())

dat_2014 <- read.csv("Data_Files/2014/Model_Ready/2014_ModelData_Ver3.csv", header=TRUE)

sum(dat_2014$Initial)  # [1] 5193

# No IDs identifying maternal full-sibling groups
# Appears to be no flubs, misc or NAs... there are some "misc" values
########################################

###### Portion i ###### clean up any bum bits
mat_2014a <- dat_2014[!dat_2014$Mat=="misc",]
mat_2014b <- mat_2014a[!mat_2014a$Mat=="misc2",]
#get rid of Na mat values in data.frame
mat_2014c <- mat_2014b[!is.na(mat_2014b$Mat),]
length(unique(mat_2014c$Mat)) # 121
length(unique(mat_2014c$Pat)) # 138
#don't need to do for pat because the anomalies are across mat and pat groups... at least for this generation
########################################
###### Portion A ######
mat_list <- list()
for(mat in unique(mat_2014c$Mat)) {
    #   Find where each Mat occurs in the data
    matloc <- mat_2014c$Mat == mat
    #   Get the corresponding Locs, as characters
    matloc <- as.character(mat_2014c$Loc[matloc])
    #   Remove NAs, append to the list
    mat_list[[mat]] <- matloc[!is.na(matloc)] 
} #length 
mat_list_sum <- data.frame(summary(mat_list)); summary(as.numeric(mat_list_sum$Freq))
#    Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
#    1.0    10.0     10.0      16.3    30.0       30.
# Maternal families are represented an average of 16.3 times 

pat_list <- list()
for(pat in unique(mat_2014c$Pat)) {
	# Find where each Pat occurs in the data
	patloc <- mat_2014c$Pat == pat
	matloc <- mat_2014c$Mat == mat# Get the corresponding Locs, as characters
	patloc <- as.character(mat_2014c$Loc[patloc])
	# Remove NAs, append to the list
	pat_list[[pat]] <- patloc[!is.na(patloc)]
} #length 77
pat_list_sum <- data.frame(summary(pat_list)); summary(as.numeric(pat_list_sum$Freq))
#   Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
#   1.00   10.00    10.00     14.29   25.00      25.00 
# Paternal families are represented an average of 14.29 times
 #########################################################################################
###### Portion B ######
#    unique Pats per block ; > 
#bk# <- mat_2014c[which(mat_2014c$Block=='#'),]; length(unique(bk#$Pat))
# Bk1=130, Bk2=104, Bk3=87, Bk4=77, Bk5=53

# unique Mats per block ; > bk5 <- mat_2014c[which(mat_2014c$Block=='5'),]; length(unique(bk5$Mat))
# Bk1=118, Bk2=87, Bk3=65, Bk4=56, Bk5=45

# unique IDs (maternal FS) per block : #Need to do this another way for this group ... but should equal the
# number of unique maternal IDs should be the number of full-sib families.. but doesn't include account for . ... 

# number of unique Mat values		
# note: in this group ID with same number but different letters are at least Maternal Half-Sibs
# could code in Quercus - 140A, 140B, 140C with mother ID of 140

# Number of unique Maternal and Paternal IDs 
length(unique(mat_2014c$Mat)) # 121 
length(unique(mat_2014c$Pat)) # 138

#Don't need to do renaming etc, this batch had the same grossing design and setting
mat_2014c$PatMat <- paste(mat_2014c$Pat, mat_2014c$Mat)
#mat_2014c$MatPat <- paste(mat_2014c$Mat, mat_2014c$Pat) # probs don't need this

Uni_mat_2014c_PatMat <- unique(mat_2014c$PatMat)

### Check values estimated here with those in Ch3_Table1_Methods.tex

##################################################################
### Written by TK# Make a character vector that gives paternal and maternal IDs. They are
###### Portion C ######

# separated by a single space
#parental_ids <- c("129B 330A", "330A 129B", "999A 888B", "888B 666D", "666D 888B", "666D 888B")

#Feed IN: Uni_FacPheno_PatMat <- unique(FacPheno$PatMat)
# Get all unique parental IDs.

# Get all unique parental IDs.
uniq_Uni_mat_2014c_PatMat <- unique(unlist(strsplit(Uni_mat_2014c_PatMat, " ")))

# Then, generate all pairwise combinations of unique parental IDs
pairwise_parental_comb <- combn(uniq_Uni_mat_2014c_PatMat, 2)

# pairwise_parental_comb is a matrix - with columns as parental ID pairs. We
# want to apply() over the columns and compare those pairs to the ones we
# observe. This will give us a numeric vector of how many times each possble
# pairwise combination of parental IDs was observed.
parent_combination_counts <- apply(pairwise_parental_comb, 2, function(x) {
    # Generate vectors that give the two directions the cross could occur
    dir1 <- paste(x[1], x[2], sep=" ")
    dir2 <- paste(x[2], x[1], sep=" ")
    # Count how many were crossed in each direction
    dir1_crosses <- sum(Uni_mat_2014c_PatMat== dir1)
    dir2_crosses <- sum(Uni_mat_2014c_PatMat== dir2)
    # dir1_crosses + dir2_crosses is the total number of families that share
    # these two parents.
    total_shared_parents <- dir1_crosses + dir2_crosses
    return(c(dir1_crosses, dir2_crosses, total_shared_parents))
    })
summary(mat_2014c$Initial)
#   Min.    1st Qu.  Median    Mean   3rd Qu.    Max. 
#   1.000   5.000    5.000    4.922   5.000      5.000 

# We will associate the parent ID combinations with the number of times that
# combination was observed. Do this in a data frame, with Parent_1 being
# the first row of pairwise_parental_comb and Parent_2 being the second row of
# pairwise_parental_comb.
# The number of crosses that is P1 x P2 is the first row of
# parent_combination_counts, the number that is P2 x P1 is the second row,
# and the total number is the third row.
parent_combinations <- data.frame(
    Parent_1=as.character(pairwise_parental_comb[1,]),
    Parent_2=as.character(pairwise_parental_comb[2,]),
    P1xP2=parent_combination_counts[1,],
    P2xP1=parent_combination_counts[2,],
    Num_Fams=parent_combination_counts[3,]
    )

# Then remove any combinations that are observed 0 times.
parent_combinations <- parent_combinations[parent_combinations$Num_Fams > 0,]
print(parent_combinations)

#number of rows in parent_combinations = 281
sum(parent_combinations$Num_Fams) # 424
424-281 # 143
#.
#.
#.
#.FS seeds per loc, mean, median, min
summary(mat_2014c$Initial)
#   Min.   1st Qu.  Median    Mean   3rd Qu.    Max. 
#  1.000   5.000    5.000     4.922   5.000     5.000 
#############################################################
#############################################################
###### Portion D ######
########   MEAN DAMS PER SIRE
uniquePat_14 <- unique(mat_2014c$Pat)
 
MatNumF_14 <- rep(0,length(uniquePat_14))
for (i in 1:length(uniquePat_14))
	{
		placeholder <- subset(mat_2014c, Pat==uniquePat_14[i])
		droplevels(placeholder)
		MatNumF_14[i] <- length(unique(placeholder$Mat))
	}
summary(MatNumF_14)
#   Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
#   1.00    2.00     2.00      3.08    4.00       8.00 
MatNumF_14

# [1] 3 2 6 5 3 4 7 7 3 6 5 1 2 4 3 3 3 6 2 5 2 3 6 4 5 5 7 5 1 7 2 2 2 2 7 1 4
# [38] 4 6 6 2 6 5 4 6 5 2 6 5 2 4 1 1 3 6 3 5 7 4 6 2 3 6 2 6 4 2 4 5 2 2 4 2 2
# [75] 3 3 1 1 4 8 2 4 2 1 2 5 2 1 7 1 1 1 2 2 3 3 2 1 1 3 2 4 1 2 2 1 4 1 2 4 2
# [112] 3 1 3 2 2 1 2 1 2 2 4 3 1 1 2 1 2 1 1 1 1 1 1 1 2 1 1
#############################################
########   MEAN INITIAL SEEDS PER SIRE (PER PAT)
star_14 <- sapply(uniquePat_14, function(x) {
		uniquePatStar_14 <- mat_2014c[mat_2014c$Pat==x,]
		return(sum(uniquePatStar_14$Initial, na.rm=TRUE))
		})
InitialPat_14 <- data.frame(Pat=uniquePat_14, star_14=as.numeric(star_14))
row.names(InitialPat_14) <- NULL
# 
summary(InitialPat_14$star_14)
#   Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
#   1.00   10.00    35.00     37.63   50.00      125.00 

length(InitialPat_14$star_14) # 138
PatMeanInitial_14 <- mean(InitialPat_14$star_14) # 37.63043
sdPatInitial_14 <- sd(InitialPat_14$star_14) # 30.77469
sePatInitial_14 <- ((sdPatInitial_14)/sqrt(length(InitialPat_14$star_14))) #[1] 2.619715
#############################################
########   MEAN INITIAL SEEDS PER DAM (PER MAT)
uniqueMat_14 <- unique(mat_2014c$Mat)

star_14m  <- sapply(uniqueMat_14, function(x) {
		uniqueMatStar_14 <- mat_2014c[mat_2014c$Mat==x,]
	return(sum(uniqueMatStar_14$Initial, na.rm=TRUE))
		})
InitialMat_14 <- data.frame(Mat=uniqueMat_14, star_14m=as.numeric(star_14m))
row.names(InitialMat_14) <- NULL
# 
summary(InitialMat_14$star_14m)
#   Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
#   1.00   10.00    30.00     42.92   60.00      170.00 
length(InitialMat_14$star_14m) #  121
MatMeanInitial_14m <- mean(InitialMat_14$star_14m) # 42.91736  
sdMatInitial_14m <- sd(InitialMat_14$star_14m) # 39.50666
seMatInitial_14m <- ((sdMatInitial_14m)/sqrt(length(InitialMat_14$star_14m))) #[1] 3.591515
seMatInitial_14m
