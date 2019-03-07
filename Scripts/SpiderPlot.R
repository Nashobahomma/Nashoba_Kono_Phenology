# Make spider plots of paternal family means for G1Y13 and G1Y14. EG and EF will
# be plotted.

# 2018-09-01
# TK and AN

# Requires: None
# Output: EG and EF spider plots

library(ggplot2)
library(reshape)
library(gtable)
library(grid)

setwd("/Users/tomkono/Dropbox/GitHub/Hoshonti_Aiona")
# setwd("/Users/ambereule-nashoba/Desktop/Dropbox/GitHubRepositories/Hoshonti_Aiona")

#   Read the raw data files. We want to create a data frame that looks like this
#
#   PatL    PatEG   PatEFP  PatSireEG   PatSireEFP
#   19C     131     190     ...         ...
#
# Where PatL is the paternal ID in 2014, and PatSire is the paternal ID in 2013

#   Read in the necessary data files.
#       First read in 2013
pat_2013 <- read.csv("Data_Files/2013/Model_Ready/2013ModelData.csv", header=TRUE)
pat_2014 <- read.csv("Data_Files/2014/Model_ready/2014_ModelData_Ver3.csv", header=TRUE)

#   "Slice" down the data sheets to just the variables of interest
#   From 2013:
#       Loc
#       Pat
#       NumSeeds
#       EarliestGerm
#       EarliestFlowPod
#   From 2014:
#       Loc
#       PatSire
#       NumSeeds
#       EarliestGerm
#       EarliestFlowPod
pat_2013 <- pat_2013[, c("Loc", "Pat", "Mat", "NumSeeds", "EarliestGerm", "EarliestFlowPod")]
pat_2014 <- pat_2014[, c("Loc", "PatSire", "MatDam", "NumSeeds", "EarliestGerm", "EarliestFlowPod")]

#   Then, find out which 2013 PatSire IDs are represented in 2014.
planted_in_both <- intersect(levels(pat_2013$Pat), levels(pat_2014$PatSire))
#   Them, for each 2013 Sire ID, calculate the mean EarliestGerm in 2013 and
#   2014, the mean EarlestFlowPod in 2013 and 2014, and the total of Numseeds
#   in 2013 and 2014.
pateg <- sapply(
    planted_in_both,
    function(x) {
        eg <- pat_2013[pat_2013$Pat == x,]
        return(mean(eg$EarliestGerm, na.rm=T))
        })
patefp <- sapply(
    planted_in_both,
    function(x) {
        efp <- pat_2013[pat_2013$Pat == x,]
        return(mean(efp$EarliestFlowPod, na.rm=T))
        })
patsireeg <- sapply(
    planted_in_both,
    function(x) {
        eg <- pat_2014[pat_2014$Pat == x,]
        return(mean(eg$EarliestGerm, na.rm=T))
        })
patsireefp <- sapply(
    planted_in_both,
    function(x) {
        efp <- pat_2014[pat_2014$Pat == x,]
        return(mean(efp$EarliestFlowPod, na.rm=T))
        })

#   Combine it all into a "wide" data frame. We will convert this into a long
#   format later, for ggplot2. For plotting, remove rows that have NA.
spider_data <- data.frame(
    PatL=planted_in_both,
    EarliestGerm_2013=pateg,
    EarliestFlowPod_2013=patefp,
    EarliestGerm_2014=patsireeg,
    EarliestFlowPod_2014=patsireefp
    )


#   Let's make the plot for EarliestGerm first
#       Remove values that are Inf for EG in 2013 or 2014
spider_data_eg <- spider_data[(is.finite(spider_data$EarliestGerm_2013) & is.finite(spider_data$EarliestGerm_2014)),]
eg_long <- melt(spider_data_eg[, c("PatL", "EarliestGerm_2013", "EarliestGerm_2014")], idvar="PatL")
eg_spider <- ggplot(eg_long, aes(x=variable, y=value, group=PatL)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    theme(
        legend.position="none",
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
        ) +
    labs(x="", y="EG (Julian Date)") +
    scale_x_discrete(name ="",
                    breaks=c("EarliestGerm_2013", "EarliestGerm_2014"),
                    labels=c(expression(paste("G"[1], "Y"[13])), expression(paste("G"[1], "Y"[14])))) +
    ylim(130, 190)

#   Then EarliestFlowPod
#       Remove values that are Inf for EF in 2013 or 2014
spider_data_efp <- spider_data[(is.finite(spider_data$EarliestFlowPod_2013) & is.finite(spider_data$EarliestFlowPod_2014)),]
efp_long <- melt(spider_data_efp[, c("PatL", "EarliestFlowPod_2013", "EarliestFlowPod_2014")], idvar="PatL")
ef_spider <- ggplot(efp_long, aes(x=variable, y=value, group=PatL)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    theme(
        legend.position="none",
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
        ) +
    labs(x="", y="EF (Julian Date)") +
    scale_x_discrete(name ="",
                    breaks=c("EarliestFlowPod_2013", "EarliestFlowPod_2014"),
                    labels=c(expression(paste("G"[1], "Y"[13])), expression(paste("G"[1], "Y"[14]))))

# Combine them into one plot, like the violins
pdf("Results/PhenFigures/G1_PatFams_SpiderPlot.pdf", width=4, height=4)
eg_spider <- ggplotGrob(eg_spider)
ef_spider <- ggplotGrob(ef_spider)
comb_spiders <- rbind(ef_spider, eg_spider, size="first")
# comb_violins$widths <- unit.pmax(eg_violin$widths, ef_violin$widths)
grid.draw(comb_spiders)
dev.off()

pdf("Results/PhenFigures/G1_PatFams_CombinedSpider.pdf", width=4, height=4)
eg_long$Var <- c("EG")
eg_long$PatL <- paste("EG", eg_long$PatL)
eg_long$variable <- as.character(eg_long$variable)
eg_long$variable[eg_long$variable == "EarliestGerm_2013"] <- "G1Y13"
eg_long$variable[eg_long$variable == "EarliestGerm_2014"] <- "G1Y14"
efp_long$Var <- c("EF")
efp_long$PatL <- paste("EF", efp_long$PatL)
efp_long$variable <- as.character(efp_long$variable)
efp_long$variable[efp_long$variable == "EarliestFlowPod_2013"] <- "G1Y13"
efp_long$variable[efp_long$variable == "EarliestFlowPod_2014"] <- "G1Y14"
comb_long <- rbind(eg_long, efp_long)
comb_long$variable <- as.factor(comb_long$variable)
print(comb_long)
spider <- ggplot(
    comb_long, aes(x=variable, y=value, group=PatL, col=Var)) +
    geom_line(size=0.25) +
    geom_point() +
    theme_bw() +
    theme(
        legend.title=element_blank(),
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()
        ) +
    labs(x="", y="Julian Date") +
    scale_x_discrete(name ="",
                    breaks=c("G1Y13", "G1Y14"),
                    labels=c(expression(paste("G"[1], "Y"[13])), expression(paste("G"[1], "Y"[14])))) +
    scale_color_manual(
        values=c("#377eb8", "#4daf4a"),
        breaks=c("EG", "EF"))
spider
dev.off()
