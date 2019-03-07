# Nashoba and Kono Phenology - Scripts

- `Quercus_Multigen/nf3mg.p`
    - Type: Pascal

        This is the Pascal source for the `nf3` program in Quercus, originally
        written by F. Shaw. This version has been modified by F. Shaw to
        analyze data from a multi-generational experiment. It will compile
        on Mac OS X with the `fpc` Pascal compiler.
- `Calculate_DTF_GDD.R`
    - Type: R
    - Depends: None

        This script will calculate days to flowering (DTF) and growing
        degree-days to flowering (GDD) for each planting location using the
        PRISM data located in the `Data/` directory.
- `Mean_Fitness_With_SE_Bars.R`
    - Type: R
    - Depends: `aster`

        This script will produce a plot of mean fitness estimates for each
        cohort with standard error bars around the mean.
- `PhenMod_G.R`
    - Type: R
    - Depends: `aster`

        This is the main analysis script. It defines the graphical model and
        fits the Aster model to the data, using the earliest observed germinant
        and earliest observed flower or fruit as predictor variables. It writes
        a R data (RDA) object that is used in other scripts listed here.
- `PhenMod_G_FitScape.R`
    - Type: R
    - Depends: `aster`, `fields`, `graphics`

        This script will produce plots of the fitness landscapes for each
        cohort, using the Aster model fit by `PhenMod_G.R`.
- `PhenMod_G_LA_Beta.R`
    - Type: R
    - Depends: `aster`

        This script will obtain estimates of beta (selection gradients) using
        Lande and Arnold (1983) methodology of ordinary least squares regression
        of relative fitness on `EG` and `EF`. It also uses the additive
        genetic variance-covariance (G matrix) estimated from Quercus to obtain
        predicted changes in mean trait values.
- `PhenMod_G_Residuals.R`
    - Type: R
    - Depends: `aster`

        This script will plot the residuals from the Lande-Arnold regression
        model and the Aster model.
- `SpiderPlot.R`
    - Type: R
    - Depends: `ggplot2`, `reshape`, `gtable`, `grid`

        This script will generate the spider plot to show paternal family means
        of G1Y13 and G1Y14 for `EG` and `EF`. The values of the same paternal
        family will be connected by a line to show genotype-by-environment
        interactions.
