
##### Martes Range-Wide Models ------------------------------------------------

# Hybrid Zone models of Martes americana and M. caurina


##### Packages -------------------------------------------------------------------------------------------------------------

# loading in the requisite packages

packs <- c('ade4', 'adehabitatMA', 'adehabitatHR', 'alphahull', 'data.table', 'dismo', 'dplyr', 'fields', 'ecospat', 
           'ggplot2', 'ggpubr', 'jsonlite', 'kuenm', 'matrixStats', 'raster', 'rgdal', 'rgeos', 'rJava', 'sf', 'sp', 
           'splitstackshape', 'stars', 'terra', 'tidyr', 'utils', 'wesanderson', 'PerformanceAnalytics', 'SDMTools')

# install.packages("remotes")
# remotes::install_version("SDMTools", "1.1-221")

# loading in the packages
lapply(packs, library, character.only=T)

# giving the package versions
packs.df <- as.data.frame(matrix(NA, nrow=length(packs), ncol=2))
colnames(packs.df) <- c('pkg.name', 'pkg.version')
for(i in 1:length(packs)){
      packs.df[i,1] <- packs[i]
      packs.df[i,2] <- as.character(packageVersion(packs[i]))
}
packs.df

# optionally writing the package versions as a csv file
# write.csv(packs.df 'package.versions.csv')

# test if the maxent.jar file is in the dismo folder
list.files( system.file('java', package = 'dismo') )
# [1] "dismo.jar"  "maxent.jar"

sort( sapply(ls(),function(x){object.size(get(x))})) 

# usethis::edit_r_environ()

# Analyses ran on 2023 16" MacBook Pro with M2 Chip with 96GB RAM





##### Data Input -------------------------------------------------------

# working directory
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP")
# sourcing all the custom functions needed to run the analyses
source('./Martes_R/NAF_Functions.R')

### reading in occurrence files

## hybrid zone occ data

# occurrence points for maam inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maam_hyzo_thin <- read.csv('./Martes_masterFiles/occ.maam_hyzo.thin.csv', header=T) %>% as_tibble()

# occurrence points for maam inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maca_hyzo_thin <- read.csv('./Martes_masterFiles/occ.maca_hyzo.thin.csv', header=T) %>% as_tibble()

# occurrence points for martes hybrids (mahy) inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.mahy_hyzo_thin <- read.csv('./Martes_masterFiles/occ.mahy_hyzo.thin.csv', header=T) %>% as_tibble()

# occurrence points for maam inside the hybrid zone
# not thinned, but reduced to one occ per  grid cell per taxa
occ.maam_hyzo_unthin <- read.csv('./Martes_masterFiles/occ.maam_hyzo.unthin.csv', header=T) %>% as_tibble()

# occurrence points for maam inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maca_hyzo_unthin <- read.csv('./Martes_masterFiles/occ.maca_hyzo.unthin.csv', header=T) %>% as_tibble()

# occurrence points for martes hybrids (mahy) inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.mahy_hyzo_unthin <- read.csv('./Martes_masterFiles/occ.mahy_hyzo.unthin.csv', header=T) %>% as_tibble()


### reading in background files

# background of the hybrid zone
# 100km buffer around a maximum convex polygon of all enotped individuals with admixed dna
back.hyzo <- read.csv('./Martes_masterFiles/back.hyzo.csv', header=T) %>% as_tibble()





##### Americana Hybrid Zone Thin Models ------------------------------------

Sys.time()
start <- Sys.time()

### finding the optimal type of models

## M. americana thinned
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maam_Thin")
maam.hyzo.thin.eval <- create.eval.df('Maam.hyzo.thin', 'Modern', 'Hyzo',
                                      f.class = c('LQ', 'LQH', 'H', 'LQT', 'T', 'LQHT'),
                                      beta.values  = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5))

head(maam.hyzo.thin.eval)
create.folders.for.maxent(maam.hyzo.thin.eval)

maam.hyzo.thin.eval <- maxent.crossval.error(eval.df=maam.hyzo.thin.eval,
                                           occs=occ.maam_hyzo_thin,
                                           background = back.hyzo,
                                           predic = 6:9,
                                           first.occ.col = 11,
                                           first.test.col = 16,
                                           omission.rate = 0.1,
                                           all.background = back.hyzo)


maam.hyzo.thin.sum <- maxent.eval.summarize(eval=maam.hyzo.thin.eval)
maam.hyzo.thin.sum
# # A tibble: 48 × 6
# # Groups:   Features [6]
#    Features Betas test.sens AUC_ratio weights standard1mSE
#    <chr>    <dbl>     <dbl>     <dbl>   <dbl>        <dbl>
#  1 H         2        0.833      1.63    1.37        0.301
#  2 LQH       2.5      0.833      1.64    1.37        0.330
#  3 LQHT      2.5      0.833      1.64    1.37        0.352
#  4 LQH       2        0.833      1.64    1.36        0.431
#  5 LQHT      2        0.833      1.63    1.36        0.477
#  6 H         2.5      0.833      1.62    1.35        0.685
#  7 H         0.25     0.733      1.85    1.34        0.871
#  8 LQH       0.25     0.733      1.85    1.34        0.872
#  9 LQ        0.1      0.833      1.73    1.44        1.00 
# 10 LQH       0.5      0.733      1.74    1.29        1.85 
# … with 38 more rows
# # ℹ Use `print(n = ...)` to see more rows

# H_2 seems to have reasonable response curves

write.csv(maam.hyzo.thin.eval$summary, 'maam.hyzo.thin.eval.csv', row.names = F)
write.csv(maam.hyzo.thin.sum, 'maam.hyzo.thin.sum.csv', row.names = F)




##### Americana Hybrid Zone Unthin Models ------------------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maam_Unthin")
maam.hyzo.unthin.eval <- create.eval.df('Maam.hyzo.unthin', 'Modern', 'Hyzo',
                                        f.class = c('LQ', 'LQH', 'H', 'LQT', 'T', 'LQHT'),
                                      beta.values  = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5))

head(maam.hyzo.unthin.eval)
create.folders.for.maxent(maam.hyzo.unthin.eval)

maam.hyzo.unthin.eval <- maxent.crossval.error(eval.df=maam.hyzo.unthin.eval,
                                             occs=occ.maam_hyzo_unthin,
                                             background = back.hyzo,
                                             predic = 6:9,
                                             first.occ.col = 11,
                                             first.test.col = 16,
                                             omission.rate = 0.1,
                                             all.background = back.hyzo)


maam.hyzo.unthin.sum <- maxent.eval.summarize(eval=maam.hyzo.unthin.eval)
maam.hyzo.unthin.sum
# # A tibble: 48 × 6
# # Groups:   Features [6]
#    Features Betas test.sens AUC_ratio weights standard1mSE
#    <chr>    <dbl>     <dbl>     <dbl>   <dbl>        <dbl>
#  1 LQHT      0.75     0.843      1.82    1.53        0.164
#  2 H         0.5      0.843      1.81    1.53        0.166
#  3 LQ        0.75     0.889      1.72    1.53        0.188
#  4 LQH       0.5      0.843      1.81    1.53        0.200
#  5 LQ        0.5      0.889      1.75    1.55        0.324
#  6 H         0.25     0.843      1.79    1.51        0.447
#  7 LQH       0.25     0.843      1.79    1.51        0.464
#  8 LQ        0.1      0.889      1.76    1.57        0.567
#  9 LQHT      2        0.839      1.77    1.50        0.751
# 10 LQ        1        0.889      1.68    1.50        0.767
# # … with 48 more rows
# # ℹ Use `print(n = ...)` to see more rows

# H_0.5 has slightly better-looking response curves that LQHT_0.75, but it's minor
# if they were both evaluated, they would probably be highly correlated with each other

write.csv(maam.hyzo.unthin.eval$summary, 'maam.hyzo.unthin.eval.csv', row.names = F)
write.csv(maam.hyzo.unthin.sum, 'maam.hyzo.unthin.sum.csv', row.names = F)





##### Caurina Hybrid Zone Thin Models ------------------------------------

### finding the optimal type of models

## M. caurina thinned
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maca_Thin")
maca.hyzo.thin.eval <- create.eval.df('Maca.hyzo.thin', 'Modern', 'Hyzo',
                                      f.class = c('LQ', 'LQH', 'H', 'LQT', 'T', 'LQHT'),
                                      beta.values  = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5))

head(maca.hyzo.thin.eval)
create.folders.for.maxent(maca.hyzo.thin.eval)

maca.hyzo.thin.eval <- maxent.crossval.error(eval.df=maca.hyzo.thin.eval,
                                             occs=occ.maca_hyzo_thin,
                                             background = back.hyzo,
                                             predic = 6:9,
                                             first.occ.col = 11,
                                             first.test.col = 16,
                                             omission.rate = 0.1,
                                             all.background = back.hyzo)


maca.hyzo.thin.sum <- maxent.eval.summarize(eval=maca.hyzo.thin.eval)
maca.hyzo.thin.sum
# # A tibble: 48 × 6
# # Groups:   Features [6]
#    Features Betas test.sens AUC_ratio weights standard1mSE
#    <chr>    <dbl>     <dbl>     <dbl>   <dbl>        <dbl>
#  1 H         2.5      0.887      1.36    1.21        0.245
#  2 LQ        0.5      0.813      1.45    1.18        0.245
#  3 H         1        0.853      1.42    1.21        0.284
#  4 LQ        0.25     0.847      1.43    1.22        0.395
#  5 H         2        0.887      1.37    1.22        0.407
#  6 LQH       1        0.813      1.42    1.15        0.746
#  7 LQ        0.75     0.813      1.41    1.15        0.858
#  8 LQ        1        0.813      1.41    1.14        0.890
#  9 H         0.75     0.813      1.40    1.14        0.961
# 10 H         1.5      0.887      1.41    1.25        1.00  
# # … with 38 more rows
# # ℹ Use `print(n = ...)` to see more rows

# H_2.5 and LQ_0.5 both seem to have reasonable response curves
# if they were both evaluated, they would probably be highly correlated with each other
# going with H since it has better testing sensitivity

write.csv(maca.hyzo.thin.eval$summary, 'maca.hyzo.thin.eval.csv', row.names = F)
write.csv(maca.hyzo.thin.sum, 'maca.hyzo.thin.sum.csv', row.names = F)



##### Caurina Hybrid Zone Unthin Models ------------------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maca_Unthin")
maca.hyzo.unthin.eval <- create.eval.df('Maca.hyzo.unthin', 'Modern', 'Hyzo',
                                        f.class = c('LQ', 'LQH', 'H', 'LQT', 'T', 'LQHT'),
                                        beta.values  = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5))

head(maca.hyzo.unthin.eval)
create.folders.for.maxent(maca.hyzo.unthin.eval)

maca.hyzo.unthin.eval <- maxent.crossval.error(eval.df=maca.hyzo.unthin.eval,
                                               occs=occ.maca_hyzo_unthin,
                                               background = back.hyzo,
                                               predic = 6:9,
                                               first.occ.col = 11,
                                               first.test.col = 16,
                                               omission.rate = 0.1,
                                               all.background = back.hyzo)

maca.hyzo.unthin.sum <- maxent.eval.summarize(eval=maca.hyzo.unthin.eval)
maca.hyzo.unthin.sum
# # A tibble: 48 × 6
# # Groups:   Features [6]
#    Features Betas test.sens AUC_ratio weights standard1mSE
#    <chr>    <dbl>     <dbl>     <dbl>   <dbl>        <dbl>
#  1 LQHT       2.5     0.879      1.42    1.25       0.0963
#  2 LQH        2       0.888      1.40    1.25       0.103 
#  3 LQH        2.5     0.888      1.41    1.25       0.247 
#  4 LQH        1       0.871      1.42    1.24       0.263 
#  5 H          2       0.879      1.41    1.24       0.280 
#  6 LQ         1.5     0.888      1.41    1.25       0.284 
#  7 LQT        2.5     0.879      1.42    1.25       0.321 
#  8 LQT        1.5     0.862      1.43    1.24       0.340 
#  9 LQ         1       0.888      1.41    1.25       0.348 
# 10 LQHT       2       0.888      1.39    1.23       0.360 
# # … with 38 more rows
# # ℹ Use `print(n = ...)` to see more rows

# LQHT_2.5 seems to have reasonable response curves

write.csv(maca.hyzo.unthin.eval$summary, 'maca.hyzo.unthin.eval.csv', row.names = F)
write.csv(maca.hyzo.unthin.sum, 'maca.hyzo.unthin.sum.csv', row.names = F)


end <- Sys.time()
end - start
Sys.time
# Time difference of 46.8073 mins





##### Evaluating Americana Thin Models -------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maam_Thin")
# making the evaluation data.frame for maam and creating the requisite folders
maam.hyzo.thin.eval1 <- create.eval.df('Maam', 'modern', 'hyzo.thin', beta.values=2, f.class='H'); maam.hyzo.thin.eval1
create.folders.for.maxent(maam.hyzo.thin.eval1)

# running the evaluation object
maam.hyzo.thin.eval1 <- maxent.crossval.error(eval.df=maam.hyzo.thin.eval1, # evaluation df
                                             occs = occ.maam_hyzo_thin,# occurrences
                                             background = back.hyzo,   # background
                                             predic = 6:9,             # column number of predictor variables
                                             first.occ.col = 11,       # first training column
                                             first.test.col = 16,      # first testing column
                                             omission.rate=0.1,        # omission rate
                                             all.background=back.hyzo
)
# View(maam.hyzo.thin.eval1$summary)
write.csv(maam.hyzo.thin.eval1$summary, 'maam.hyzo.thin.eval_H_2.csv', row.names=F)

### variable contribution & importance
get.var.contrib.importance(eval=maam.hyzo.thin.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10  1.3228143 2.4389234
# bio12 97.8914753 3.3544793
# bio15  0.3101002 0.4325547
# bio4   0.4756086 0.9134818
# 
# $VarImportance
#          w.mean     w.sd
# bio10  1.332059 1.589064
# bio12 94.471642 7.235590
# bio15  1.198966 1.308334
# bio4   2.997358 5.756902

### variable ranges
get.var.ranges(maam.hyzo.thin.eval1)
#          bio4     bio10 bio12    bio15
# min  657.4054  5.976667   191 12.66888
# max 1167.8830 22.607500  1230 74.82084

### response curves
maam.hyzo.thin.rc <- get.maxent.response.curves(eval=maam.hyzo.thin.eval1, expand=0)
maam.hyzo.thin.rc.poly <- process.maxent.response.curves(maam.hyzo.thin.rc)

# weighting the cross-validation reps and projecting to all occ and background points
maam.hyzo.means.thin1 <- maxent.eval(eval=maam.hyzo.thin.eval1)
nrow(maam.hyzo.means.thin1$occ)
# [1] 13
# picking a quantile of the occs to define a threshold
# 10% 
quantile(maam.hyzo.means.thin1$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))
#        0%        1%        5%       10%       20%       25% 
# 0.2359541 0.2841105 0.4767361 0.6414816 0.7040258 0.7724975 

# the lowest 10% occurrence threshold
maam.thin.thresh1 <- 0.6414816

suit.uncert.plot(maam.hyzo.means.thin1)
ggsave('maam.hyzo.thin.means.pdf')


### projecting model to hybrid zone background
maam.hyzo.thin.everything <- maxent.everything(eval=maam.hyzo.thin.eval1,
                                               means=maam.hyzo.means.thin1,
                                               everything=back.hyzo,
                                               predic=6:9,
                                               thresh=maam.thin.thresh1,
                                               name='maam_hyzo.thin')
write.csv(maam.hyzo.thin.everything, 'maam.hyzo.thin.everything.csv', row.names=F)



##### Evaluating Americana Unthin Models -------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maam_Unthin")
# making the evaluation data.frame for maam and creating the requisite folders
maam.hyzo.unthin.eval1 <- create.eval.df('Maam', 'modern', 'hyzo.unthin', beta.values=0.75, f.class='LQHT'); maam.hyzo.unthin.eval1
create.folders.for.maxent(maam.hyzo.unthin.eval1)

# running the evaluation object
maam.hyzo.unthin.eval1 <- maxent.crossval.error(eval.df=maam.hyzo.unthin.eval1, # evaluation df
                                               occs = occ.maam_hyzo_unthin,# occurrences
                                               background = back.hyzo,   # background
                                               predic = 6:9,             # column number of predictor variables
                                               first.occ.col = 11,       # first training column
                                               first.test.col = 16,      # first testing column
                                               omission.rate=0.1,        # omission rate
                                               all.background=back.hyzo
)
# View(maam.hyzo.unthin.eval1$summary)
write.csv(maam.hyzo.unthin.eval1$summary, 'maam.hyzo.unthin.eval_LQHT_0.75.csv', row.names=F)

### variable contribution & importance
get.var.contrib.importance(eval=maam.hyzo.unthin.eval1)
# $VarContribution
#        w.mean      w.sd
# bio10 20.89383 2.966083
# bio12 48.48206 1.130057
# bio15 17.06699 5.848159
# bio4  13.55716 5.302238
# 
# $VarImportance
#          w.mean     w.sd
# bio10 22.930292 3.847195
# bio12 28.733872 2.419967
# bio15 40.090430 7.594261
# bio4   8.245405 5.730571


### variable ranges
get.var.ranges(maam.hyzo.unthin.eval1)
#          bio4     bio10 bio12    bio15
# min  657.4054  5.976667   191 12.66888
# max 1167.8830 22.607500  1230 74.82084

### response curves
maam.hyzo.unthin.rc <- get.maxent.response.curves(eval=maam.hyzo.unthin.eval1, expand=0)
maam.hyzo.unthin.rc.poly <- process.maxent.response.curves(maam.hyzo.unthin.rc)


# weighting the cross-validation reps and projecting to all occ and background points
maam.hyzo.means.unthin1 <- maxent.eval(eval=maam.hyzo.unthin.eval1)
nrow(maam.hyzo.means.unthin1$occ)
# [1] 36
# picking a quantile of the occs to define a threshold
# 10% 
quantile(maam.hyzo.means.unthin1$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))
#         0%         1%         5%        10%        20%        25% 
# 0.02008053 0.07524930 0.22664840 0.34361619 0.47344464 0.58011977

# the lowest 10% occurrence threshold
maam.unthin.thresh1 <- 0.34361619

suit.uncert.plot(maam.hyzo.means.unthin1)
ggsave('maam.hyzo.unthin.means.pdf')


### projecting model to hybrid zone background
maam.hyzo.unthin.everything <- maxent.everything(eval=maam.hyzo.unthin.eval1,
                                                 means=maam.hyzo.means.unthin1,
                                                 everything=back.hyzo,
                                                 predic=6:9,
                                                 thresh=maam.unthin.thresh1,
                                                 name='maam_hyzo.unthin')
write.csv(maam.hyzo.unthin.everything, 'maam.hyzo.unthin.everything.csv', row.names=F)





##### Evaluating Caurina Thin Models -------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maca_Thin")
# making the evaluation data.frame for maca and creating the requisite folders
maca.hyzo.thin.eval1 <- create.eval.df('Maca', 'modern', 'hyzo.thin', beta.values=2.5, f.class='H'); maca.hyzo.thin.eval1
create.folders.for.maxent(maca.hyzo.thin.eval1)

# running the evaluation object
maca.hyzo.thin.eval1 <- maxent.crossval.error(eval.df=maca.hyzo.thin.eval1, # evaluation df
                                             occs = occ.maca_hyzo_thin,# occurrences
                                             background = back.hyzo,   # background
                                             predic = 6:9,             # column number of predictor variables
                                             first.occ.col = 11,       # first training column
                                             first.test.col = 16,      # first testing column
                                             omission.rate=0.1,        # omission rate
                                             all.background=back.hyzo
)
# View(maca.hyzo.thin.eval1$summary)
write.csv(maca.hyzo.thin.eval1$summary, 'maca.hyzo.thin.eval_H_2.5.csv', row.names=F)

### variable contribution & importance
get.var.contrib.importance(eval=maca.hyzo.thin.eval1)
# $VarContribution
#            w.mean       w.sd
# bio10 70.16812557 17.4248501
# bio12  0.05739845  0.1376624
# bio15 16.65823521 22.0118264
# bio4  13.11624077 18.3343318
# 
# $VarImportance
#         w.mean     w.sd
# bio10 67.64131 25.17765
# bio12  0.00000  0.00000
# bio15 13.25295 17.68871
# bio4  19.10575 28.30390

### variable ranges
get.var.ranges(maca.hyzo.thin.eval1)
#          bio4     bio10 bio12    bio15
# min  657.4054  5.976667   191 12.66888
# max 1167.8830 22.607500  1230 74.82084

### response curves
maca.hyzo.thin.rc <- get.maxent.response.curves(eval=maca.hyzo.thin.eval1, expand=0)
maca.hyzo.thin.rc.poly <- process.maxent.response.curves(maca.hyzo.thin.rc)


# weighting the cross-validation reps and projecting to all occ and background points
maca.hyzo.means.thin1 <- maxent.eval(eval=maca.hyzo.thin.eval1)
nrow(maca.hyzo.means.thin1$occ)
# [1] 27
# picking a quantile of the occs to define a threshold
# 10% 
quantile(maca.hyzo.means.thin1$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))
#        0%        1%        5%       10%       20%       25% 
# 0.4053938 0.4147657 0.4659394 0.5240544 0.5788954 0.6264689 

# the lowest 10% occurrence threshold
maca.thin.thresh1 <- 0.5240544

suit.uncert.plot(maca.hyzo.means.thin1)
ggsave('maca.hyzo.thin.means.pdf')


### projecting model to hybrid zone background
maca.hyzo.thin.everything <- maxent.everything(eval=maca.hyzo.thin.eval1,
                                               means=maca.hyzo.means.thin1,
                                               everything=back.hyzo,
                                               predic=6:9,
                                               thresh=maca.thin.thresh1,
                                               name='maca_hyzo.thin')
write.csv(maca.hyzo.thin.everything, 'maca.hyzo.thin.everything.csv', row.names=F)



##### Evaluating Caurina Unthin Models -------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/Hyzo/Maca_Unthin")
# making the evaluation data.frame for maca and creating the requisite folders
maca.hyzo.unthin.eval1 <- create.eval.df('Maca', 'modern', 'hyzo.unthin', beta.values=2.5, f.class='LQHT'); maca.hyzo.unthin.eval1
create.folders.for.maxent(maca.hyzo.unthin.eval1)

# running the evaluation object
maca.hyzo.unthin.eval1 <- maxent.crossval.error(eval.df=maca.hyzo.unthin.eval1, # evaluation df
                                               occs = occ.maca_hyzo_unthin,# occurrences
                                               background = back.hyzo,   # background
                                               predic = 6:9,             # column number of predictor variables
                                               first.occ.col = 11,       # first training column
                                               first.test.col = 16,      # first testing column
                                               omission.rate=0.1,        # omission rate
                                               all.background=back.hyzo
)
# View(maca.hyzo.unthin.eval1$summary)
write.csv(maca.hyzo.unthin.eval1$summary, 'maca.hyzo.unthin.eval_LQHT_2.5.csv', row.names=F)

### variable contribution & importance
get.var.contrib.importance(eval=maca.hyzo.unthin.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10 50.518351 6.915238
# bio12 25.854412 3.442250
# bio15 19.835000 5.642115
# bio4   3.792236 1.495884
# 
# $VarImportance
#          w.mean     w.sd
# bio10 53.322839 4.875316
# bio12 33.284165 3.465755
# bio15 11.108212 3.239661
# bio4   2.284765 3.597911

### variable ranges
get.var.ranges(maca.hyzo.unthin.eval1)
#          bio4     bio10 bio12    bio15
# min  657.4054  5.976667   191 12.66888
# max 1167.8830 22.607500  1230 74.82084

### response curves
maca.hyzo.unthin.rc <- get.maxent.response.curves(eval=maca.hyzo.unthin.eval1, expand=0)
maca.hyzo.unthin.rc.poly <- process.maxent.response.curves(maca.hyzo.unthin.rc)

# weighting the cross-validation reps and projecting to all occ and background points
maca.hyzo.means.unthin1 <- maxent.eval(eval=maca.hyzo.unthin.eval1)
nrow(maca.hyzo.means.unthin1$occ)
# [1] 116
# picking a quantile of the occs to define a threshold
# 10% 
quantile(maca.hyzo.means.unthin1$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))
#         0%         1%         5%        10%        20%        25% 
# 0.06356826 0.09854981 0.23221264 0.30067011 0.39441418 0.44397863

# the lowest 10% occurrence threshold
maca.unthin.thresh1 <- 0.30067011

suit.uncert.plot(maca.hyzo.means.unthin1)
ggsave('maca.hyzo.unthin.means.pdf')


### projecting model to hybrid zone background
maca.hyzo.unthin.everything <- maxent.everything(eval=maca.hyzo.unthin.eval1,
                                                 means=maca.hyzo.means.unthin1,
                                                 everything=back.hyzo,
                                                 predic=6:9,
                                                 thresh=maca.unthin.thresh1,
                                                 name='maca_hyzo.unthin')
write.csv(maca.hyzo.unthin.everything, 'maca.hyzo.unthin.everything.csv', row.names=F)








##### HyZo Response Curve Plots - Thinned --------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")

hyzo.ranges <- get.var.ranges(maca.hyzo.unthin.eval1)
get.var.ranges(maca.hyzo.unthin.eval1)


### variable contribution & importance
get.var.contrib.importance(eval=maam.hyzo.thin.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10  1.3228143 2.4389234
# bio12 97.8914753 3.3544793
# bio15  0.3101002 0.4325547
# bio4   0.4756086 0.9134818
# 
# $VarImportance
#          w.mean     w.sd
# bio10  1.332059 1.589064
# bio12 94.471642 7.235590
# bio15  1.198966 1.308334
# bio4   2.997358 5.756902

get.var.contrib.importance(eval=maca.hyzo.unthin.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10 50.518351 6.915238
# bio12 25.854412 3.442250
# bio15 19.835000 5.642115
# bio4   3.792236 1.495884
# 
# $VarImportance
#          w.mean     w.sd
# bio10 53.322839 4.875316
# bio12 33.284165 3.465755
# bio15 11.108212 3.239661
# bio4   2.284765 3.597911


# bio4
bio4 <- ggplot() +
      geom_polygon(data=maca.hyzo.thin.rc.poly$bio4, aes(x=x/100, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.thin.rc.poly$bio4, aes(x=x/100, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.thin.rc$bio4$input/100, y=maca.hyzo.thin.rc$bio4$bio4_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.thin.rc$bio4$input/100, y=maam.hyzo.thin.rc$bio4$bio4_mean), col='red', size=1) + # caurina mean
      geom_vline(xintercept = hyzo.ranges[1,1]/100, linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,1]/100, linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + xlim(0,18) + ylim(0,1) + 
      xlab('Temp. Seasonality (°C)') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=9, y=.9, label=paste0('3.8 ± 1.4%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=9, y=.25, label=paste0('13.6 ± 5.3%'), color='red', size=3.5 ) # americana variable contribution
bio4

# bio10
bio10 <- ggplot() +
      geom_polygon(data=maca.hyzo.thin.rc.poly$bio10, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.thin.rc.poly$bio10, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.thin.rc$bio10$input, y=maca.hyzo.thin.rc$bio10$bio10_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.thin.rc$bio10$input, y=maam.hyzo.thin.rc$bio10$bio10_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = hyzo.ranges[1,2], linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,2], linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + xlim(-2.5, 37.5) + ylim(0,1) + 
      xlab('Mean Summer Temp. (°C)') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=15, y=.8, label=paste0('50.5 ± 6.9%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=12, y=.3, label=paste0('20.9 ± 3.0%'), color='red', size=3.5 ) # americana variable contribution
bio10

# bio12
bio12 <- ggplot() +
      geom_polygon(data=maca.hyzo.thin.rc.poly$bio12, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.thin.rc.poly$bio12, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.thin.rc$bio12$input, y=maca.hyzo.thin.rc$bio12$bio12_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.thin.rc$bio12$input, y=maam.hyzo.thin.rc$bio12$bio12_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = hyzo.ranges[1,3], linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,3], linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + xlim(0,5400) + ylim(0,1) + 
      xlab(expression(paste("Annual Precipitation mm*", yr^{-1}))) + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=2000, y=.45, label=paste0('25.8 ± 3.4%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=2000, y=.85, label=paste0('48.4 ± 1.1%'), color='red', size=3.5 ) # americana variable contribution
bio12

# bio15
bio15 <- ggplot() +
      geom_polygon(data=maca.hyzo.thin.rc.poly$bio15, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.thin.rc.poly$bio15, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.thin.rc$bio15$input, y=maca.hyzo.thin.rc$bio15$bio15_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.thin.rc$bio15$input, y=maam.hyzo.thin.rc$bio15$bio15_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = hyzo.ranges[1,4], linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,4], linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + ylim(0,1) + 
      xlab('Precip. Seasonality') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=65, y=.65, label=paste0('16.6 ± 22.0%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=25, y=.3, label=paste0('0.3 ± 0.4%'), color='red', size=3.5 ) # americana variable contribution
bio15


ggpubr::ggarrange(bio4, bio10, bio12, bio15, 
                  nrow=2, ncol=2, labels = c('A', 'B', 'C', 'D'), label.x=0.8, label.y=0.95, align='hv',
                  font.label=list(size=24, color="black", face="bold", family=NULL) ) + bgcolor('White')
ggsave('MartesResponseCurves_HyzoThin.png', width=20, height=20, units='cm')
ggsave('FigSI8.eps', width=20, height=20, units='cm', dpi=1200)

##### HyZo Response Curve Plots - Unthinned --------------------------


### variable contribution & importance
get.var.contrib.importance(eval=maam.hyzo.unthin.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10 20.89246 2.966742
# bio12 48.48217 1.129944
# bio15 17.07067 5.851185
# bio4  13.55474 5.303256
# 
# $VarImportance
#          w.mean     w.sd
# bio10 22.928227 3.847172
# bio12 28.732862 2.419766
# bio15 40.095255 7.595564
# bio4   8.243655 5.731005

get.var.contrib.importance(eval=maca.hyzo.unthin.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10 50.518351 6.915238
# bio12 25.854412 3.442250
# bio15 19.835000 5.642115
# bio4   3.792236 1.495884
# 
# $VarImportance
#          w.mean     w.sd
# bio10 53.322839 4.875316
# bio12 33.284165 3.465755
# bio15 11.108212 3.239661
# bio4   2.284765 3.597911

# bio4
bio4 <- ggplot() +
      geom_polygon(data=maca.hyzo.unthin.rc.poly$bio4, aes(x=x/100, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.unthin.rc.poly$bio4, aes(x=x/100, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.unthin.rc$bio4$input/100, y=maca.hyzo.unthin.rc$bio4$bio4_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.unthin.rc$bio4$input/100, y=maam.hyzo.unthin.rc$bio4$bio4_mean), col='red', size=1) + # caurina mean
      geom_vline(xintercept = hyzo.ranges[1,1]/100, linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,1]/100, linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + xlim(0,18) + ylim(0,1) + 
      xlab('Temp. Seasonality (°C)') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=9, y=.75, label=paste0('13.0 ± 18.3%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=4, y=.25, label=paste0('0.5 ± 0.9%'), color='red', size=3.5 ) # americana variable contribution
bio4

# bio10
bio10 <- ggplot() +
      geom_polygon(data=maca.hyzo.unthin.rc.poly$bio10, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.unthin.rc.poly$bio10, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.unthin.rc$bio10$input, y=maca.hyzo.unthin.rc$bio10$bio10_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.unthin.rc$bio10$input, y=maam.hyzo.unthin.rc$bio10$bio10_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = hyzo.ranges[1,2], linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,2], linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + xlim(-2.5, 37.5) + ylim(0,1) + 
      xlab('Mean Summer Temp. (°C)') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=15, y=.8, label=paste0('70.2 ± 17.4%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=12, y=.3, label=paste0('1.3 ± 2.4%'), color='red', size=3.5 ) # americana variable contribution
bio10

# bio12
bio12 <- ggplot() +
      geom_polygon(data=maca.hyzo.unthin.rc.poly$bio12, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.unthin.rc.poly$bio12, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.unthin.rc$bio12$input, y=maca.hyzo.unthin.rc$bio12$bio12_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.unthin.rc$bio12$input, y=maam.hyzo.unthin.rc$bio12$bio12_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = hyzo.ranges[1,3], linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,3], linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + xlim(0,5400) + ylim(0,1) + 
      xlab(expression(paste("Annual Precipitation mm*", yr^{-1}))) + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=2000, y=.45, label=paste0('0.1 ± 0.1%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=2000, y=.85, label=paste0('97.9 ± 3.3%'), color='red', size=3.5 ) # americana variable contribution
bio12

# bio15
bio15 <- ggplot() +
      geom_polygon(data=maca.hyzo.unthin.rc.poly$bio15, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.hyzo.unthin.rc.poly$bio15, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.hyzo.unthin.rc$bio15$input, y=maca.hyzo.unthin.rc$bio15$bio15_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.hyzo.unthin.rc$bio15$input, y=maam.hyzo.unthin.rc$bio15$bio15_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = hyzo.ranges[1,4], linetype="dashed", color = "black", size=.5) + # hybrid zone lower bound
      geom_vline(xintercept = hyzo.ranges[2,4], linetype="dashed", color = "black", size=.5) + # hybrid zone upper bound
      theme_classic() + ylim(0,1) + 
      xlab('Precip. Seasonality') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=60, y=.35, label=paste0('19.8 ± 5.6%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=30, y=.75, label=paste0('17.1 ± 5.9%'), color='red', size=3.5 ) # americana variable contribution
bio15


ggpubr::ggarrange(bio4, bio10, bio12, bio15, 
                  nrow=2, ncol=2, labels = c('A', 'B', 'C', 'D'), label.x=0.8, label.y=0.95, align='hv',
                  font.label=list(size=24, color="black", face="bold", family=NULL) ) + bgcolor('White')
ggsave('MartesResponseCurves_HyzoUnthin.png', width=20, height=20, units='cm')
ggsave('FigSI7.eps', width=20, height=20, units='cm', dpi=1200)

##### Map Figures -----------------------------------


# loading in raster data just for re-projecting purposes
setwd("~/Dropbox/Marten_ENMs/Raster Layers/wc2.0_5m_bio")
raster.names <- list.files( pattern='\\.tif$')

raster.names <- raster.names[c(4,10,12,15)]

bc <- terra::rast(paste0('./', raster.names))
bc <- crop(bc, y=extent(-170, -50, 30, 77))


n.am.rast <- raster(xmn=-4000000, xmx=4000000, ymn=-3000000, ymx=3000000, res=10000, crs=n.am.crs)
n.am.rast[n.am.rast] <- rnorm(ncell(n.am.rast))

### projecting the thinned hybrid zone-based models to the entire north american extent
## will be reprojected to focus on the hybrid zone
# projecting model to hybrid zone background
maam.hyzo.thin.n.am <- maxent.everything(eval=maam.hyzo.thin.eval1,
                                         means=maam.hyzo.means1,
                                         everything=back.all,
                                         predic=6:9,
                                         thresh=maam.thin.thresh1,
                                         name='maam_hyzo.thin')
# projecting model to hybrid zone background
maca.hyzo.thin.n.am <- maxent.everything(eval=maca.hyzo.unthin.eval1,
                                         means=maca.hyzo.means1,
                                         everything=back.all,
                                         predic=6:9,
                                         thresh=maca.unthin.thresh1,
                                         name='maca_hyzo.unthin')

### projecting the unthinned hybrid zone-basedmodels to the entire north american extent
# americana
maam.hyzo.unthin.n.am <- maxent.everything(eval=maam.hyzo.unthin.eval1,
                                           means=maam.hyzo.means1,
                                           everything=back.all,
                                           predic=6:9,
                                           thresh=maam.unthin.thresh1,
                                           name='maam_hyzo.unthin')

# caurina
maca.hyzo.unthin.n.am <- maxent.everything(eval=maca.hyzo.unthin.eval1,
                                           means=maca.hyzo.means1,
                                           everything=back.all,
                                           predic=6:9,
                                           thresh=maca.unthin.thresh1,
                                           name='maca_hyzo.unthin')



hyzo.models <- do.call('cbind', list(
      back.all[,1:9],
      maam.hyzo.thin.n.am,
      maca.hyzo.thin.n.am,
      maam.hyzo.unthin.n.am,
      maca.hyzo.unthin.n.am
      )
      )
colnames(hyzo.models) <- c(colnames(back.all[,1:9]),
                           c('Maam_contin.t',
                             'Maam_thresh.t',
                             'Maca_contin.t',
                             'Maca_thresh.t',
                             'Maam_contin.u',
                             'Maam_thresh.u',
                             'Maca_contin.u',
                             'Maca_thresh.u'
                           )
                           )
nrow(hyzo.models)
# [1] 317406

# cropping for speed performancesfor making projections to the hybrid zone
hyzo.models <- hyzo.models %>% filter(longitude > -120 & longitude < -105 & latitude > 40 & latitude < 53)
nrow(hyzo.models)
# [1] 27295
hyzo.models <- as_tibble(hyzo.models)

# americana thinned
hm.a.t <- hyzo.models %>% dplyr::filter(Maam_thresh.t > 0) %>% dplyr::select(longitude, latitude, Maam_contin.t)
hm.a.t <- vect(hm.a.t, crs=wgs84, geom=c('longitude', 'latitude'))
hm.a.t <- rasterize(hm.a.t, bc[[1]], field='Maam_contin.t')
hm.a.t <- terra::project(hm.a.t, n.am.rast)
hm.a.t[hm.a.t > 1] <- 1; hm.a.t[hm.a.t < maam.thin.thresh1] <- 0; plot(hm.a.t)
# caurina thinned
hm.c.t <- hyzo.models %>% dplyr::filter(Maca_thresh.t > 0) %>% dplyr::select(longitude, latitude, Maca_contin.t)
hm.c.t <- vect(hm.c.t, crs=wgs84, geom=c('longitude', 'latitude'))
hm.c.t <- rasterize(hm.c.t, bc[[1]], field='Maca_contin.t')
hm.c.t <- terra::project(hm.c.t, n.am.rast)
hm.c.t[hm.c.t > 1] <- 1; hm.c.t[hm.c.t < maca.thin.thresh1] <- 0; plot(hm.c.t)
# americana unthinned
hm.a.u <- hyzo.models %>% dplyr::filter(Maam_thresh.u > 0) %>% dplyr::select(longitude, latitude, Maam_contin.u)
hm.a.u <- vect(hm.a.u, crs=wgs84, geom=c('longitude', 'latitude'))
hm.a.u <- rasterize(hm.a.u, bc[[1]], field='Maam_contin.u')
hm.a.u <- terra::project(hm.a.u, n.am.rast)
hm.a.u[hm.a.u > 1] <- 1; hm.a.u[hm.a.u < maam.unthin.thresh1] <- 0; plot(hm.a.u)
# caurina unthinned
hm.c.u <- hyzo.models %>% dplyr::filter(Maca_thresh.u > 0) %>% dplyr::select(longitude, latitude, Maca_contin.u)
hm.c.u <- vect(hm.c.u, crs=wgs84, geom=c('longitude', 'latitude'))
hm.c.u <- rasterize(hm.c.u, bc[[1]], field='Maca_contin.u')
hm.c.u <- terra::project(hm.c.u, n.am.rast)
hm.c.u[hm.c.u > 1] <- 1; hm.c.u[hm.c.u < maca.unthin.thresh1] <- 0; plot(hm.c.u)



suit.leg.pos <- c(0.75, 0.9)

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")

# americana suitability map
maam.hyzothin.model.hyzo <- 
      ggplot() + 
      geom_spatraster(data = hm.a.t) +
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'GnBu'), guide="colorbar", 
                           na.value="white",  name='', limits=0:1) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      # geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = suit.leg.pos, 
            axis.text = element_blank(),
            legend.box = "vertical",
            legend.justification = "top",
            legend.margin = margin(0.25, 0.9, 0.25, 0.25, unit = "cm"), # top, right, bottom, left
            legend.box.background = element_rect(fill = "white", colour = 'black'),
            legend.key.width = unit(1.2, 'cm' ), 
            legend.spacing.y = unit(0.2, "lines"),
            legend.text = element_text(size=16), 
            axis.ticks = element_blank(),
            legend.title = element_text(size = 20, colour = "black", face = "bold"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA)
      ) + 
      labs(x = NULL, y = NULL)
maam.hyzothin.model.hyzo
ggsave('maam.hyzothin.model.hyzo.png', height=20, width=19, unit='cm', dpi=320)



maca.hyzothin.model.hyzo <- 
      ggplot() + 
      geom_spatraster(data = hm.c.t) +
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'GnBu'), guide="colorbar", 
                           na.value="white",  name='', limits=0:1) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      # geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = suit.leg.pos, 
            axis.text = element_blank(),
            legend.box = "vertical",
            legend.justification = "top",
            legend.margin = margin(0.25, 0.9, 0.25, 0.25, unit = "cm"), # top, right, bottom, left
            legend.box.background = element_rect(fill = "white", colour = 'black'),
            legend.key.width = unit(1.2, 'cm' ), 
            legend.spacing.y = unit(0.2, "lines"),
            legend.text = element_text(size=16), 
            axis.ticks = element_blank(),
            legend.title = element_text(size = 20, colour = "black", face = "bold"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA)
      ) + 
      labs(x = NULL, y = NULL)
maca.hyzothin.model.hyzo
ggsave('maca.hyzothin.model.hyzo.png', height=20, width=19, unit='cm', dpi=320)



maam.hyzounthin.model.hyzo <- 
      ggplot() + 
      geom_spatraster(data = hm.a.u) +
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'GnBu'), guide="colorbar", 
                           na.value="white",  name='', limits=0:1) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      # geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = suit.leg.pos, 
            axis.text = element_blank(),
            legend.box = "vertical",
            legend.justification = "top",
            legend.margin = margin(0.25, 0.9, 0.25, 0.25, unit = "cm"), # top, right, bottom, left
            legend.box.background = element_rect(fill = "white", colour = 'black'),
            legend.key.width = unit(1.2, 'cm' ), 
            legend.spacing.y = unit(0.2, "lines"),
            legend.text = element_text(size=16), 
            axis.ticks = element_blank(),
            legend.title = element_text(size = 20, colour = "black", face = "bold"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA)
      ) + 
      labs(x = NULL, y = NULL)
maam.hyzounthin.model.hyzo
ggsave('maam.hyzounthin.model.hyzo.png', height=20, width=19, unit='cm', dpi=320)



maca.hyzounthin.model.hyzo <- 
      ggplot() + 
      geom_spatraster(data = hm.c.u) +
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'GnBu'), guide="colorbar", 
                           na.value="white",  name='', limits=0:1) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      # geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = suit.leg.pos, 
            axis.text = element_blank(),
            legend.box = "vertical",
            legend.justification = "top",
            legend.margin = margin(0.25, 0.9, 0.25, 0.25, unit = "cm"), # top, right, bottom, left
            legend.box.background = element_rect(fill = "white", colour = 'black'),
            legend.key.width = unit(1.2, 'cm' ), 
            legend.spacing.y = unit(0.2, "lines"),
            legend.text = element_text(size=16), 
            axis.ticks = element_blank(),
            legend.title = element_text(size = 20, colour = "black", face = "bold"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA)
      ) + 
      labs(x = NULL, y = NULL)
maca.hyzounthin.model.hyzo
ggsave('maca.hyzounthin.model.hyzo.png', height=20, width=19, unit='cm', dpi=320)


ggarrange(maam.hyzothin.model.hyzo, maca.hyzothin.model.hyzo, nrow=1, ncol=2, 
          labels=c('A', 'B'), font.label=list(size=24), label.x=0.9, label.y=0.99)
ggsave('martes.hyzomodel.thin.png', height=20, width=40, unit='cm', dpi=320)
ggsave('FigSI4.eps', height=20, width=40, unit='cm', dpi=1200)
 
ggarrange(maam.hyzounthin.model.hyzo, maca.hyzounthin.model.hyzo, nrow=1, ncol=2, 
          labels=c('A', 'B'), font.label=list(size=24), label.x=0.9, label.y=0.99)
ggsave('martes.hyzomodel.unthin.png', height=20, width=40, unit='cm', dpi=320)
ggsave('FigSI5.eps', height=20, width=40, unit='cm', dpi=1200)


### reading in the data from the range-wide models to make 
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/NoHyzo")

maam.rw.ev <- fread('Maam_NoBio12/maam.nohyzo.everything.nam1.csv', header=T)
maca.rw.ev <- fread('Maca/maca.nohyzo.everything.nam.csv', header=T)

rw.to.hyzo <- do.call('cbind', list(back.all[,1:9], maam.rw.ev, maca.rw.ev)) %>% as_tibble()
rw.to.hyzo <- rw.to.hyzo %>% filter(longitude > -120 & longitude < -105 & latitude > 40 & latitude < 53)


# americana
maam.rw <- rw.to.hyzo %>% dplyr::filter(Maam_thresh > 0) %>% dplyr::select(longitude, latitude, Maam_contin)
maam.rw <- vect(maam.rw, crs=wgs84, geom=c('longitude', 'latitude'))
maam.rw <- rasterize(maam.rw, bc[[1]], field='Maam_contin')
maam.rw <- terra::project(maam.rw, n.am.rast)
maam.rw[maam.rw > 1] <- 1; maam.rw[maam.rw < 0.41542521] <- 0; plot(maam.rw)
# caurina
maca.rw <- rw.to.hyzo %>% dplyr::filter(Maca_thresh > 0) %>% dplyr::select(longitude, latitude, Maca_contin)
maca.rw <- vect(maca.rw, crs=wgs84, geom=c('longitude', 'latitude'))
maca.rw <- rasterize(maca.rw, bc[[1]], field='Maca_contin')
maca.rw <- terra::project(maca.rw, n.am.rast)
maca.rw[maca.rw > 1] <- 1; maca.rw[maca.rw < 0.24855051] <- 0; plot(maca.rw)

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")


# americana suitability map
maam.rw.plot <- 
      ggplot() + 
      geom_spatraster(data = maam.rw) +
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'GnBu'), guide="colorbar", 
                           na.value="white",  name='', limits=0:1) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      # geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = suit.leg.pos, 
            axis.text = element_blank(),
            legend.box = "vertical",
            legend.justification = "top",
            legend.margin = margin(0.25, 0.9, 0.25, 0.25, unit = "cm"), # top, right, bottom, left
            legend.box.background = element_rect(fill = "white", colour = 'black'),
            legend.key.width = unit(1.2, 'cm' ), 
            legend.spacing.y = unit(0.2, "lines"),
            legend.text = element_text(size=16), 
            axis.ticks = element_blank(),
            legend.title = element_text(size = 20, colour = "black", face = "bold"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA)
      ) + 
      labs(x = NULL, y = NULL)
maam.rw.plot
ggsave('maam.rw.hyzo.png', height=20, width=19, unit='cm', dpi=320)



maca.rw.plot <- 
      ggplot() + 
      geom_spatraster(data = maca.rw) +
      scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'GnBu'), guide="colorbar", 
                           na.value="white",  name='', limits=0:1) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      # geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      # geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      # geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = suit.leg.pos, 
            axis.text = element_blank(),
            legend.box = "vertical",
            legend.justification = "top",
            legend.margin = margin(0.25, 0.9, 0.25, 0.25, unit = "cm"), # top, right, bottom, left
            legend.box.background = element_rect(fill = "white", colour = 'black'),
            legend.key.width = unit(1.2, 'cm' ), 
            legend.spacing.y = unit(0.2, "lines"),
            legend.text = element_text(size=16), 
            axis.ticks = element_blank(),
            legend.title = element_text(size = 20, colour = "black", face = "bold"),
            plot.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA)
      ) + 
      labs(x = NULL, y = NULL)
maca.rw.plot
ggsave('maca.rw.hyzo.png', height=20, width=19, unit='cm', dpi=320)


ggarrange(maam.rw.plot, maca.rw.plot, nrow=1, ncol=2, 
          labels=c('A', 'B'), font.label=list(size=24), label.x=0.9, label.y=0.99)
ggsave('martes.rw.png', height=20, width=40, unit='cm', dpi=320)
ggsave('FigSI1.eps', height=20, width=40, unit='cm', dpi=1200)

##### Hotelling's T^2 tests -------------------



### projecting the THINNED hybrid zone models to each of the americana and caurina occurrences in the hybrid zone
# using thinned models to be consistent with the range-wide models

## projections to americana
# from caurina
maca.hyzo.to.maam.thin <- maxent.everything(eval=maca.hyzo.thin.eval1, 
                                            means=maca.hyzo.means.thin1,
                                            everything=occ.maam_hyzo_thin, 
                                            predic=6:9, 
                                            thresh=maca.thin.thresh1, 
                                            name='maca_to_maam')
# from americana
maam.hyzo.to.maam.thin <- maxent.everything(eval=maam.hyzo.thin.eval1, 
                                            means=maam.hyzo.means.thin1,
                                            everything=occ.maam_hyzo_thin, 
                                            predic=6:9, 
                                            thresh=maam.thin.thresh1, 
                                            name='maam_to_maam')

## projections to caurina
# from caurina
maca.hyzo.to.maca.thin <- maxent.everything(eval=maca.hyzo.thin.eval1, 
                                            means=maca.hyzo.means.thin1,
                                            everything=occ.maca_hyzo_thin, 
                                            predic=6:9, 
                                            thresh=maca.thin.thresh1, 
                                            name='maca_to_maca')
# from americana
maam.hyzo.to.maca.thin <- maxent.everything(eval=maam.hyzo.thin.eval1, 
                                            means=maam.hyzo.means.thin1,
                                            everything=occ.maca_hyzo_thin, 
                                            predic=6:9, 
                                            thresh=maam.thin.thresh1, 
                                            name='maam_to_maca')

## projections to hybrids
# from caurina
maca.hyzo.to.mahy.thin <- maxent.everything(eval=maca.hyzo.thin.eval1, 
                                            means=maca.hyzo.means.thin1,
                                            everything=occ.mahy_hyzo_thin, 
                                            predic=6:9, 
                                            thresh=maca.thin.thresh1, 
                                            name='maca_to_mahy')
# from americana
maam.hyzo.to.mahy.thin <- maxent.everything(eval=maam.hyzo.thin.eval1, 
                                            means=maam.hyzo.means.thin1,
                                            everything=occ.mahy_hyzo_thin, 
                                            predic=6:9, 
                                            thresh=maam.thin.thresh1, 
                                            name='maam_to_mahy')




## binding together
# to americana
to.maam.thin <- do.call('cbind', list(maca.hyzo.to.maam.thin, maam.hyzo.to.maam.thin))

# to caurina
to.maca.thin <- do.call('cbind', list(maca.hyzo.to.maca.thin, maam.hyzo.to.maca.thin))

# to hybrids
to.mahy.thin <- do.call('cbind', list(maca.hyzo.to.mahy.thin, maam.hyzo.to.mahy.thin))

hotel.thin.hyzo.test <- ICSNP::HotellingsT2(to.maam.thin[,c(1,3)], 
                                            to.maca.thin[,c(1,3)])
hotel.thin.hyzo.test
#         Hotelling's two sample T2-test
# data:  to.maam.thin[, c(1, 3)] and to.maca.thin[, c(1, 3)]
# T.2 = 5.0425, df1 = 2, df2 = 37, p-value = 0.01157
# alternative hypothesis: true location difference is not equal to c(0,0)

## fairly reasonable evidence that the multivariate niches of the species are not the same
# sample sizes for maam/maca/mahy are 13/27/19 respectively



### projecting the UNTHINNED hybrid zone models to each of the americana and caurina occurrences in the hybrid zone
# using thinned models to be consistent with the range-wide models

## projections to americana
# from caurina
maca.hyzo.to.maam.unthin <- maxent.everything(eval=maca.hyzo.unthin.eval1, 
                                              means=maca.hyzo.means.unthin1,
                                              everything=occ.maam_hyzo_unthin, 
                                              predic=6:9, 
                                              thresh=maca.unthin.thresh1, 
                                              name='maca_to_maam')
# from americana
maam.hyzo.to.maam.unthin <- maxent.everything(eval=maam.hyzo.unthin.eval1, 
                                              means=maam.hyzo.means.unthin1,
                                              everything=occ.maam_hyzo_unthin, 
                                              predic=6:9, 
                                              thresh=maam.unthin.thresh1, 
                                              name='maam_to_maam')

## projections to caurina
# from caurina
maca.hyzo.to.maca.unthin <- maxent.everything(eval=maca.hyzo.unthin.eval1, 
                                              means=maca.hyzo.means.unthin1,
                                              everything=occ.maca_hyzo_unthin, 
                                              predic=6:9, 
                                              thresh=maca.unthin.thresh1, 
                                              name='maca_to_maca')
# from americana
maam.hyzo.to.maca.unthin <- maxent.everything(eval=maam.hyzo.unthin.eval1, 
                                              means=maam.hyzo.means.unthin1,
                                              everything=occ.maca_hyzo_unthin, 
                                              predic=6:9, 
                                              thresh=maam.unthin.thresh1, 
                                              name='maam_to_maca')


## projections to hybrids
# from caurina
maca.hyzo.to.mahy.unthin <- maxent.everything(eval=maca.hyzo.unthin.eval1, 
                                              means=maca.hyzo.means.unthin1,
                                              everything=occ.mahy_hyzo_unthin, 
                                              predic=6:9, 
                                              thresh=maca.unthin.thresh1, 
                                              name='maca_to_mahy')
# from americana
maam.hyzo.to.mahy.unthin <- maxent.everything(eval=maam.hyzo.unthin.eval1, 
                                              means=maam.hyzo.means.unthin1,
                                              everything=occ.mahy_hyzo_unthin, 
                                              predic=6:9, 
                                              thresh=maam.unthin.thresh1, 
                                              name='maam_to_mahy')


## binding together
# to americana
to.maam.unthin <- do.call('cbind', list(maca.hyzo.to.maam.unthin, maam.hyzo.to.maam.unthin))

# to caurina
to.maca.unthin <- do.call('cbind', list(maca.hyzo.to.maca.unthin, maam.hyzo.to.maca.unthin))

# to hybrids
to.mahy.unthin <- do.call('cbind', list(maca.hyzo.to.mahy.unthin, maam.hyzo.to.mahy.unthin))

hotel.unthin.hyzo.test <- ICSNP::HotellingsT2(to.maam.unthin[,c(1,3)], 
                                              to.maca.unthin[,c(1,3)])
hotel.unthin.hyzo.test
#         Hotelling's two sample T2-test
# data:  to.maam.unthin[, c(1, 3)] and to.maca.unthin[, c(1, 3)]
# T.2 = 68.882, df1 = 2, df2 = 149, p-value < 2.2e-16
# alternative hypothesis: true location difference is not equal to c(0,0)

## reasonable evidence that the multivariate niches of the species are not the same
# sample sizes for maam/maca/mahy are 36/116/44 respectively

##### Bivariate Suitability Plots -------------------------------------


# making a data frame of points projected into the hybrid zone
# to be used for plotting later
hyzo.hotelling.df <- do.call('cbind', list(
      back.hyzo[,1:9], 
      maam.hyzo.thin.everything, maca.hyzo.thin.everything,
      maam.hyzo.unthin.everything, maca.hyzo.unthin.everything
))

colnames(hyzo.hotelling.df) <- c(colnames(hyzo.hotelling.df[,1:9]),
                                 c('Maam_contin.t',
                                   'Maam_thresh.t',
                                   'Maca_contin.t',
                                   'Maca_thresh.t',
                                   'Maam_contin.u',
                                   'Maam_thresh.u',
                                   'Maca_contin.u',
                                   'Maca_thresh.u' )
)
hyzo.hotelling.df <- as_tibble(hyzo.hotelling.df)

colnames(hyzo.hotelling.df)
#  [1] "Name"          "longitude"     "latitude"      "time.bin"      "region"        "bio4"          "bio10"         "bio12"         "bio15"         "Maam_contin.t"
# [11] "Maam_thresh.t" "Maca_contin.t" "Maca_thresh.t" "Maam_contin.u" "Maam_thresh.u" "Maca_contin.u" "Maca_thresh.u"

PerformanceAnalytics::chart.Correlation(hyzo.hotelling.df[,10:17])




setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")



# thinned models
bvs.maam.hyzo.thin <- ggplot() + 
      geom_point(data=hyzo.hotelling.df, aes(x=Maam_contin.t, y=Maca_contin.t), size=0.5) + # background points
      geom_point(data=to.maam.thin, aes(x=Maam_contin, y=Maca_contin), col='red', size=3) + # americana thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.maca.hyzo.thin <- ggplot() + 
      geom_point(data=hyzo.hotelling.df, aes(x=Maam_contin.t, y=Maca_contin.t), size=0.5) + # background points
      geom_point(data=to.maca.thin, aes(x=Maam_contin, y=Maca_contin), col='blue', size=3) + # caurina thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.mahy.hyzo.thin <- ggplot() + 
      geom_point(data=hyzo.hotelling.df, aes(x=Maam_contin.t, y=Maca_contin.t), size=0.5) + # background points
      geom_point(data=to.mahy.thin, aes(x=Maam_contin, y=Maca_contin), col='purple', size=3) + # hybrids thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

ggarrange(bvs.maam.hyzo.thin, bvs.maca.hyzo.thin, bvs.mahy.hyzo.thin, nrow=1, ncol=3, labels=c('A', 'B', 'C'), label.x = 0.2)

ggsave('BivarSuitHyzoThin.eps', height=10, width=30, units='cm', dpi=1200)

# unthinned models
bvs.maam.hyzo.unthin <- ggplot() + 
      geom_point(data=hyzo.hotelling.df, aes(x=Maam_contin.u, y=Maca_contin.u), size=0.5) + # background points
      geom_point(data=to.maam.unthin, aes(x=Maam_contin, y=Maca_contin), col='red', size=3) + # americana unthinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - unthinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - unthinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.maca.hyzo.unthin <- ggplot() + 
      geom_point(data=hyzo.hotelling.df, aes(x=Maam_contin.u, y=Maca_contin.u), size=0.5) + # background points
      geom_point(data=to.maca.unthin, aes(x=Maam_contin, y=Maca_contin), col='blue', size=3) + # caurina unthinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - unthinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - unthinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.mahy.hyzo.unthin <- ggplot() + 
      geom_point(data=hyzo.hotelling.df, aes(x=Maam_contin.u, y=Maca_contin.u), size=0.5) + # background points
      geom_point(data=to.mahy.unthin, aes(x=Maam_contin, y=Maca_contin), col='purple', size=3) + # hybrids unthinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - unthinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - unthinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

ggarrange(bvs.maam.hyzo.unthin, bvs.maca.hyzo.unthin, bvs.mahy.hyzo.unthin, nrow=1, ncol=3, labels=c('A', 'B', 'C'), label.x = 0.2)
ggsave('FigSI6.eps', height=10, width=30, units='cm', dpi=1200)
ggsave('BivarSuitHyzoUnthin.png', height=10, width=30, units='cm')








### showing the bivariate suitability space for areas North or South of I-90

hyzo.north <- hyzo.hotelling.df |> filter(region == 'HybridZoneNorth')
hyzo.south <- hyzo.hotelling.df |> filter(region == 'HybridZoneSouth')



ggnorth.u <- ggplot() + 
      geom_point(data=hyzo.north, aes(x=Maam_contin.u, y=Maca_contin.u), size=0.5) + # background points
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - unthinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - unthinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

ggsouth.u <- ggplot() + 
      geom_point(data=hyzo.south, aes(x=Maam_contin.u, y=Maca_contin.u), size=0.5) + # background points
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - unthinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - unthinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )


ggarrange(ggnorth.u, ggsouth.u, nrow=1, ncol=2, labels=c('A', 'B'), label.x = 0.2)

ggsave('BivarSuitHyzoNorthSouth.u.png', height=10, width=20, units='cm')






ggnorth.t <- ggplot() + 
      geom_point(data=hyzo.north, aes(x=Maam_contin.t, y=Maca_contin.t), size=0.5) + # background points
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

ggsouth.t <- ggplot() + 
      geom_point(data=hyzo.south, aes(x=Maam_contin.t, y=Maca_contin.t), size=0.5) + # background points
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )


ggarrange(ggnorth.t, ggsouth.t, nrow=1, ncol=2, labels=c('A', 'B'), label.x = 0.2)

ggsave('BivarSuitHyzoNorthSouth.t.png', height=10, width=20, units='cm')

 ##### END ------------------------------------------------------------