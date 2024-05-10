
##### Martes Range-Wide Models ------------------------------------------------

# Range-Wide Models of both Martes americana and M. caurina


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

library(RColorBrewer)
# optionally writing the package versions as a csv file
# write.csv(packs.df 'package.versions.csv')

# test if the maxent.jar file is in the dismo folder
list.files( system.file('java', package = 'dismo') )
# [1] "dismo.jar"  "maxent.jar"

# sort( sapply(ls(),function(x){object.size(get(x))})) 

# usethis::edit_r_environ()

# Analyses ran on 2023 16" MacBook Pro with M2 Chip with 96GB RAM

##### Spatial Data Prep ----------------------------


# WGS-1984 CRS
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

n.am.crs <- '+proj=aea +lon_0=-107.5 +lat_1=37.5 +lat_2=67.5 +lat_0=52.5 +datum=WGS84 +units=m +no_defs'
n.am.wkt <- 'PROJCS["ProjWiz_Custom_Albers",
                   GEOGCS["GCS_WGS_1984",
                          DATUM["D_WGS_1984",
                                SPHEROID["WGS_1984",6378137.0,298.257223563]
                                ],
                          PRIMEM["Greenwich",0.0],
                          UNIT["Degree",0.0174532925199433]],
                   PROJECTION["Albers"],
                   PARAMETER["False_Easting",0.0],
                   PARAMETER["False_Northing",0.0],
                   PARAMETER["Central_Meridian",-107.5],
                   PARAMETER["Standard_Parallel_1",37.5],
                   PARAMETER["Standard_Parallel_2",67.5],
                   PARAMETER["Latitude_Of_Origin",52.5],
                   UNIT["Meter",1.0] ]'

library(rnaturalearth)


### ggraster

## converts a raster object into a ggplot 

ggraster <- function(rast, x='long', y='lat', z='z', na.value = NULL){
      require(sp)
      require(raster)
      require(ggplot2)
      
      # making the data.frame
      r.values <- as.data.frame(rast@data@values)
      if( nrow(r.values) == 0){
            r.values <- as.data.frame(sampleRandom(rast, ncell(rast) ) )
      }
      df <- as.data.frame( cbind( coordinates(rast), r.values ) )
      colnames(df) <- c('x', 'y', 'z')
      if(!is.null(na.value)){
            df <- df %>% filter(df[,3] != na.value )
      }
      
      ggplot(df, aes(x=x, y=y) ) + geom_raster(aes(fill=z)) + xlab(x) + ylab(y) + labs(fill=z)
      
}

#  making a north american

n.am.rast <- raster(xmn=-4000000, xmx=4000000, ymn=-3000000, ymx=3000000, res=10000, crs=n.am.crs)
n.am.rast[n.am.rast] <- rnorm(ncell(n.am.rast))
plot(n.am.rast)


conts <- ne_coastline(50, "sf")
conts <- st_transform(conts, crs = n.am.crs)




##### Data Input -------------------------------------------------------

# working directory
setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP")
# sourcing all the custom functions needed to run the analyses
source('./Martes_R/NAF_Functions.R')

### reading in occurrence files

## range-wide occ data

# occurrence points for range-wide martes americana (maam)
# thinned to 25 km over 500 reps using spThin
occ.maam_nohyzo <- read.csv('./Martes_masterFiles/occ.maam_nohyzo.csv', header=T)

# occurrence points for range-wide martes caurina (maca)
# thinned to 25 km over 500 reps using spThin
occ.maca_nohyzo <- read.csv('./Martes_masterFiles/occ.maca_nohyzo.csv', header=T)

## hybrid zone occ data

# occurrence points for maam inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maam_hyzo_thin <- read.csv('./Martes_masterFiles/occ.maam_hyzo.thin.csv', header=T)

# occurrence points for maam inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maca_hyzo_thin <- read.csv('./Martes_masterFiles/occ.maca_hyzo.thin.csv', header=T)

# occurrence points for martes hybrids (mahy) inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.mahy_hyzo_thin <- read.csv('./Martes_masterFiles/occ.mahy_hyzo.thin.csv', header=T)

# occurrence points for maam inside the hybrid zone
# not thinned, but reduced to one occ per  grid cell per taxa
occ.maam_hyzo_unthin <- read.csv('./Martes_masterFiles/occ.maam_hyzo.unthin.csv', header=T)

# occurrence points for maam inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maca_hyzo_unthin <- read.csv('./Martes_masterFiles/occ.maca_hyzo.unthin.csv', header=T)

# occurrence points for martes hybrids (mahy) inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.mahy_hyzo_unthin <- read.csv('./Martes_masterFiles/occ.mahy_hyzo.unthin.csv', header=T)


### reading in background files

# background of all of n. america - used for projecting models after model construction
back.all <- read.csv('./Martes_masterFiles/back.all.csv', header=T)

# background of the hybrid zone
# 100km buffer around a maximum convex polygon of all enotped individuals with admixed dna
back.hyzo <- read.csv('./Martes_masterFiles/back.hyzo.csv', header=T)

# background used for maam models
# 250km buffer around thinned maam points
# excludes grid cells inside the hybrid zone
back.maam <- read.csv('./Martes_masterFiles/back.maam_nohyzo.csv', header=T)

# background used for maca models
# 250km buffer around thinned maca points
# excludes grid cells inside the hybrid zone
back.maca <- read.csv('./Martes_masterFiles/back.maca_nohyzo.csv', header=T)


##### Americana range-wide model --------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/NoHyzo/Maam")
# creating null data frame
maam.nohyzo.null <- create.null.df('Maam', 'modern', 'nohyzo')

# creating summary data.frame
maam.nohyzo.summary  <- create.summary.df('Maam', 'modern', 'nohyzo')

# creating folders where the model output will be stored
create.folders.for.maxent(maam.nohyzo.summary)

# run null model
maam.nohyzo.null <- null.aic(null.df = maam.nohyzo.null,
                        occs = occ.maam_nohyzo,
                        background = back.maam,
                        first.occ.col = 10); maam.nohyzo.null
write.csv(maam.nohyzo.null, 'maam.nohyzo.null.csv', row.names=F)

# running the range-wide models for maam
maam.nohyzo.optim <- optimize.maxent.likelihood(maam.nohyzo.summary,    # name of the summary file
                                                occs = occ.maam_nohyzo, # occurrences
                                                background = back.maam, # background
                                                predic = 6:9,           # column number of predictor variables
                                                first.occ.col = 10)     # first occurrence column

View(maam.nohyzo.optim)
beep(2)
write.csv(maam.nohyzo.optim, 'maam.nohyzo.optim.csv', row.names = F)

# the best performing model (as inferred from AICc; LQP-0.05) has reasonable-looking response curves
# and has lower AICc values than the null model


# making the evaluation data.frame for maam and creating the requisite folders
maam.nohyzo.eval <- create.eval.df('Maam', 'modern', 'nohyzo', beta.values=0.05, f.class='LQP'); maam.nohyzo.eval
create.folders.for.maxent(maam.nohyzo.eval)

# running the evaluation object
maam.nohyzo.eval <- maxent.crossval.error(eval.df=maam.nohyzo.eval, # evaluation df
                                          occs = occ.maam_nohyzo,   # occurrences
                                          background = back.maam,   # background
                                          predic = 6:9,             # column number of predictor variables
                                          first.occ.col = 11,       # first training column
                                          first.test.col = 16,      # first testing column
                                          omission.rate=0.1,        # omission rate
                                          all.background=back.maam
                                          )
View(maam.nohyzo.eval$summary)
write.csv(maam.nohyzo.eval$summary, 'maam.nohyzo.eval.csv', row.names=F)

beep(2)
Sys.time()

# weighting the cross-validation reps and projecting to all occ and background points
maam.nohyzo.means <- maxent.eval(eval=maam.nohyzo.eval)
nrow(maam.nohyzo.means$occ)
# [1] 200
# picking a quantile of the occs to define a threshold
# 0% = lowest training presence (LTP; sensu Pearson et al. 2007)
quantile(maam.nohyzo.means$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))

# the lowest 10% occurrence threshold
maam.thresh <- 0.33826309

suit.uncert.plot(maam.nohyzo.means)
ggsave('maam.nohyzo.means.pdf')


### projecting model to all extents
# all of north america
maam.nohyzo.everything.nam <- maxent.everything(eval=maam.nohyzo.eval,
                                               means=maam.nohyzo.means,
                                               everything=back.all,
                                               predic=6:9,
                                               thresh=maam.thresh,
                                               name='maam_nohyzo')
write.csv(maam.nohyzo.everything.nam, 'maam.nohyzo.everything.nam.csv', row.names=F)

# projecting to the hybrid zone
maam.nohyzo.everything.hyzo <- maxent.everything(eval=maam.nohyzo.eval,
                                                means=maam.nohyzo.means,
                                                everything=back.hyzo,
                                                thresh=maam.thresh,
                                                predic=6:9,
                                                name='maam_nohyzo')
write.csv(maam.nohyzo.everything.hyzo, 'maam.nohyzo.everything.hyzo.csv', row.names=F)



### mess analysis
# tolarance.vector
maam.nohyzo.tol <- c(1,1,
                     0,1,
                     1,1,
                     1,1)
# all of north america
maam.nohyzo.mess.nam <- informed.mess(ref.extent=back.maam,
                                      mess.extent=back.all,
                                      coord.cols=2:3,
                                      predic=6:9,
                                      tolerance=maam.nohyzo.tol)
write.csv(maam.nohyzo.mess.nam, 'maam.nohyzo.mess.nam.csv', row.names=F)
# hybrid zone
maam.nohyzo.mess.hyzo <- informed.mess(ref.extent=back.maam,
                                       mess.extent=back.hyzo,
                                       coord.cols=2:3,
                                       predic=6:9,
                                       tolerance=maam.nohyzo.tol)
write.csv(maam.nohyzo.mess.hyzo, 'maam.nohyzo.mess.hyzo.csv', row.names=F)

# ### mop analysis
# # reference extent
# maam.nohyzo.mop.ref <- calculate.mop(ref.extent=back.maam,
#                                      predic=6:9)
# # all of north america
# maam.nohyzo.mop.nam <- calculate.mop(ref.extent=back.maam,
#                                      new.extent=back.all,
#                                      predic=6:9)
# # all of north america
# maam.nohyzo.mop.hyzo <- calculate.mop(ref.extent=back.maam,
#                                       new.extent=back.hyzo,
#                                       predic=6:9)



##### Americana range-wide model w/o bio12 --------------------------

### when plotting the response curves with uncertainties, there is a problem with bio12
##  the response curve plot appears inverted, with high suitabilities at low and high MAP, but . . .
##  . . . low suitability in the middle (min suit at ~1500mm/year)
##  this is because there are only a few americana occurrences in cells with MAP ≥1500mm/yr, and . . .
#   . . . those cells are predominantly in southeast Alaska and coastal British Columbia
##  because there are only a small proportion of grid cells with MAP ≥1500mm/yr, this causes . . .
#   . . . an inverted response curve
##  the caurina response curve for bio12 is flat and has <2% contribution and importance to the model, so . . .
#   . . . this behavior is almost certainly a statistical artifact
##  americana likely doesn't really care what precipitation is
##  based on this, this section re-runs the americana model, but without bio12 in the model

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/NoHyzo/Maam_NoBio12")
# creating null data frame
maam.nohyzo.null1 <- create.null.df('Maam', 'modern', 'nohyzo')

# creating summary data.frame
maam.nohyzo.summary1  <- create.summary.df('Maam', 'modern', 'nohyzo')

# creating folders where the model output will be stored
create.folders.for.maxent(maam.nohyzo.summary1)

# run null model
maam.nohyzo.null1 <- null.aic(null.df = maam.nohyzo.null1,
                             occs = occ.maam_nohyzo,
                             background = back.maam,
                             first.occ.col = 10); maam.nohyzo.null
write.csv(maam.nohyzo.null1, 'maam.nohyzo.null.csv', row.names=F)

# running the range-wide models for maam
maam.nohyzo.optim1 <- optimize.maxent.likelihood(maam.nohyzo.summary1,  # name of the summary file
                                                occs = occ.maam_nohyzo, # occurrences
                                                background = back.maam, # background
                                                predic = c(6,7,9),      # column number of predictor variables
                                                first.occ.col = 10)     # first occurrence column

View(maam.nohyzo.optim1)
beep(2)
write.csv(maam.nohyzo.optim1, 'maam.nohyzo.optim1.csv', row.names = F)

# the best performing model (as inferred from AICc; LQP-0.25) has reasonable-looking response curves
# and has lower AICc values than the null model


# making the evaluation data.frame for maam and creating the requisite folders
maam.nohyzo.eval1 <- create.eval.df('Maam', 'modern', 'nohyzo', beta.values=0.25, f.class='LQP'); maam.nohyzo.eval1
create.folders.for.maxent(maam.nohyzo.eval1)

# running the evaluation object
maam.nohyzo.eval1 <- maxent.crossval.error(eval.df=maam.nohyzo.eval1, # evaluation df
                                          occs = occ.maam_nohyzo,   # occurrences
                                          background = back.maam,   # background
                                          predic = c(6,7,9),        # column number of predictor variables
                                          first.occ.col = 11,       # first training column
                                          first.test.col = 16,      # first testing column
                                          omission.rate=0.1,        # omission rate
                                          all.background=back.maam
)
View(maam.nohyzo.eval1$summary)
write.csv(maam.nohyzo.eval1$summary, 'maam.nohyzo.eval1.csv', row.names=F)

beep(2)
Sys.time()

# weighting the cross-validation reps and projecting to all occ and background points
maam.nohyzo.means1 <- maxent.eval(eval=maam.nohyzo.eval1)
nrow(maam.nohyzo.means1$occ)
# [1] 200
# picking a quantile of the occs to define a threshold
# 0% = lowest training presence (LTP; sensu Pearson et al. 2007)
quantile(maam.nohyzo.means1$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))

# the lowest 10% occurrence threshold
maam.thresh1 <- 0.41542521

suit.uncert.plot(maam.nohyzo.means1)
ggsave('maam.nohyzo.means1.pdf')


### projecting model to all extents
# all of north america
maam.nohyzo.everything.nam1 <- maxent.everything(eval=maam.nohyzo.eval1,
                                                means=maam.nohyzo.means1,
                                                everything=back.all,
                                                predic=c(6,7,9),
                                                thresh=maam.thresh1,
                                                name='maam_nohyzo1')
write.csv(maam.nohyzo.everything.nam1, 'maam.nohyzo.everything.nam1.csv', row.names=F)

# projecting to the hybrid zone
maam.nohyzo.everything.hyzo1 <- maxent.everything(eval=maam.nohyzo.eval1,
                                                 means=maam.nohyzo.means1,
                                                 everything=back.hyzo,
                                                 thresh=maam.thresh1,
                                                 predic=c(6,7,9),
                                                 name='maam_nohyzo1')
write.csv(maam.nohyzo.everything.hyzo1, 'maam.nohyzo.everything.hyzo1.csv', row.names=F)



### mess analysis
# tolarance.vector
maam.nohyzo.tol1 <- c(1,1,
                      1,1,
                      1,1)
# all of north america
maam.nohyzo.mess.nam1 <- informed.mess(ref.extent=back.maam,
                                      mess.extent=back.all,
                                      coord.cols=2:3,
                                      predic=c(6,7,9),
                                      tolerance=maam.nohyzo.tol1)
write.csv(maam.nohyzo.mess.nam1, 'maam.nohyzo.mess.nam1.csv', row.names=F)
# hybrid zone
maam.nohyzo.mess.hyzo1 <- informed.mess(ref.extent=back.maam,
                                       mess.extent=back.hyzo,
                                       coord.cols=2:3,
                                       predic=c(6,7,9),
                                       tolerance=maam.nohyzo.tol1)
write.csv(maam.nohyzo.mess.hyzo1, 'maam.nohyzo.mess.hyzo1.csv', row.names=F)



##### Caurina range-wide model --------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Models/NoHyzo/Maca")
# creating null data frame
maca.nohyzo.null <- create.null.df('Maca', 'modern', 'nohyzo')

# creating summary data.frame
maca.nohyzo.summary  <- create.summary.df('Maca', 'modern', 'nohyzo')

# creating folders where the model output will be stored
create.folders.for.maxent(maca.nohyzo.summary)

# run null model
maca.nohyzo.null <- null.aic(null.df = maca.nohyzo.null,
                             occs = occ.maca_nohyzo,
                             background = back.maca,
                             first.occ.col = 10); maca.nohyzo.null
write.csv(maca.nohyzo.null, 'maca.nohyzo.null.csv', row.names=F)

# running the range-wide models for maca
maca.nohyzo.optim <- optimize.maxent.likelihood(maca.nohyzo.summary,    # name of the summary file
                                                occs = occ.maca_nohyzo, # occurrences
                                                background = back.maca, # background
                                                predic = 6:9,           # column number of predictor variables
                                                first.occ.col = 10)     # first occurrence column

View(maca.nohyzo.optim)
beep(2)
write.csv(maca.nohyzo.optim, 'maca.nohyzo.optim.csv', row.names = F)

# the best performing model (as inferred from AICc; Q-0.025) has reasonable-looking response curves
# and has lower AICc values than the null model


# making the evaluation data.frame for maca and creating the requisite folders
maca.nohyzo.eval <- create.eval.df('Maca', 'modern', 'nohyzo', beta.values=0.025, f.class='Q'); maca.nohyzo.eval
create.folders.for.maxent(maca.nohyzo.eval)

# running the evaluation object
maca.nohyzo.eval <- maxent.crossval.error(eval.df=maca.nohyzo.eval, # evaluation df
                                          occs = occ.maca_nohyzo,   # occurrences
                                          background = back.maca,   # background
                                          predic = 6:9,             # column number of predictor variables
                                          first.occ.col = 11,       # first training column
                                          first.test.col = 16,      # first testing column
                                          omission.rate=0.1,        # omission rate
                                          all.background=back.maca
)
View(maca.nohyzo.eval$summary)
write.csv(maca.nohyzo.eval$summary, 'maca.nohyzo.eval.csv', row.names=F)
maca.nohyzo.eval$summary


beep(2)
Sys.time()

# weighting the cross-validation reps and projecting to all occ and background points
maca.nohyzo.means <- maxent.eval(eval=maca.nohyzo.eval)
nrow(maca.nohyzo.means$occ)
# [1] 70
# picking a quantile of the occs to define a threshold
# 0% = lowest training presence (LTP; sensu Pearson et al. 2007)
quantile(maca.nohyzo.means$occ$w.mean, probs=c(0, 0.01, 0.05, 0.1, 0.20, 0.25))

# the lowest 10% occurrence threshold
maca.thresh <- 0.24855051

suit.uncert.plot(maca.nohyzo.means)
ggsave('maca.nohyzo.means.pdf')


### projecting model to all extents
# all of north america
maca.nohyzo.everything.nam <- maxent.everything(eval=maca.nohyzo.eval,
                                                means=maca.nohyzo.means,
                                                everything=back.all,
                                                predic=6:9,
                                                thresh=maca.thresh,
                                                name='maca_nohyzo')
write.csv(maca.nohyzo.everything.nam, 'maca.nohyzo.everything.nam.csv', row.names=F)

# projecting to he hybrid zone
maca.nohyzo.everything.hyzo <- maxent.everything(eval=maca.nohyzo.eval,
                                                 means=maca.nohyzo.means,
                                                 everything=back.hyzo,
                                                 thresh=maca.thresh,
                                                 predic=6:9,
                                                 name='maca_nohyzo')
write.csv(maca.nohyzo.everything.hyzo, 'maca.nohyzo.everything.hyzo.csv', row.names=F)

### mess analysis
# tolarance.vector
maca.nohyzo.tol <- c(1,1,
                     0,1,
                     1,1,
                     1,0)
# all of north america
maca.nohyzo.mess.nam <- informed.mess(ref.extent=back.maca,
                                      mess.extent=back.all,
                                      coord.cols=2:3,
                                      predic=6:9,
                                      tolerance=maca.nohyzo.tol)
write.csv(maca.nohyzo.mess.nam, 'maca.nohyzo.mess.nam.csv', row.names=F)
# hybrid zone
maca.nohyzo.mess.hyzo <- informed.mess(ref.extent=back.maca,
                                       mess.extent=back.hyzo,
                                       coord.cols=2:3,
                                       predic=6:9,
                                       tolerance=maca.nohyzo.tol)
write.csv(maca.nohyzo.mess.hyzo, 'maca.nohyzo.mess.hyzo.csv', row.names=F)

### mop analysis
# reference extent
# maca.nohyzo.mop.ref <- calculate.mop(ref.extent=back.maca,
#                                      predic=6:9)
# Time to calculate distance matrix: 39.45038 secs
# Time to run MOP: 12.40025 mins
# Time to run function: 15.69368 mins
# Sys.time()
# beep(2)
# write.csv(maca.nohyzo.mop.ref, 'maca.nohyzo.mop.ref.csv', row.names=F)
# summary(maca.nohyzo.mop.ref)
#      pt1                pt5               pt10            pt1perc          pt5perc      
# Min.   :0.002994   Min.   :0.01073   Min.   :0.01359   Min.   :0.1344   Min.   :0.2827  
# 1st Qu.:0.037992   1st Qu.:0.06095   1st Qu.:0.07596   1st Qu.:0.2510   1st Qu.:0.4499  
# Median :0.057604   Median :0.08532   Median :0.10443   Median :0.3185   Median :0.5487  
# Mean   :0.065752   Mean   :0.09471   Mean   :0.11537   Mean   :0.3648   Mean   :0.6495  
# 3rd Qu.:0.084369   3rd Qu.:0.11964   3rd Qu.:0.14535   3rd Qu.:0.4485   3rd Qu.:0.7970  
# Max.   :0.580860   Max.   :0.74909   Max.   :0.85474   Max.   :1.6391   Max.   :2.6322  
# all of north america
# maca.nohyzo.mop.nam <- calculate.mop(ref.extent=back.maca,
#                                      new.extent=back.all,
#                                      predic=6:9)
# all of north america
# maca.nohyzo.mop.hyzo <- calculate.mop(ref.extent=back.maca,
#                                       new.extent=back.hyzo,
#                                       predic=6:9)
# write.csv(maca.nohyzo.mop.hyzo, 'maca.nohyzo.mop.hyzo.csv', row.names=F)
# beep(2)
# summary(maca.nohyzo.mop.hyzo)
#      pt1               pt5               pt10            pt1perc          pt5perc      
# Min.   :0.01052   Min.   :0.04132   Min.   :0.05456   Min.   :0.2172   Min.   :0.3534  
# 1st Qu.:0.09403   1st Qu.:0.12714   1st Qu.:0.14862   1st Qu.:0.3560   1st Qu.:0.5614  
# Median :0.14454   Median :0.17834   Median :0.20140   Median :0.4266   Median :0.6512  
# Mean   :0.17012   Mean   :0.20318   Mean   :0.22555   Mean   :0.4475   Mean   :0.6679  
# 3rd Qu.:0.22285   3rd Qu.:0.25399   3rd Qu.:0.27703   3rd Qu.:0.5178   3rd Qu.:0.7552  
# Max.   :0.75996   Max.   :0.79220   Max.   :0.82382   Max.   :1.0163   Max.   :1.4008  


##### Response Curve Plots --------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")

### americana response curves
get.var.contrib.importance(eval=maam.nohyzo.eval)
# $VarContribution
#          w.mean     w.sd
# bio10 54.593271 3.110879
# bio12 26.882622 2.714219
# bio15  5.948312 2.282906
# bio4  12.575792 2.605966
# 
# $VarImportance
#         w.mean     w.sd
# bio10 41.26852 2.728782
# bio12 18.00279 1.315015
# bio15 14.36389 2.858370
# bio4  26.36484 2.797013
maam.nohyzo.ranges <- get.var.ranges(maam.nohyzo.eval)
maam.nohyzo.ranges
#         bio4      bio10 bio12     bio15
# min    0.000 -0.2371667   108  5.470546
# max 1786.079 23.1412792  5353 96.510367
maam.nohyzo.rc <- get.maxent.response.curves(eval=maam.nohyzo.eval, expand=0)
maam.nohyzo.rc.poly <- process.maxent.response.curves(maam.nohyzo.rc)



### americana1 response curves
get.var.contrib.importance(eval=maam.nohyzo.eval1)
# $VarContribution
#          w.mean      w.sd
# bio10 71.749284 3.0314393
# bio15  8.565683 0.8264535
# bio4  19.685073 2.7942611
# 
# $VarImportance
#         w.mean     w.sd
# bio10 54.28250 3.748905
# bio15 18.24920 1.595225
# bio4  27.46828 3.310155
maam.nohyzo.ranges1 <- get.var.ranges(maam.nohyzo.eval1)
maam.nohyzo.ranges1
#         bio4      bio10     bio15
# min    0.000 -0.2371667  5.470546
# max 1786.079 23.1412792 96.510367
maam.nohyzo.rc1 <- get.maxent.response.curves(eval=maam.nohyzo.eval1, expand=0)
maam.nohyzo.rc.poly1 <- process.maxent.response.curves(maam.nohyzo.rc1)



### caurina response curves
get.var.contrib.importance(eval=maca.nohyzo.eval)
# $VarContribution
#          w.mean     w.sd
# bio10 89.614530 2.107234
# bio12  1.865814 1.211975
# bio15  4.356362 1.727013
# bio4   4.163272 2.368148
# 
# $VarImportance
#         w.mean     w.sd
# bio10 91.0791016 2.3202787
# bio12  0.4678093 0.9238896
# bio15  7.3835017 2.8715792
# bio4   1.0695823 1.5464794
maca.nohyzo.ranges <- get.var.ranges(maca.nohyzo.eval)
maca.nohyzo.ranges
#          bio4     bio10 bio12      bio15
# min  160.5096  4.660833    56   5.609513
# max 1170.2495 36.032500  3426 108.915160
maca.nohyzo.rc <- get.maxent.response.curves(eval=maca.nohyzo.eval, expand=0)
maca.nohyzo.rc.poly <- process.maxent.response.curves(maca.nohyzo.rc)




# bio4
bio4 <- ggplot() +
      geom_polygon(data=maca.nohyzo.rc.poly$bio4, aes(x=x/100, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.nohyzo.rc.poly1$bio4, aes(x=x/100, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.nohyzo.rc$bio4$input/100, y=maca.nohyzo.rc$bio4$bio4_mean), col='blue', linewidth=1) + # caurina mean
      geom_line(aes(x=maam.nohyzo.rc1$bio4$input/100, y=maam.nohyzo.rc1$bio4$bio4_mean), col='red', linewidth=1) + # caurina mean
      geom_vline(xintercept = maca.nohyzo.ranges[1,1]/100, linetype="dashed", color = "blue", linewidth=.5) + # caurina lower bound
      geom_vline(xintercept = maca.nohyzo.ranges[2,1]/100, linetype="dashed", color = "blue", linewidth=.5) + # caurina upper bound
      geom_vline(xintercept = maam.nohyzo.ranges[1,1]/100, linetype="dashed", color = "red", linewidth=.5) + # americana lower bound
      geom_vline(xintercept = maam.nohyzo.ranges[2,1]/100, linetype="dashed", color = "red", linewidth=.5) + # americana upper bound
      theme_classic() + xlim(0,18) + ylim(0,1) + 
      xlab('Temp. Seasonality (°C)') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=5, y=.3, label=paste0('1.1 ± 1.5%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=8, y=.9, label=paste0('27.5 ± 3.3%'), color='red', size=3.5 ) # americana variable contribution
bio4

# bio10
bio10 <- ggplot() +
      geom_polygon(data=maca.nohyzo.rc.poly$bio10, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.nohyzo.rc.poly1$bio10, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.nohyzo.rc$bio10$input, y=maca.nohyzo.rc$bio10$bio10_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.nohyzo.rc1$bio10$input, y=maam.nohyzo.rc1$bio10$bio10_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = maca.nohyzo.ranges[1,2], linetype="dashed", color = "blue", size=.5) + # caurina lower bound
      geom_vline(xintercept = maca.nohyzo.ranges[2,2], linetype="dashed", color = "blue", size=.5) + # caurina upper bound
      geom_vline(xintercept = maam.nohyzo.ranges[1,2], linetype="dashed", color = "red", size=.5) + # americana lower bound
      geom_vline(xintercept = maam.nohyzo.ranges[2,2], linetype="dashed", color = "red", size=.5) + # americana upper bound
      theme_classic() + xlim(-2.5, 37.5) + ylim(0,1) + 
      xlab('Mean Summer Temp. (°C)') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=16, y=.975, label=paste0('91.1 ± 2.3%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=29.5, y=.35, label=paste0('54.2 ± 3.7%'), color='red', size=3.5 ) # americana variable contribution
bio10

# bio12
bio12 <- ggplot() +
      geom_polygon(data=maca.nohyzo.rc.poly$bio12, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      # geom_polygon(data=maam.nohyzo.rc.poly$bio12, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.nohyzo.rc$bio12$input, y=maca.nohyzo.rc$bio12$bio12_mean), col='blue', size=1) + # caurina mean
      # geom_line(aes(x=maam.nohyzo.rc$bio12$input, y=maam.nohyzo.rc$bio12$bio12_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = maca.nohyzo.ranges[1,3], linetype="dashed", color = "blue", size=.5) + # caurina lower bound
      geom_vline(xintercept = maca.nohyzo.ranges[2,3], linetype="dashed", color = "blue", size=.5) + # caurina upper bound
      geom_vline(xintercept = maam.nohyzo.ranges[1,3], linetype="dashed", color = "red", size=.5) + # americana lower bound
      geom_vline(xintercept = maam.nohyzo.ranges[2,3], linetype="dashed", color = "red", size=.5) + # americana upper bound
      theme_classic() + xlim(0,5400) + ylim(0,1) + 
      xlab(expression(paste("Annual Precipitation mm*", yr^{-1}))) + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=1500, y=.3, label=paste0('0.5 ± 0.9%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=4500, y=.3, label=paste0('0%'), color='red', size=3.5 ) # americana variable contribution
bio12

ggsave('Bio12_RC_Problem.pdf', height=15, width=15, units='cm', dpi=320)
# showing the distribution of bio12 environmental and occurrence data for range-wide Americana models 
# the areas that have "high Mean Annual Precipitation" (≥1500mm/yr) are areas in SE Alaska and coastal BC
ggplot() + 
      geom_density(data=back.maam, mapping=aes(x=bio12), color='black', lty=1 ) + 
      geom_density(data=occ.maam_nohyzo, mapping=aes(x=bio12), color='firebrick1' )  + 
      theme_classic() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab('Density') +
      scale_x_continuous(limits=c(0,5500), breaks=seq(0,5000,1000), labels=seq(0,5000,1000), 
                         name=expression(paste("Annual Precipitation cm*", yr^{-1})) )
ggsave('Bio12_Distrib.pdf', height=15, width=15, units='cm', dpi=320)
countries <- map_data("world")
highprecip <- 1500
high.precip.maam <- occ.maam_nohyzo %>% filter(bio12 >= highprecip)
high.precip.back <- back.all %>% filter(bio12 >= highprecip)
ggplot() + 
      geom_polygon(data=countries, mapping=aes(x=long, y=lat, group = group), col=NA, lwd=3, fill = "grey50") +
      xlim(-147.5,-117.5) + ylim(40,65) +
      coord_fixed(ratio=2) + theme_classic() + theme(legend.position = 'none') + 
      geom_raster(data=high.precip.back, mapping=aes(x=longitude, y=latitude, fill=bio12)) + 
      geom_point(data=high.precip.maam, mapping=aes(x=longitude, y=latitude, color='firebrick1' ))
ggsave('Bio12_HighPrecipAreas.png', height=15, width=15, units='cm', dpi=320)

# bio15
bio15 <- ggplot() +
      geom_polygon(data=maca.nohyzo.rc.poly$bio15, aes(x=x, y=y), alpha=0.25, fill='blue') + # caurina uncertainty
      geom_polygon(data=maam.nohyzo.rc.poly1$bio15, aes(x=x, y=y), alpha=0.25, fill='red') + # americana uncertainty
      geom_line(aes(x=maca.nohyzo.rc$bio15$input, y=maca.nohyzo.rc$bio15$bio15_mean), col='blue', size=1) + # caurina mean
      geom_line(aes(x=maam.nohyzo.rc1$bio15$input, y=maam.nohyzo.rc1$bio15$bio15_mean), col='red', size=1) + # americana mean
      geom_vline(xintercept = maca.nohyzo.ranges[1,4], linetype="dashed", color = "blue", size=.5) + # caurina lower bound
      geom_vline(xintercept = maca.nohyzo.ranges[2,4], linetype="dashed", color = "blue", size=.5) + # caurina upper bound
      geom_vline(xintercept = maam.nohyzo.ranges[1,4], linetype="dashed", color = "red", size=.5) + # americana lower bound
      geom_vline(xintercept = maam.nohyzo.ranges[2,4], linetype="dashed", color = "red", size=.5) + # americana upper bound
      theme_classic() + ylim(0,1) + 
      xlab('Precip. Seasonality') + ylab('Suitability') + # basic theme
      theme(aspect.ratio = 1, # more advanced theme
            axis.text = element_text(size=15, color='black'),
            axis.title = element_text(size=15, color='black') ) + 
      annotate('text', x=70, y=.05, label=paste0('7.4 ± 2.9%'), color='blue', size=3.5 ) + # caurina variable contribution
      annotate('text', x=70, y=.75, label=paste0('18.2 ± 1.6%'), color='red', size=3.5 ) # americana variable contribution
bio15


ggpubr::ggarrange(bio4, bio10, bio12, bio15, 
                  nrow=2, ncol=2, labels = c('A', 'B', 'C', 'D'), label.x=0.8, label.y=0.95, align='hv',
                  font.label=list(size=24, color="black", face="bold", family=NULL) ) + bgcolor('White')
ggsave('MartesResponseCurves.eps', width=20, height=20, units='cm', dpi=1200)


##### projecting into the hybrid zone ------------------------

### projecting the caurina range-wide model into various aspects of the hybrid zone
# hybrid zone
maca.nohyzo.everything.hyzo <- maxent.everything(eval=maca.nohyzo.eval,
                                                 means=maca.nohyzo.means,
                                                 everything=back.hyzo,
                                                 thresh=maca.thresh,
                                                 predic=6:9,
                                                 name='maca_rw_hyzo')
# thinned caurina occs in the hybrid zone
maca.nohyzo.everything.maca.thin <- maxent.everything(eval=maca.nohyzo.eval,
                                                      means=maca.nohyzo.means,
                                                      everything=occ.maca_hyzo_thin,
                                                      thresh=maca.thresh,
                                                      predic=6:9,
                                                      name='maca_rw_maca_thin')
# unthinned caurina occs in the hybrid zone
maca.nohyzo.everything.maca.unthin <- maxent.everything(eval=maca.nohyzo.eval,
                                                        means=maca.nohyzo.means,
                                                        everything=occ.maca_hyzo_thin,
                                                        thresh=maca.thresh,
                                                        predic=6:9,
                                                        name='maca_rw_maca_unthin')
# thinned americana occs in the hybrid zone
maca.nohyzo.everything.maam.thin <- maxent.everything(eval=maca.nohyzo.eval,
                                                      means=maca.nohyzo.means,
                                                      everything=occ.maam_hyzo_thin,
                                                      thresh=maca.thresh,
                                                      predic=6:9,
                                                      name='maca_rw_maam_thin')
# unthinned americana occs in the hybrid zone
maca.nohyzo.everything.maam.unthin <- maxent.everything(eval=maca.nohyzo.eval,
                                                        means=maca.nohyzo.means,
                                                        everything=occ.maam_hyzo_unthin,
                                                        thresh=maca.thresh,
                                                        predic=6:9,
                                                        name='maca_rw_maam_unthin')
# thinned hybrid occs
maca.nohyzo.everything.mahy.thin <- maxent.everything(eval=maca.nohyzo.eval,
                                                      means=maca.nohyzo.means,
                                                      everything=occ.mahy_hyzo_thin,
                                                      thresh=maca.thresh,
                                                      predic=6:9,
                                                      name='maca_rw_mahy_unthin')
# unthinned hyrbid occs
maca.nohyzo.everything.mahy.unthin <- maxent.everything(eval=maca.nohyzo.eval,
                                                        means=maca.nohyzo.means,
                                                        everything=occ.mahy_hyzo_unthin,
                                                        thresh=maca.thresh,
                                                        predic=6:9,
                                                        name='maca_rw_mahy_thin')




### projecting the americana range-wide model into various aspects of the hybrid zone
# hybrid zone
maam.nohyzo.everything.hyzo <- maxent.everything(eval=maam.nohyzo.eval1,
                                                 means=maam.nohyzo.means1,
                                                 everything=back.hyzo,
                                                 thresh=maam.thresh1,
                                                 predic=c(6,7,9),
                                                 name='maam_rw_hyzo')
# thinned caurina occs in the hybrid zone
maam.nohyzo.everything.maca.thin <- maxent.everything(eval=maam.nohyzo.eval1,
                                                      means=maam.nohyzo.means1,
                                                      everything=occ.maca_hyzo_thin,
                                                      thresh=maam.thresh1,
                                                      predic=c(6,7,9),
                                                      name='maam_rw_maca_thin')
# unthinned caurina occs in the hybrid zone
maam.nohyzo.everything.maca.unthin <- maxent.everything(eval=maam.nohyzo.eval1,
                                                        means=maam.nohyzo.means1,
                                                        everything=occ.maca_hyzo_thin,
                                                        thresh=maam.thresh1,
                                                        predic=c(6,7,9),
                                                        name='maam_rw_maca_unthin')
# thinned americana occs in the hybrid zone
maam.nohyzo.everything.maam.thin <- maxent.everything(eval=maam.nohyzo.eval1,
                                                      means=maam.nohyzo.means1,
                                                      everything=occ.maam_hyzo_thin,
                                                      thresh=maam.thresh1,
                                                      predic=c(6,7,9),
                                                      name='maam_rw_maam_thin')
# unthinned americana occs in the hybrid zone
maam.nohyzo.everything.maam.unthin <- maxent.everything(eval=maam.nohyzo.eval1,
                                                        means=maam.nohyzo.means1,
                                                        everything=occ.maam_hyzo_unthin,
                                                        thresh=maam.thresh1,
                                                        predic=c(6,7,9),
                                                        name='maam_rw_maam_unthin')
# thinned hybrid occs
maam.nohyzo.everything.mahy.thin <- maxent.everything(eval=maam.nohyzo.eval1,
                                                      means=maam.nohyzo.means1,
                                                      everything=occ.mahy_hyzo_thin,
                                                      thresh=maam.thresh1,
                                                      predic=c(6,7,9),
                                                      name='maam_rw_mahy_unthin')
# unthinned hyrbid occs
maam.nohyzo.everything.mahy.unthin <- maxent.everything(eval=maam.nohyzo.eval1,
                                                        means=maam.nohyzo.means1,
                                                        everything=occ.mahy_hyzo_unthin,
                                                        thresh=maam.thresh1,
                                                        predic=c(6,7,9),
                                                        name='maam_rw_mahy_thin')


### making data.frames of all these projections
# projection to the whole hybrid zone
proj.hyzo <- data.frame(maam=maam.nohyzo.everything.hyzo[,1], maca=maca.nohyzo.everything.hyzo[,1])
colnames(proj.hyzo) <- c('maam', 'maca')
# projection to americana occs - thinned
proj.maam.thin <- data.frame(maam=maam.nohyzo.everything.maam.thin[,1], maca=maca.nohyzo.everything.maam.thin[,1])
colnames(proj.maam.thin) <- c('maam', 'maca')
# projection to caurina occs - thinned
proj.maca.thin <- data.frame(maam=maam.nohyzo.everything.maca.thin[,1], maca=maca.nohyzo.everything.maca.thin[,1])
colnames(proj.maca.thin) <- c('maam', 'maca')
# projection to caurina occs - thinned
proj.mahy.thin <- data.frame(maam=maam.nohyzo.everything.mahy.thin[,1], maca=maca.nohyzo.everything.mahy.thin[,1])
colnames(proj.mahy.thin) <- c('maam', 'maca')
# projection to americana occs - unthinned
proj.maam.unthin <- data.frame(maam=maam.nohyzo.everything.maam.unthin[,1], maca=maca.nohyzo.everything.maam.unthin[,1])
colnames(proj.maam.unthin) <- c('maam', 'maca')
# projection to caurina occs - unthinned
proj.maca.unthin <- data.frame(maam=maam.nohyzo.everything.maca.unthin[,1], maca=maca.nohyzo.everything.maca.unthin[,1])
colnames(proj.maca.unthin) <- c('maam', 'maca')
# projection to caurina occs - unthinned
proj.mahy.unthin <- data.frame(maam=maam.nohyzo.everything.mahy.unthin[,1], maca=maca.nohyzo.everything.mahy.unthin[,1])
colnames(proj.mahy.unthin) <- c('maam', 'maca')


##### Hotelling's T^2 tests -------------------

# running Hotelling's T^2 tests to see how models projected into the hybrid zone 

## thinned data
hotel.thin.test <- ICSNP::HotellingsT2(proj.maam.thin, 
                                       proj.maca.thin)
hotel.thin.test
#         Hotelling's two sample T2-test
# data:  proj.maam.thin and proj.maca.thin
# T.2 = 1.4671, df1 = 2, df2 = 37, p-value = 0.2437
# alternative hypothesis: true location difference is not equal to c(0,0)

## unthinned data
hotel.unthin.test <- ICSNP::HotellingsT2(proj.maam.unthin, 
                                         proj.maca.unthin)
hotel.unthin.test
#         Hotelling's two sample T2-test
# data:  proj.maam.unthin and proj.maca.unthin
# 5.9183, df1 = 2, df2 = 60, p-value = 0.00451
# alternative hypothesis: true location difference is not equal to c(0,0)


### the test using models projected to the thinned occs shows that there doesn't appear to be any . . .
#   . . . broad scale pattern in terms of causing hybridization
##  but the test using models projected to the "unthinned" occs (but reduced to one occ/species.grid cell) . . .
#   suggests that there is some evidence that hybridization is important


setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")

# thinned models
ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maam.thin, aes(x=maam, y=maca), col='red') + # americana thinned
      geom_point(data=proj.maca.thin, aes(x=maam, y=maca), col='blue') + # caurina thinned
      geom_point(data=proj.mahy.thin, aes(x=maam, y=maca), col='purple') + # hybrids thinned
      xlim(c(0,1)) + ylim(c(0,1)) + theme_classic() + coord_equal() + 
      labs(x = 'Americana Suitability', y = 'Caurina Suitability')
ggsave('Proj_Hyzo_Thin.pdf', height=15, width=15, units='cm')

# unthinned models
ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maam.unthin, aes(x=maam, y=maca), col='red') + # americana thinned
      geom_point(data=proj.maca.unthin, aes(x=maam, y=maca), col='blue') + # caurina thinned
      geom_point(data=proj.mahy.unthin, aes(x=maam, y=maca), col='purple') + # hybrids thinned
      xlim(c(0,1)) + ylim(c(0,1)) + theme_classic() + coord_equal() + 
      labs(x = 'Americana Suitability', y = 'Caurina Suitability')
ggsave('Proj_Hyzo_Unthin.pdf', height=15, width=15, units='cm')


setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")



# thinned models
bvs.maam.rw.thin <- ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maam.thin, aes(x=maam, y=maca), col='red', size=3) + # americana thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.maca.rw.thin <- ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maca.thin, aes(x=maam, y=maca), col='blue', size=3) + # caurina thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.mahy.rw.thin <- ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maca.thin, aes(x=maam, y=maca), col='purple', size=3) + # hybrids thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

ggarrange(bvs.maam.rw.thin, bvs.maca.rw.thin, bvs.mahy.rw.thin, nrow=1, ncol=3, labels=c('A', 'B', 'C'), label.x = 0.2)
ggsave('FigSI2.eps', height=10, width=30, units='cm', dpi=1200)
ggsave('BivarSuitRangeWideThin.png', height=10, width=30, units='cm')

# unthinned models
bvs.maam.rw.unthin <- ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maam.unthin, aes(x=maam, y=maca), col='red', size=3) + # americana thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.maca.rw.unthin <- ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maca.unthin, aes(x=maam, y=maca), col='blue', size=3) + # caurina thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

bvs.mahy.rw.unthin <- ggplot() + 
      geom_point(data=proj.hyzo, aes(x=maam, y=maca), size=0.5) + # background points
      geom_point(data=proj.maca.unthin, aes(x=maam, y=maca), col='purple', size=3) + # hybrids thinned
      xlim(c(0,1)) + ylim(c(0,1)) + 
      xlab(expression(~italic('M. americana')~'suitability - thinned')) + 
      ylab(expression(~italic('M. caurina')~'suitability - thinned')) + 
      theme_classic() + coord_equal() + 
      theme(
            axis.text=element_text(size=12),
            axis.title=element_text(size=12)
      )

ggarrange(bvs.maam.rw.unthin, bvs.maca.rw.unthin, bvs.mahy.rw.unthin, nrow=1, ncol=3, labels=c('A', 'B', 'C'), label.x = 0.2)
ggsave('FigSI3.eps', height=10, width=30, units='cm', dpi=1200)
ggsave('BivarSuitRangeWideUnthin.png', height=10, width=30, units='cm')

##### Map Figures -----------------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_Results")

library(USA.state.boundaries)
library(canadianmaps)
library(tidyterra)

# occurrence points for maam inside the hybrid zone
# not thinned, but reduced to one occ per  grid cell per taxa
occ.maam_hyzo_unthin <- read.csv('../Martes_masterFiles/occ.maam_hyzo.unthin.csv', header=T)
occ.maam_hyzo_unthin_n <- occ.maam_hyzo_unthin %>% filter(region == 'HybridZoneNorth')
occ.maam_hyzo_unthin_n_sf <- st_as_sf(occ.maam_hyzo_unthin_n, coords=c('longitude', 'latitude'), crs=wgs84)
occ.maam_hyzo_unthin_n_sf <- st_transform(occ.maam_hyzo_unthin_n_sf, n.am.crs)
occ.maam_hyzo_unthin_s <- occ.maam_hyzo_unthin %>% filter(region == 'HybridZoneSouth')
occ.maam_hyzo_unthin_s_sf <- st_as_sf(occ.maam_hyzo_unthin_s, coords=c('longitude', 'latitude'), crs=wgs84)
occ.maam_hyzo_unthin_s_sf <- st_transform(occ.maam_hyzo_unthin_s_sf, n.am.crs)
table(occ.maam_hyzo_unthin$region)
# HybridZoneNorth HybridZoneSouth 
#              34               2 
# 5.56% of americana individuals are south of I-90


# occurrence points for martes hybrids (mahy) inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.mahy_hyzo_unthin <- read.csv('../Martes_masterFiles/occ.mahy_hyzo.unthin.csv', header=T)
occ.mahy_hyzo_unthin_n <- occ.mahy_hyzo_unthin %>% filter(region == 'HybridZoneNorth')
occ.mahy_hyzo_unthin_n_sf <- st_as_sf(occ.mahy_hyzo_unthin_n, coords=c('longitude', 'latitude'), crs=wgs84)
occ.mahy_hyzo_unthin_n_sf <- st_transform(occ.mahy_hyzo_unthin_n_sf, n.am.crs)
occ.mahy_hyzo_unthin_s <- occ.mahy_hyzo_unthin %>% filter(region == 'HybridZoneSouth')
occ.mahy_hyzo_unthin_s_sf <- st_as_sf(occ.mahy_hyzo_unthin_s, coords=c('longitude', 'latitude'), crs=wgs84)
occ.mahy_hyzo_unthin_s_sf <- st_transform(occ.mahy_hyzo_unthin_s_sf, n.am.crs)
table(occ.mahy_hyzo_unthin$region)
# HybridZoneNorth HybridZoneSouth 
#              19              25 
# the proportion of hybrids on either side of I-90 is much closer to 50% (43.18 north and 56.81% south)


# occurrence points for maca inside the hybrid zone
# thinned to 25 km over 500 reps using spThin
occ.maca_hyzo_unthin <- read.csv('../Martes_masterFiles/occ.maca_hyzo.unthin.csv', header=T)
occ.maca_hyzo_unthin_n <- occ.maca_hyzo_unthin %>% filter(region == 'HybridZoneNorth')
occ.maca_hyzo_unthin_n_sf <- st_as_sf(occ.maca_hyzo_unthin_n, coords=c('longitude', 'latitude'), crs=wgs84)
occ.maca_hyzo_unthin_n_sf <- st_transform(occ.maca_hyzo_unthin_n_sf, n.am.crs)
occ.maca_hyzo_unthin_s <- occ.maca_hyzo_unthin %>% filter(region == 'HybridZoneSouth')
occ.maca_hyzo_unthin_s_sf <- st_as_sf(occ.maca_hyzo_unthin_s, coords=c('longitude', 'latitude'), crs=wgs84)
occ.maca_hyzo_unthin_s_sf <- st_transform(occ.maca_hyzo_unthin_s_sf, n.am.crs)
table(occ.maca_hyzo_unthin$region)
# HybridZoneNorth HybridZoneSouth 
#              13             103 
# 11.21% of caurina individuals are north of I-90



# sf polygon of hybrid zone
hyzo.mask <- rast('../Martes_Shapefiles/hybrid.zone.mask.asc'); hyzo.mask[is.na(hyzo.mask)] <- 0
hyzo.mask <- st_contour(st_as_stars(hyzo.mask), contour_lines = TRUE)
hyzo.mask <- st_cast(st_transform(hyzo.mask, crs = n.am.crs), 'POLYGON')

# us interstates shapefile
# downloaded from https://hub.arcgis.com/maps/esri::usa-freeway-system/about
interstate <- read_sf('../Martes_Shapefiles/USA_Freeway_System/USA_Freeway_System.shp')
# interstate
interstate <- st_transform(interstate[interstate$ROUTE_NUM %in% c('I90', 'I15'),], n.am.crs)
interstate <- st_crop(interstate, st_buffer(hyzo.mask, 50000))
i90 <- interstate  %>% filter(ROUTE_NUM == 'I90')
# 
# hyzo.mask.split <- lwgeom::st_split(hyzo.mask, interstate)
# plot(hyzo.mask.split)
# hyzo.mask.split <- st_collection_extract(hyzo.mask.split, 'POLYGON')
# plot(hyzo.mask.split, col=factor(1:3))
# hyzo.mask.split[1,] <- st_union(hyzo.mask.split[1,], hyzo.mask.split[2,])
# hyzo.mask.split <- hyzo.mask.split[c(1,3),]
# # there's a ~10km long stretch in Butte, MT where I-90 is split by I-15
# # connecting the two ends of
# # identifying where the two ends 
# heads <- as.data.frame(do.call('rbind', lapply(i90$geometry[[1]], head,1)))
# heads[,3] <- LETTERS[1:10]; colnames(heads) <- c('x', 'y', 'line')
# tails <- as.data.frame(do.call('rbind', lapply(i90$geometry[[1]], tail,1)))
# tails[,3] <- LETTERS[1:10]; colnames(tails) <- c('x', 'y', 'line')
# heads; tails
# fuzzyjoin::distance_full_join(heads, tails, by=c('x', 'y'))
# #           x.x       y.x line.x        x.y       y.y line.y
# #  1  -67869.54 -763764.3      A  -67869.54 -763764.3      B
# #  2 -156865.00 -775959.0      B -156865.00 -775959.0      C
# #  3 -265568.02 -776038.1      C -265568.02 -776038.1      D
# #  4 -474216.62 -624622.5      E -474216.62 -624622.5      F
# #  5 -558489.19 -568388.0      F -558489.19 -568388.0      G
# #  6 -622570.03 -536849.5      G -622570.03 -536849.5      I
# #  7 -391373.92 -684328.9      H -391373.92 -684328.9      E
# #  8 -712780.80 -506901.6      I -712780.80 -506901.6      J
# #  9 -376493.32 -734235.5      D         NA        NA   <NA>
# # 10 -839812.72 -556437.7      J         NA        NA   <NA>
# # 11         NA        NA   <NA>  -19737.71 -773645.6      A
# # 12         NA        NA   <NA> -388231.75 -730274.8      H
# 
# # the linestrings are arranged from east to west, but the points within a linestring go west to east
# # linestrings D and E appear to be where the break is
# nrow(i90$geometry[[1]][[6]]); nrow(i90$geometry[[1]][[5]]); nrow(i90$geometry[[1]][[4]]); nrow(i90$geometry[[1]][[3]])
# # [1] 355
# # [1] 369
# # [1] 358
# # [1] 343
# i90$geometry[[1]][[6]][c(1,355),]
# #           [,1]      [,2]
# # [1,] -558489.2 -568388.0
# # [2,] -474216.6 -624622.5
# i90$geometry[[1]][[5]][c(1,369),]
# #           [,1]      [,2]
# # [1,] -474216.6 -624622.5
# # [2,] -391373.9 -684328.9
# i90$geometry[[1]][[4]][c(1,358),]
# #           [,1]      [,2]
# # [1,] -376493.3 -734235.5
# # [2,] -265568.0 -776038.1
# i90$geometry[[1]][[3]][c(1,343),]
# #         [,1]      [,2]
# # [1,] -265568 -776038.1
# # [2,] -156865 -775959.0
# 
# 
# 
# # the end of string E does not connect to the beginning of string D
# # connecting the end of string E to the beginning of string D
# i90$geometry[[1]][[5]][369,] <- i90$geometry[[1]][[4]][1,]


# rotated points in hybrid zone 

# sf polygon of americana M
maam.mask <- rast('../Martes_Shapefiles/maam.temp.asc'); maam.mask[is.na(maam.mask)] <- 0
maam.mask <- st_contour(st_as_stars(maam.mask), contour_lines = TRUE)
maam.mask <- st_transform(maam.mask, crs=n.am.crs)

# sf polygon of americana M
maca.mask <- rast('../Martes_Shapefiles/maca.temp.asc'); maca.mask[is.na(maca.mask)] <- 0
maca.mask <- st_contour(st_as_stars(maca.mask), contour_lines = TRUE)
maca.mask <- st_transform(maca.mask, crs=n.am.crs)

# sf polygon of coastlines
land <- ne_countries(50, returnclass = "sf")
land$Merge <- 1
land <- land %>% dplyr::group_by(Merge)
land <- st_transform(land, crs = n.am.crs) %>% filter(admin == 'Canada' | admin == 'United States of America')

rw.plots <- do.call('cbind', list(back.all[,1:9], maam.nohyzo.everything.nam1, maam.nohyzo.mess.nam1,
                                  maca.nohyzo.everything.nam, maca.nohyzo.mess.nam))
colnames(rw.plots) <- c("Name", "longitude", "latitude","time.bin", "region", "bio4", "bio10", "bio12", "bio15",
                        "Maam_contin", "Maam_thresh", "Maam_mess.raw", "Maam_mess.thresh", 
                        "Maca_contin", "Maca_thresh", "Maca_mess.raw", "Maca_mess.thresh")

# digital elevation model of study region
# from MERIT DEM - Yamazaki et al. 2017 GRL
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL072874
dem <- rast('../Martes_Shapefiles/DEM/NA_DEM.asc')
dem <- terra::project(dem, n.am.crs)
dem <- crop(x=dem, y=st_buffer(hyzo.mask, 50000))
dem[dem > 3500] <- 3500; dem[dem<500] <- 500
dem.df <- as.data.table(as.data.frame(dem, xy=T))
# dem.df$NA_DEM[dem.df$NA_DEM > 3000] <- 3000; dem.df$NA_DEM[dem.df$NA_DEM < 500] <- 500
# slope, aspect, and hillshade rasters
slope <- terrain(dem, 'slope', unit='radians')
aspect <- terrain(dem, 'aspect', unit='radians')
hill <- shade(slope, aspect, 30, 270)
# normalize names of hillshade
names(hill) <- "shades"
# Hillshading, but we need a palette
pal_greys <- hcl.colors(100, "Grays")[10:100]


# relevant US states and Canadian provinces
states <- state_boundaries_wgs84 %>% dplyr::filter(NAME %in% c('Colorado', 'Idaho', 'Montana', 
                                                               'Oregon', 'Washington', 'Wyoming') )
states <- st_transform(states, crs=n.am.crs)
provinces <- canadianmaps::PROV %>% filter(PRENAME %in% c('British Columbia', 'Alberta', 'Saskatchewan') )
provinces <- st_transform(provinces, crs=n.am.crs)

## projected occ corrdinates
occ.maam_nohyzo.proj <- PlotSvalbard::transform_coord(occ.maam_nohyzo, lon='longitude', lat='latitude',
                                                      proj.og=wgs84, proj.out=n.am.crs)
occ.maca_nohyzo.proj <- PlotSvalbard::transform_coord(occ.maca_nohyzo, lon='longitude', lat='latitude',
                                                      proj.og=wgs84, proj.out=n.am.crs)

# americana range-wide data.frame
maam.rw <- rw.plots %>% filter(Maam_mess.thresh > 0)
maam.rw <- vect(maam.rw, crs=wgs84, geom=c('longitude', 'latitude'))
maam.rw <- rasterize(maam.rw, bc[[1]], field='Maam_contin')
maam.rw <- terra::project(maam.rw, n.am.rast)
# thresholding
maam.rw[maam.rw > 1] <- 1; maam.rw[maam.rw < maam.thresh1] <- 0; plot(maam.rw)
# converting to data table
maam.rw <- as.data.table(as.data.frame(maam.rw, xy = TRUE, long = TRUE))
maam.rw <- maam.rw[complete.cases(maam.rw),]
maam.rw
# # cropping the range-wide projections to the hybrid zone
# maam.rw.to.hyzo <- terra::crop(x=maam.rw, y=st_buffer(hyzo.mask, 50000))

# americana suitability map
a <- ggplot() + geom_raster(data=maam.rw, aes(x=x, y=y, fill=last)) + 
      scale_fill_gradientn(colors=brewer.pal(9, 'GnBu'), guide="colorbar", na.value="white", 
                           name='', limits=0:1) + 
      coord_equal() + theme_classic() + 
      geom_sf(data=maam.mask, col = "black", fill = grey(level=0, alpha=0.2), linewidth=0.5, lty = 1) + 
      geom_sf(data=hyzo.mask, col = "black", fill = NA, linewidth=0.5, lty = 2) + 
      geom_sf(data=land, fill='NA', col=grey(level=0, alpha=1)) + 
      geom_point(data=occ.maam_nohyzo.proj, aes(x=lon.utm, y=lat.utm), fill='red', col='black', size=1, pch=21) + 
      xlim(-3250000, 3750000) + ylim(-2250000, 2750000) + 
      theme(legend.position = c(0.15, 0.4),
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.line=element_blank(), axis.title=element_blank())


# caurina range-wide data.frame
maca.rw <- rw.plots %>% filter(Maca_mess.thresh > 0)
maca.rw <- vect(maca.rw, crs=wgs84, geom=c('longitude', 'latitude'))
maca.rw <- rasterize(maca.rw, bc[[1]], field='Maca_contin')
maca.rw <- terra::project(maca.rw, n.am.rast)
# thresholding
maca.rw[maca.rw > 1] <- 1; maca.rw[maca.rw < maca.thresh] <- 0; plot(maca.rw)
# converting to data table
maca.rw <- as.data.table(as.data.frame(maca.rw, xy = TRUE, long = TRUE))
maca.rw <- maca.rw[complete.cases(maca.rw),]
maca.rw

# caurina suitability map
b <- ggplot() + geom_raster(data=maca.rw, aes(x=x, y=y, fill=last)) + 
      scale_fill_gradientn(colors=brewer.pal(9, 'GnBu'), guide="colorbar", na.value="white", 
                           name='', limits=0:1) + 
      coord_equal() + theme_classic() + 
      geom_sf(data=maca.mask, col = "black", fill = grey(level=0, alpha=0.2), linewidth=0.5, lty = 1) + 
      geom_sf(data=hyzo.mask, col = "black", fill = NA, linewidth=0.5, lty = 2) + 
      geom_sf(data=land, fill='NA', col=grey(level=0, alpha=1)) + 
      geom_point(data=occ.maca_nohyzo.proj, aes(x=lon.utm, y=lat.utm), fill='red', col='black', size=1, pch=21) + 
      xlim(-3250000, 3750000) + ylim(-2250000, 2750000) + 
      theme(legend.position = c(0.15, 0.4), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.line=element_blank(), axis.title=element_blank())

ggpubr::ggarrange(a, b, labels=c('A', 'B'), ncol = 2) # + bgcolor("White")
ggsave('MartesRangeWide2.eps', height=10, width=30, units='cm', dpi=1200)





# ## cropping the range-wide projections to the hybrid zone
# maam.rw.to.hyzo <- terra::crop(x=maam.rw, y=hyzo.mask)


# ggplot() + 
#       geom_spatraster(data = hill) +
#       scale_fill_gradientn(colors = pal_greys, na.value = NA) + 
#       geom_sf(data=states, color='black', linewidth = 0.4, fill=NA, lty = 1) + 
#       geom_sf(data=provinces, color='black', linewidth = 0.4, fill=NA, lty = 1) + 
#       geom_sf(data=hyzo.mask, color = "white", fill = NA, linewidth = 1, lty = 2) + 
#       geom_sf(data=i90, color='yellow',  linewidth=0.5, fill=NA, lty=1) + 
#       geom_sf(data=occ.maam_hyzo_unthin_sf, col='red') + 
#       geom_sf(data=occ.maca_hyzo_unthin_sf, col='blue') + 
#       geom_sf(data=occ.mahy_hyzo_unthin_sf, col='purple') + 
#       coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
#       theme(legend.direction = 'horizontal',
#             legend.position = 'none',
#             axis.text = element_blank(),
#             axis.ticks = element_blank(),
#             plot.title = element_text(size = 12, colour = "black", face = "bold"),
#             panel.border = element_rect(colour = "black", fill=NA)
#       ) + labs(x = NULL, y = NULL)
# ggsave('Hyzo_DEM1.png', height=20, width=20, units='cm', dpi=320)

hyzo.demplot <- ggplot() + 
      geom_spatraster(data = dem) +
      scale_fill_gradientn(colors = pal_greys, 
                           na.value = NA, name='Elevation\n(meters asl.)   ', 
                           limits=c(500,3500), breaks=seq(500,3500,1000), 
                           labels=c('≤500', '1500', '2500', '≥3500' ),
                           guide = guide_colorbar(frame.colour = "black", 
                                                  ticks.colour = "black") ) + 
      geom_sf(data=states, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=provinces, color='black', linewidth = 0.75, fill=NA, lty = 1) + 
      geom_sf(data=i90, color='black', linewidth=2.1, fill=NA, lty=1) +
      geom_sf(data=i90, color='#173E82', linewidth=2, fill=NA, lty=1) + 
      geom_sf(data=i90, color='white', linewidth=1, fill=NA, lty=1) + 
      geom_sf(data=i90, color='#A12E32', linewidth=0.6, fill=NA, lty=1) + 
      geom_sf(data=hyzo.mask, color = "black", fill = NA, linewidth = 1, lty = 2) + 
      geom_sf(data=occ.maam_hyzo_unthin_n_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      geom_sf(data=occ.mahy_hyzo_unthin_n_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      geom_sf(data=occ.maca_hyzo_unthin_n_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      geom_sf(data=occ.maca_hyzo_unthin_s_sf, pch=21, fill='dodgerblue', col='black', size=3.5 ) +
      geom_sf(data=occ.mahy_hyzo_unthin_s_sf, pch=21, fill='darkorchid1', col='black', size=3.5 ) +
      geom_sf(data=occ.maam_hyzo_unthin_s_sf, pch=21, fill='firebrick1', col='black', size=3.5 ) +
      coord_sf(xlim=ext(dem)[1:2], ylim=ext(dem)[3:4], expand=FALSE ) + 
      theme(legend.direction = 'horizontal',
            legend.position = c(0.35, 0.11), 
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
# hyzo.demplot
# ggsave('Hyzo_DEM2.png', height=20, width=20, units='cm', dpi=320)

ggarrange(hyzo.demplot, labels='A', font.label = list(color='white', size=24), label.x=0.915, label.y=0.99) #+ bgcolor("White")
ggsave('./DEM_Final1.eps', height=20, width=19, unit='cm', dpi=1200)
#




 hyzo.poly


# sphere:
library(s2)
g <- as_s2_geography(TRUE) # Earth
co <- s2_data_countries()
oc <- s2_difference(g, s2_union_agg(co)) # oceans
b <- s2_buffer_cells(as_s2_geography("POINT(-107.5 52.5)"), 9800000) # visible half
i <- s2_intersection(b, oc) # visible ocean
co <- s2_intersection(b, co)

plot(st_transform(st_as_sfc(i), "+proj=ortho +lon_0=-107.5 +lat_0=52.5"), col = 'lightblue')
plot(st_transform(st_as_sfc(co), "+proj=ortho +lon_0=-107.5 +lat_0=52.5"), col = NA, add = TRUE)
plot(st_transform(st_as_sfc(hyzo.poly), "+proj=ortho +lon_0=-107.5 +lat_0=52.5"), col='white', add=TRUE)
plot(st_transform(st_as_sfc(hyzo.poly), "+proj=ortho +lon_0=-107.5 +lat_0=52.5"), col=rgb(1,0,0,1/3), add=TRUE)

##### DEM Plots  ----------------------------------------------

### caurina
# north of i-90
dem.maca.n <- terra::extract(dem, occ.maca_hyzo_unthin_n_sf)
dem.maca.n$y  <- runif(nrow(dem.maca.n), 0.55, 0.95)
# south of i-90
dem.maca.s <- terra::extract(dem, occ.maca_hyzo_unthin_s_sf)
dem.maca.s$y  <- runif(nrow(dem.maca.s), 0.05, 0.45)

### ameriana
# north of i-90
dem.maam.n <- terra::extract(dem, occ.maam_hyzo_unthin_n_sf)
dem.maam.n$y  <- runif(nrow(dem.maam.n), 0.55, 0.95)
# south of i-90
dem.maam.s <- terra::extract(dem, occ.maam_hyzo_unthin_s_sf)
dem.maam.s$y  <- runif(nrow(dem.maam.s), 0.05, 0.45)

### hybrids
# north of i-90
dem.mahy.n <- terra::extract(dem, occ.mahy_hyzo_unthin_n_sf)
dem.mahy.n$y  <- runif(nrow(dem.mahy.n), 0.55, 0.95)
# south of i-90
dem.mahy.s <- terra::extract(dem, occ.mahy_hyzo_unthin_s_sf)
dem.mahy.s$y  <- runif(nrow(dem.mahy.s), 0.05, 0.45)

### dem plot
# transforming hyzo.poly.split that was made in MartesRFLP_1_DataPrep
hyzo.poly.split.proj <- st_transform(hyzo.poly.split, n.am.crs)
# rasterizing
hyzo.poly.split.rast <- terra::rasterize(hyzo.poly.split.proj, dem, field='ID', background=NA)
# spplot(hyzo.poly.split.rast)
# hybrid zones
# 2:3 = southern hybrid zone
# 3 = northern hybrid zone
names(hyzo.poly.split.rast) <- 'zone'
dem.rasters1 <- rast(list(dem, hyzo.poly.split.rast))
dem.rasters.df <- as.data.frame(dem.rasters1)
dem.rasters.df <- dem.rasters.df  %>% dplyr::filter(!is.na(zone))
# north
dem.n <- dem.rasters.df %>% dplyr::filter(zone == 3)
# south
dem.s <- dem.rasters.df %>% dplyr::filter(zone %in% 1:2)



demplot.maam <- ggplot() + 
      geom_hline(yintercept=0.5, linewidth=2, color='#173E82') + 
      geom_hline(yintercept=0.5, linewidth=1, color='white') + 
      geom_hline(yintercept=0.5, linewidth=0.5, color='#A12E32') + 
      geom_point(data=dem.maam.n, mapping=aes(x=NA_DEM, y=y), pch=21, fill='firebrick1', color='black', size=3.5) + 
      geom_point(data=dem.maam.s, mapping=aes(x=NA_DEM, y=y), pch=21, fill='firebrick1', color='black', size=3.5) + 
      annotate('text', x=3000, y=.75, label=paste0('n = 34'), color='firebrick1', size=5, hjust=0 ) + 
      annotate('text', x=3000, y=.25, label=paste0('n = 2'), color='firebrick1', size=5, hjust=0 ) + 
      annotate('text', x=500, y=.75, label='N', color='black', size=5, hjust=0 ) + 
      annotate('text', x=500, y=.25, label='S', color='black', size=5, hjust=0 ) + 
      scale_x_continuous(limits=c(500,3500), breaks=seq(500,3500,500), labels=seq(500,3500,500) ) + 
      scale_y_continuous(limits=0:1, name=expression(italic('M. americana')) ) + 
      theme_classic() + coord_fixed(ratio=600) + 
      theme(axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_blank())
demplot.maam

mahy.lab <- expression(paste(italic('Martes'), ' hybrids') )

demplot.mahy <- ggplot() + 
      geom_hline(yintercept=0.5, linewidth=2, color='#173E82') + 
      geom_hline(yintercept=0.5, linewidth=1, color='white') + 
      geom_hline(yintercept=0.5, linewidth=0.5, color='#A12E32') + 
      geom_point(data=dem.mahy.n, mapping=aes(x=NA_DEM, y=y), pch=21, fill='darkorchid1', color='black', size=3.5) + 
      geom_point(data=dem.mahy.s, mapping=aes(x=NA_DEM, y=y), pch=21, fill='darkorchid1', color='black', size=3.5) + 
      annotate('text', x=3000, y=.75, label=paste0('n = 19'), color='darkorchid1', size=5, hjust=0 ) + 
      annotate('text', x=3000, y=.25, label=paste0('n = 25'), color='darkorchid1', size=5, hjust=0 ) + 
      annotate('text', x=500, y=.75, label='N', color='black', size=5, hjust=0 ) + 
      annotate('text', x=500, y=.25, label='S', color='black', size=5, hjust=0 ) + 
      scale_x_continuous(limits=c(500,3500), breaks=seq(500,3500,500), labels=seq(500,3500,500) ) + 
      scale_y_continuous(limits=0:1, name=mahy.lab ) + 
      theme_classic() + coord_fixed(ratio=600) + 
      theme(axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_blank() )
demplot.mahy

demplot.maca <- ggplot() + 
      geom_hline(yintercept=0.5, linewidth=2, color='#173E82') + 
      geom_hline(yintercept=0.5, linewidth=1, color='white') + 
      geom_hline(yintercept=0.5, linewidth=0.5, color='#A12E32') + 
      geom_point(data=dem.maca.n, mapping=aes(x=NA_DEM, y=y), pch=21, fill='dodgerblue1', color='black', size=3.5) + 
      geom_point(data=dem.maca.s, mapping=aes(x=NA_DEM, y=y), pch=21, fill='dodgerblue1', color='black', size=3.5) + 
      annotate('text', x=3000, y=.75, label=paste0('n = 13'), color='dodgerblue1', size=5, hjust=0 ) + 
      annotate('text', x=3000, y=.25, label=paste0('n = 103'), color='dodgerblue1', size=5, hjust=0 ) + 
      annotate('text', x=500, y=.75, label='N', color='black', size=5, hjust=0 ) + 
      annotate('text', x=500, y=.25, label='S', color='black', size=5, hjust=0 ) + 
      scale_x_continuous(limits=c(500,3500), breaks=seq(500,3500,500), labels=seq(500,3500,500) ) + 
      scale_y_continuous(limits=0:1, name=expression(italic('M. caurina')) ) + 
      theme_classic() + coord_fixed(ratio=600) + 
      theme(axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_blank() )
demplot.maca

demplot.elev <- ggplot() + 
      geom_density(data=dem.n, aes(x=NA_DEM), lty=2) + 
      geom_density(data=dem.s, aes(x=NA_DEM)) + 
      annotate('text', x=775, y=1e-3, label='N', size=5) + 
      annotate('text', x=2600, y=5.75e-4, label='S', size=5) + 
      scale_x_continuous(limits=c(500,3500), breaks=seq(500,3500,500), 
                         labels=seq(500,3500,500), name='Elevation (meters above sea level)' ) + 
      scale_y_continuous(name='Density' ) + 
      theme_classic() + coord_fixed(ratio=420000) + 
      theme( axis.ticks.y = element_blank(), 
             axis.text.y = element_blank(), 
             axis.title.y = element_blank(),
             axis.text.x = element_text(size=16),
             axis.title.x = element_text(size=20) )
demplot.elev

demplots <- ggarrange(demplot.maam, demplot.mahy, demplot.maca, demplot.elev, font.label = list(size=24), 
                      ncol=1, labels=c('B', 'C', 'D', 'E'), label.x=0.925, label.y=0.9)# + bgcolor("White")
demplots
ggsave('./DEM_Final2.eps', plot=demplots, height=20, width=18, unit='cm', dpi=1200)





demfinal <- ggarrange(hyzo.demplot, demplots, nrow=1, labels=c('A',  ''),  widths=c(4,1))
demfinal

ggsave('./DEM_Final.eps', height=20, width=25, unit='cm', dpi=1200)
