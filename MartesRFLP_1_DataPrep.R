##### Martes Data Prep ----------------------------

### Prepping the environmental and occurrence data for running Marten ENMs

##### libraries -------------------------

library(data.table)
library(DescTools)
library(dplyr)
library(ggplot2)
library(sf)
library(stringr)
library(raster)
library(stars)
library(terra)
# library(spDataLarge)
library(beepr)
library(spThin)
library(janitor)
library(alphahull)

# Analyses ran on 2023 16" MacBook Pro with M2 Chip with 96GB RAM

##### occ data -----------------

# loading and basic cleaning of occurrence data - making sure there are no duplicate occurrences
setwd("~/Desktop/MartesRFLP")

# sourcing all the custom functions needed to run the analyses
source('./Martes_R/NAF_Functions.R')

# north american crs
n.am.crs <- '+proj=aea +lon_0=-107.5 +lat_1=37.5 +lat_2=67.5 +lat_0=52.5 +datum=WGS84 +units=m +no_defs'

### occ data

# martes americana from GBIF + Arctos
maam <- read.csv('./Martes_OccData/Maam_Final.csv')[,-(10:11)]
# martes caurina from GBIF + Arctos
maca <- read.csv('./Martes_OccData/MacaML_Final.csv')[,-(10:11)]
# genotyped martes samples from Jocie
geno <- read.csv('./Martes_OccData/martes_genotype.csv')[,-(10:11)]
# removing a genotyped caurina from southeast Alaska - ID MSB:Mamm:197437
geno <- geno %>% filter(catalogNumber != 'MSB:Mamm:197437')

# non-genotyped specimens
nogeno <- full_join(maam, maca)

# joining genotyped and non-genotyped samples
martes <- full_join(geno, nogeno, by="catalogNumber") %>% 
      mutate(species_ID=coalesce(species_ID.x,species_ID.y), 
             longitude=coalesce(longitude.x,longitude.y), 
             latitude=coalesce(latitude.x,latitude.y), 
             nuclear_ID=coalesce(nuclear_ID.x,nuclear_ID.y), 
             mito_ID=coalesce(mito_ID.x,mito_ID.y), 
             hasGenetic=coalesce(hasGenetic.x,hasGenetic.y), 
             sex=coalesce(sex.x,sex.y), 
             year=coalesce(year.x,year.y) ) %>%
      dplyr::select(c(longitude, latitude, species_ID, nuclear_ID, mito_ID, catalogNumber, hasGenetic, sex, year)) %>% 
      distinct()
martes[,1:2] <- round(martes[,1:2], 3)
martes <- unique(martes)
head(martes)

# # removing occs that are duplicates except for catalog number, sex, and year
# # makes for more straightforward spatial thinning
# martes$n <- 1:nrow(martes)
# dups <- janitor::get_dupes(martes, -c(catalogNumber, sex, year, n))$n
# martes <- martes[-dups,]
# martes$n <- NULL

# check just to make sure there are no duplicated IDs - there are 1453 unique occs
# View(table(martes$catalogNumber))
# View(martes)

##### environmental data -----------------------

setwd("~/Dropbox/Marten_ENMs/Raster Layers/wc2.0_5m_bio")
raster.names <- list.files( pattern='\\.tif$')

raster.names <- raster.names[c(4,10,12,15)]

rasters <- raster::stack(paste0('./', raster.names))
names(rasters) <- paste0('bc',c(4,10,12,15) )

rasters <- raster::crop(rasters, y=extent(-170, -50, 30, 77))

r <- rasters[[1]]

bc <- terra::rast(paste0('./', raster.names))
bc <- crop(bc, y=extent(-170, -50, 30, 77))

### shapefile of lakes and rivers
setwd("~/Desktop/MartesRFLP/Martes_Shapefiles/Rivers_and_Lakes_Shapefile/NA_Lakes_and_Rivers/data/riverslakes_p")
list.files()
# [1] "hydrography_p_lakes_v2.cpg"     "hydrography_p_lakes_v2.dbf"     "hydrography_p_lakes_v2.prj"     "hydrography_p_lakes_v2.sbn"    
# [5] "hydrography_p_lakes_v2.sbx"     "hydrography_p_lakes_v2.shp"     "hydrography_p_lakes_v2.shp.xml" "hydrography_p_lakes_v2.shx"

lakes <- st_read("hydrography_p_lakes_v2.shp")
lakes <- lakes %>% filter(grepl('Lake', NAMEEN))
lakes
lakes  <- st_transform(lakes, crs(bc))


### shapefile of continental north america
setwd("~/Desktop/MartesRFLP/Martes_Shapefiles")
na.crop <- st_read('N_Am_Crop.shp')
na.crop  <- st_transform(na.crop, crs(bc))

bc <- mask(bc, na.crop)

# defining a temporary raster that will not have lakes masked out
# this will make the range-wide plots look cleaner
temp <- bc[[1]]

bc <- mask(bc, lakes, inverse=TRUE)
names(bc) <- c('bio4', 'bio10', 'bio12', 'bio15')



mask.raster <- bc[[1]]
mask.raster[!is.na(mask.raster)] = 1
mask.raster <- terra::classify(mask.raster, cbind(-Inf, Inf,1))
plot(mask.raster)

# writing mask.raster
writeRaster(mask.raster, 'mask.raster.asc', overwrite=TRUE)

##### hybrid zone boundary ---------------

### defining the boundary of the hybrid zone based on a distance buffer from hybrid individuals


# making data frame of hybrids

mahy <- martes %>% filter(species_ID != 'americana') %>% filter(species_ID != 'caurina')
head(mahy)
mahy.pts <- vect(mahy[,1:2], geom=c('longitude', 'latitude'), crs=crs(bc))
hyhull <- convHull(mahy.pts)

# cutoff distance in meters for the boundary of the hybrid zone
cutoff.dist <- 150000
# hybrid zone extent
hy.dist <- rasterize(x=hyhull, y=mask.raster)
hy.dist <- terra::buffer(x=hy.dist, width=cutoff.dist, background=NA)
plot(hy.dist)
names(hy.dist) <- 'hyzo150'
writeRaster(hy.dist, 'hybrid.zone.mask.asc', overwrite=TRUE)

# hybrid zone rasters
hyzo.rasters <- hy.dist*bc
names(hyzo.rasters) <- names(bc)
# inverse hybrid zone raster for defining the training area for the range-wide models
hyzo.inv <- hy.dist
hyzo.inv[is.na(hyzo.inv)] <- 0
hyzo.inv[hyzo.inv==1] <- NA
hyzo.inv[!is.na(hyzo.inv)] <- 1
names(hyzo.inv) <- 'hyzo150inv'
temp <- hyzo.inv*temp
hyzo.inv <- hyzo.inv*mask.raster

### defining occurrences that are in the hybrid zone
hy.dist[is.na(hy.dist)] <- 0
names(hy.dist) <- 'in_hyzo'
in_hyzo <- terra::extract(x=hy.dist, y=martes[,1:2])
martes <- cbind(martes, in_hyzo)[,-10]; head(martes)
#    longitude latitude species_ID nuclear_ID   mito_ID   catalogNumber hasGenetic    sex year in_hyzo
# 4   -113.239   47.317  americana  americana americana MSB:Mamm:214493          1   male 2001       1
# 13  -115.010   48.007  americana  americana americana  UAM:Mamm:75657          1 female 1998       1
# 14  -113.974   48.065  americana  americana americana  UAM:Mamm:75651          1   male 1998       1
# 15  -114.596   48.083  americana  americana americana  UAM:Mamm:75607          1 female 1997       1
# 16  -113.697   48.092  americana  americana americana  UAM:Mamm:69033          1 female 1989       1
# 33  -113.893   48.437  americana  americana americana MSB:Mamm:195908          1   male 2000       1

 
table(martes$in_hyzo)
#    0    1 
# 1026  426

##### hybrid zone north and south ---------------------------------

# making a polygon of the hybrid zone boundary
hyzo.poly <- st_contour(st_as_stars(hy.dist), contour_lines = TRUE)
hyzo.poly <- st_cast(hyzo.poly, 'POLYGON')

# us interstates shapefile
# downloaded from https://hub.arcgis.com/maps/esri::usa-freeway-system/about
interstate <- read_sf('../Martes_Shapefiles/USA_Freeway_System/USA_Freeway_System.shp')
# interstate
interstate <- st_crop(interstate, st_buffer(hyzo.poly, 0.01))

# 
hyzo.poly.split <- lwgeom::st_split(hyzo.poly, interstate)
plot(hyzo.poly.split)
hyzo.poly.split <- st_collection_extract(hyzo.poly.split, 'POLYGON')
hyzo.poly.split$ID <- 1:3
plot(hyzo.poly.split, col=factor(hyzo.poly.split$ID))

# hybrid zones
# 0 = southern hybrid zone
# 1 = northern hybrid zone
hyzo.zones <- terra::rasterize(hyzo.poly.split, hy.dist, field = 'ID', background=0)
hyzo.zones[hyzo.zones == 0] <- NA
hyzo.zones[hyzo.zones < 3] <- 0
hyzo.zones[hyzo.zones > 0] <- 1
names(hyzo.zones) <- 'zone'
hyzo.rasters1 <- rast(list(hyzo.rasters, hyzo.zones))

##### area of hybrid zone + ranges -----------------------------------

# finding out how much area the hybrid zone occupies
hyzo.poly |> st_transform(n.am.crs) |> st_area()/1000^2
# 402149.3 [km^2]


#### loading in IUCN range maps for Martes
# iucn's range maps clearly include M. americana and M. caurina as a single taxonomic unit
# I made two additional polygons in QGIS that will be used to crop/mask out parts of the Martes range size for each species
# "americana_crop.shp" excludes the caurina habitat S of 45.5ºN
# "caurina.shp" includes the caurina habitat up to 51ºN

americana_crop.poly <- read_sf('MartesIUCN/americana_crop.shp') |> st_transform(n.am.crs)
caurina_tocrop_iucn.poly <- read_sf('MartesIUCN/caurina.shp') |> st_transform(n.am.crs)
# martes iucn
martes_iucn.poly <- read_sf('MartesIUCN/data_0.shp')[,c(3,16)] |> 
      st_transform(n.am.crs) |> st_simplify(dTolerance=1000)

# size of the americana range
(st_difference(martes_iucn.poly, americana_crop.poly) |> st_area()/1000^2) |> sum()
# 7121803 [km^2]

# percent size of the americana range that is within the hyrbrid zone
100*402149.3/(7121803+402149.3)
# [1] 5.344921

# size of the caurina range
st_intersection(caurina_tocrop_iucn.poly, americana_crop.poly) |> st_area()/1000^2
# 1853628 [km^2]

# percent size of the caurina range that is within the hyrbrid zone
100*402149.3/(1853628+402149.3)
# [1] 17.82753

##### range-wide sampling area ------------------------------

# each parent species
maca <- martes %>% filter(species_ID == 'caurina')
maca_hyzo <- maca %>% filter(in_hyzo == 1)
maca_nohyzo <- maca %>% filter(in_hyzo != 1)
maam <- martes %>% filter(species_ID == 'americana')
maam_hyzo <- maam %>% filter(in_hyzo == 1)
maam_nohyzo <- maam %>% filter(in_hyzo != 1)
mahy <- martes %>% filter(species_ID != 'americana') %>% filter(species_ID != 'caurina')


# making vectors of points
maca_hyzo.pts <- vect(maca_hyzo[,1:2], geom=c('longitude', 'latitude'), crs=crs(bc))
maca_nohyzo.pts <- vect(maca_nohyzo[,1:2], geom=c('longitude', 'latitude'), crs=crs(bc))

maam_hyzo.pts <- vect(maam_hyzo[,1:2], geom=c('longitude', 'latitude'), crs=crs(bc))
maam_nohyzo.pts <- vect(maam_nohyzo[,1:2], geom=c('longitude', 'latitude'), crs=crs(bc))

mahy.pts <- vect(mahy[,1:2], geom=c('longitude', 'latitude'), crs=crs(bc))

# cutoff distance M in meters
m.dist <- 250000

# distance rasters
maca.dist <- rasterize(x=maca_nohyzo.pts, y=mask.raster, field=1)
maca.dist <- terra::buffer(x=maca.dist, width=m.dist, background=NA)
maca.dist <- maca.dist*hyzo.inv; plot(maca.dist)
points(maca_nohyzo.pts, col='red', cex=0.5); points(maca_hyzo.pts, col='blue', cex=0.5); points(mahy.pts, col='purple', cex=0.5)

maam.dist <- rasterize(x=maam_nohyzo.pts, y=mask.raster, field=1)
maam.dist <- terra::buffer(x=maam.dist, width=m.dist, background=NA)
maam.dist <- maam.dist*hyzo.inv; plot(maam.dist)
points(maam_nohyzo.pts, col='red', cex=0.5);points(maam_hyzo.pts, col='blue', cex=0.5); points(mahy.pts, col='purple', cex=0.5)

writeRaster(maca.dist, 'maca.back.mask.asc', overwrite=TRUE)
writeRaster(maam.dist, 'maam.back.mask.asc', overwrite=TRUE)


# temp files for easier plotting
maca.temp <- rasterize(x=maca_nohyzo.pts, y=temp, field=1)
maca.temp <- terra::buffer(x=maca.temp, width=m.dist, background=NA)
maca.temp <- maca.temp*temp
maca.temp[!is.na(maca.temp)] <- 1; plot(maca.temp)
writeRaster(maca.temp, 'maca.temp.asc', overwrite=TRUE)

maam.temp <- rasterize(x=maam_nohyzo.pts, y=temp, field=1)
maam.temp <- terra::buffer(x=maam.temp, width=m.dist, background=NA)
maam.temp <- maam.temp*temp
maam.temp[!is.na(maam.temp)] <- 1; plot(maam.temp)
writeRaster(maam.temp, 'maam.temp.asc', overwrite=TRUE)



##### spatial thinning --------------------

setwd("~/Desktop/MartesRFLP/Martes_OccData/Martes_spThin")

# americana non hybrid zone
maam_nohyzo.thin <- thin(loc.data = maam_nohyzo, 
                  lat.col = "latitude", long.col = "longitude", 
                  spec.col = "species_ID", 
                  thin.par = 50, reps = 500, 
                  locs.thinned.list.return = TRUE, 
                  write.files = TRUE, 
                  max.files = 5, 
                  out.dir = "maam_nohyzo.thin/", out.base = "maam_nohyzo.thin", 
                  write.log.file = TRUE,
                  log.file = "maam_nohyzo.thin.txt")

# lat.long.thin.count
# 192 194 195 196 197 198 199 200 
#   1   2  15  58 146 164  95  19 
nrow(maam_nohyzo)
# [1] 806

maam_nohyzo.numbers <- as.numeric(rownames(maam_nohyzo.thin[[1]]))
maam_nohyzo.final <- maam_nohyzo[maam_nohyzo.numbers,]
length(maam_nohyzo.numbers)
# [1] 200
write.csv(maam_nohyzo.final, 'maam_nohyzo.thin.csv', row.names = F)




# americana hybrid zone
maam_hyzo.thin <- thin(loc.data = maam_hyzo, 
                       lat.col = "latitude", long.col = "longitude", 
                       spec.col = "species_ID", 
                       thin.par = 50, reps = 500, 
                       locs.thinned.list.return = TRUE, 
                       write.files = TRUE, 
                       max.files = 5, 
                       out.dir = "maam_hyzo.thin/", out.base = "maam_hyzo.thin", 
                       write.log.file = TRUE,
                       log.file = "maam_hyzo.thin.txt")
# lat.long.thin.count
#  13 
# 500
nrow(maam_hyzo)
# [1] 108

# thinning to one occurrence per grid cell in the hybrid zone
maam_hyzo1 <- thin_one_occ_per_cell(maam_hyzo, x='longitude', y='latitude', mask.raster)
# Unthinned occs: 108
# Thinned occs: 36

maam_hyzo.numbers <- as.numeric(rownames(maam_hyzo.thin[[1]]))
maam_hyzo.final <- maam_hyzo[maam_hyzo.numbers,]
# writing the thinned americana hybrid zone
write.csv(maam_hyzo.final, 'maam_hyzo.thin.csv', row.names = F)
# writing the unthinned americana hybrid zone
write.csv(maam_hyzo1, 'maam.hyzo.unthin.csv', row.names=F)


# caurina non hybrid zone
maca_nohyzo.thin <- thin(loc.data = maca_nohyzo, 
                         lat.col = "latitude", long.col = "longitude", 
                         spec.col = "species_ID", 
                         thin.par = 50, reps = 500, 
                         locs.thinned.list.return = TRUE, 
                         write.files = TRUE, 
                         max.files = 5, 
                         out.dir = "maca_nohyzo.thin/", out.base = "maca_nohyzo.thin", 
                         write.log.file = TRUE,
                         log.file = "maca_nohyzo.thin.txt")

# lat.long.thin.count
#  69  70 
# 224 276 
nrow(maca_nohyzo)
# [1] 220

maca_nohyzo.numbers <- as.numeric(rownames(maca_nohyzo.thin[[1]]))
length(maca_nohyzo.numbers)
# [1] 70
maca_nohyzo.final <- maca_nohyzo[maca_nohyzo.numbers,]
write.csv(maca_nohyzo.final, 'maca_nohyzo.final.csv', row.names = F)




# caurina hybrid zone
maca_hyzo.thin <- thin(loc.data = maca_hyzo, 
                       lat.col = "latitude", long.col = "longitude", 
                       spec.col = "species_ID", 
                       thin.par = 50, reps = 500, 
                       locs.thinned.list.return = TRUE, 
                       write.files = TRUE, 
                       max.files = 5, 
                       out.dir = "maca_hyzo.thin/", out.base = "maca_hyzo.thin", 
                       write.log.file = TRUE,
                       log.file = "maca_hyzo.thin.txt")
# lat.long.thin.count
#  27 
# 500 
nrow(maca_hyzo)
# [1] 265

# thinning to one occurrence per grid cell in the hybrid zone
maca_hyzo1 <- thin_one_occ_per_cell(maca_hyzo, x='longitude', y='latitude', mask.raster)
# Unthinned occs: 265
# Thinned occs: 116

maca_hyzo.numbers <- as.numeric(rownames(maca_hyzo.thin[[1]]))
maca_hyzo.final <- maca_hyzo[maca_hyzo.numbers,]
# writing the thinned caurina hybrid zone
write.csv(maca_hyzo.final, 'maca_hyzo.thin.csv', row.names = F)
# writing the unthinned caurina hybrid zone
write.csv(maca_hyzo1, 'maca_hyzo.unthin.csv', row.names = F)



### organizing hybrids
table(mahy$species_ID)
# AC CA XA XC 
# 11 19  9 14
mahy.all <- mahy
mahy.all$species_ID <- 'hybrid'
nrow(mahy.all)
# [1] 53

# hybrids in the hybrid zone
mahy_hyzo.thin <- thin(loc.data = mahy.all, 
                       lat.col = "latitude", long.col = "longitude", 
                       spec.col = "species_ID", 
                       thin.par = 50, reps = 500, 
                       locs.thinned.list.return = TRUE, 
                       write.files = TRUE, 
                       max.files = 5, 
                       out.dir = "maca_hyzo.thin/", out.base = "maca_hyzo.thin", 
                       write.log.file = TRUE,
                       log.file = "maca_hyzo.thin.txt")
# lat.long.thin.count
#  18  19  20 
# 131 239 130 
nrow(mahy.all)
# [1] 53

# thinning to one occurrence per grid cell in the hybrid zone
mahy_hyzo1 <- thin_one_occ_per_cell(mahy.all, x='longitude', y='latitude', mask.raster)
# Unthinned occs: 53
# Thinned occs: 44

mahy_hyzo.numbers <- as.numeric(rownames(mahy_hyzo.thin[[1]]))
mahy_hyzo.final <- maca[mahy_hyzo.numbers,]
write.csv(mahy_hyzo.final, 'mahy_hyzo.thin.csv', row.names = F)
write.csv(mahy_hyzo1, 'mahy_hyzo.unthin.csv', row.names = F)
# because we aren't actually using the hybrids for any models, we can try with the thinned and unthinned data


##### matching env data -------------------------

setwd("~/Desktop/MartesRFLP/Martes_MasterFiles")

### all north american data
# making binary raster for continuity
rast.nam <- bc
back.all <- make.background.df.terra(r.stack=rast.nam, 
                                     x.col='longitude', y.col='latitude', 
                                     time.bin='modern', region="NorthAmerica")
write.csv(back.all, 'back.all.csv', row.names = FALSE)

### hybrid zone
hy.dist[hy.dist==0] <- NA
rast.hyzo.all <- hyzo.rasters1
hyzo.all <- make.background.df.terra(r.stack=rast.hyzo.all, 
                                     x.col='longitude', y.col='latitude', 
                                     time.bin='modern', region='HybridZone')
hyzo.all$region <- ifelse(hyzo.all$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
hyzo.all$zone <- NULL
write.csv(hyzo.all, 'back.hyzo.csv', row.names = FALSE)

### americana non-hybrid zone training region
rast.maam_nohyzo.train <- bc*maam.dist
maam_nohyzo.train <- make.background.df.terra(r.stack=rast.maam_nohyzo.train, 
                                              x.col='longitude', y.col='latitude', 
                                              time.bin='modern', region='americanaTrain')
write.csv(maam_nohyzo.train, 'back.maam_nohyzo.csv', row.names = FALSE)

### caurina non-hybrid zone training region
rast.maca_nohyzo.train <- bc*maca.dist
maca_nohyzo.train <- make.background.df.terra(r.stack=rast.maca_nohyzo.train, 
                                              x.col='longitude', y.col='latitude', 
                                              time.bin='modern', region='caurinaTrain')
write.csv(maca_nohyzo.train, 'back.maca_nohyzo.csv', row.names = FALSE)


##### matching occ data -------------------------

setwd("~/Desktop/MartesRFLP/Martes_MasterFiles")

#### range-wide occurrences
### americana non-hybrid zone occs
occ.maam_nohyzo <- make.occurrence.df.terra(r.stack=rast.maam_nohyzo.train, taxa.df=maam_nohyzo.final[,c(3,1:2)],
                                            x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                            time.bin='modern', region='americanaTrain')
# Warning messages:
# Snapping 22 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.maam_nohyzo)
# [1] 200
write.csv(occ.maam_nohyzo, 'occ.maam_nohyzo.csv', row.names=FALSE)

### caurina non-hybrid zone occs
occ.maca_nohyzo <- make.occurrence.df.terra(r.stack=rast.maca_nohyzo.train, taxa.df=maca_nohyzo.final[,c(3,1:2)],
                                            x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                            time.bin='modern', region='caurinaTrain')
# Warning message:
# Snapping 1 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.maca_nohyzo)
# [1] 70
write.csv(occ.maca_nohyzo, 'occ.maca_nohyzo.csv', row.names=FALSE)

#### hybrid zone occs

### americana hybrid zone occs - thinned
occ.maam_hyzo <- make.occurrence.df.terra(r.stack=rast.hyzo.all, taxa.df=maam_hyzo.final[,c(3,1:2)],
                                          x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                          time.bin='modern', region='HybridZone')
occ.maam_hyzo$region <- ifelse(occ.maam_hyzo$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
occ.maam_hyzo$zone <- NULL
# Warning message:
# Snapping 1 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.maam_hyzo)
# [1] 13
write.csv(occ.maam_hyzo, 'occ.maam_hyzo.thin.csv', row.names=FALSE)

### caurina hybrid zone occs - thinned
occ.maca_hyzo.thin <- make.occurrence.df.terra(r.stack=rast.hyzo.all, taxa.df=maca_hyzo.final[,c(3,1:2)],
                                          x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                          time.bin='modern', region='HybridZone')
occ.maca_hyzo.thin$region <- ifelse(occ.maca_hyzo.thin$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
occ.maca_hyzo.thin$zone <- NULL
# Warning message:
# Snapping 1 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.maca_hyzo.thin)
# [1] 27
write.csv(occ.maca_hyzo.thin, 'occ.maca_hyzo.thin.csv', row.names=FALSE)

### thinned hybrids hybrid zone occs - thinned
occ.mahy_hyzo.thin <- make.occurrence.df.terra(r.stack=rast.hyzo.all, taxa.df=mahy_hyzo.final[,c(3,1:2)],
                                          x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                          time.bin='modern', region='HybridZone')
occ.mahy_hyzo.thin$region <- ifelse(occ.mahy_hyzo.thin$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
occ.mahy_hyzo.thin$zone <- NULL
# Warning messages:
# 1: Snapping 1 points to nearest grid cell based on EPSG.snap and snap.dist
# 2: 1 points were beyond snap.dist and were not matched
# Use row.names(taxa.df) to get indices of matched points
nrow(occ.mahy_hyzo.thin)
# [1] 19
write.csv(occ.mahy_hyzo.thin, 'occ.mahy_hyzo.thin.csv', row.names=FALSE)



### americana hybrid zone occs - unthinned
occ.maam_hyzo <- make.occurrence.df.terra(r.stack=rast.hyzo.all, taxa.df=maam_hyzo1[,c(3,1:2)],
                                          x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                          time.bin='modern', region='HybridZone')
occ.maam_hyzo$region <- ifelse(occ.maam_hyzo$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
occ.maam_hyzo$zone <- NULL
# Warning message:
# Snapping 1 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.maam_hyzo)
# [1] 36
write.csv(occ.maam_hyzo, 'occ.maam_hyzo.unthin.csv', row.names=FALSE)

### caurina hybrid zone occs - unthinned
occ.maca_hyzo <- make.occurrence.df.terra(r.stack=rast.hyzo.all, taxa.df=maca_hyzo1[,c(3,1:2)],
                                          x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                          time.bin='modern', region='HybridZone')
occ.maca_hyzo$region <- ifelse(occ.maca_hyzo$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
occ.maca_hyzo$zone <- NULL
# Warning message:
# Snapping 2 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.maca_hyzo)
# [1] 116
write.csv(occ.maca_hyzo, 'occ.maca_hyzo.unthin.csv', row.names=FALSE)

### thinned hybrids hybrid zone occs - unthinned
occ.mahy_hyzo.thin <- make.occurrence.df.terra(r.stack=rast.hyzo.all, taxa.df=mahy_hyzo1[,c(3,1:2)],
                                               x.col='longitude', y.col='latitude', EPSG.snap = n.am.crs, 
                                               time.bin='modern', region='HybridZone')
occ.mahy_hyzo.thin$region <- ifelse(occ.mahy_hyzo.thin$zone==0, 'HybridZoneSouth', 'HybridZoneNorth')
occ.mahy_hyzo.thin$zone <- NULL
# Warning message:
# Snapping 1 points to nearest grid cell based on EPSG.snap and snap.dist
nrow(occ.mahy_hyzo.thin)
# [1] 44
write.csv(occ.mahy_hyzo.thin, 'occ.mahy_hyzo.unthin.csv', row.names=FALSE)


##### Plotting in Bivariate E-Space ----------------------------------------

setwd("~/Dropbox/My Mac (MacBook-Pro.local)/Desktop/MartesRFLP/Martes_EnvData")

ggplot() + 
      geom_point(data=back.maam, mapping=aes(x=bio4, y=bio10), col='red', size=0.5, alpha=0.5) + 
      geom_point(data=back.maca, mapping=aes(x=bio4, y=bio10), col='blue', size=0.5, alpha=0.5) + 
      geom_point(data=back.hyzo, mapping=aes(x=bio4, y=bio10), col='black', size=0.25, alpha=0.5) + 
      xlab('Bio4 - Temp. Seasonality') + ylab('Bio10 - Mean Temp. Warmest Quarter') + xlim(0,1800) + ylim(-2,36) + 
      theme_classic()
ggsave('Bio4_Bio10.pdf', width=20, height=20, unit='cm')

ggplot() + 
      geom_point(data=back.maam, mapping=aes(x=bio4, y=bio12), col='red', size=0.5, alpha=0.5) + 
      geom_point(data=back.maca, mapping=aes(x=bio4, y=bio12), col='blue', size=0.5, alpha=0.5) + 
      geom_point(data=back.hyzo, mapping=aes(x=bio4, y=bio12), col='black', size=0.25, alpha=0.5) + 
      xlab('Bio4 - Temp. Seasonality') + ylab('Bio12 - Annual Precip.') + xlim(0,1800) + ylim(0,5400) + 
      theme_classic()
ggsave('Bio4_Bio12.pdf', width=20, height=20, unit='cm')

ggplot() + 
      geom_point(data=back.maam, mapping=aes(x=bio4, y=bio15), col='red', size=0.5, alpha=0.5) + 
      geom_point(data=back.maca, mapping=aes(x=bio4, y=bio15), col='blue', size=0.5, alpha=0.5) + 
      geom_point(data=back.hyzo, mapping=aes(x=bio4, y=bio15), col='black', size=0.25, alpha=0.5) + 
      xlab('Bio4 - Temp. Seasonality') + ylab('Bio15 - Precip. Seasonality') + xlim(0,1800) + ylim(0,110) + 
      theme_classic()
ggsave('Bio4_Bio15.pdf', width=20, height=20, unit='cm')

ggplot() + 
      geom_point(data=back.maam, mapping=aes(x=bio10, y=bio12), col='red', size=0.5, alpha=0.5) + 
      geom_point(data=back.maca, mapping=aes(x=bio10, y=bio12), col='blue', size=0.5, alpha=0.5) + 
      geom_point(data=back.hyzo, mapping=aes(x=bio10, y=bio12), col='black', size=0.25, alpha=0.5) + 
      xlab('Bio10 - Mean Temp. Warmest Quarter') + ylab('Bio12 - Annual Precip.') + xlim(-2,36) + ylim(0,5400) + 
      theme_classic()
ggsave('Bio10_Bio12.pdf', width=20, height=20, unit='cm')

ggplot() + 
      geom_point(data=back.maam, mapping=aes(x=bio10, y=bio15), col='red', size=0.5, alpha=0.5) + 
      geom_point(data=back.maca, mapping=aes(x=bio10, y=bio15), col='blue', size=0.5, alpha=0.5) + 
      geom_point(data=back.hyzo, mapping=aes(x=bio10, y=bio15), col='black', size=0.25, alpha=0.5) + 
      xlab('Bio10 - Mean Temp. Warmest Quarter') + ylab('Bio15 - Precip. Seasonality') + xlim(-2,36) + ylim(0,110) + 
      theme_classic()
ggsave('Bio10_Bio15.pdf', width=20, height=20, unit='cm')

ggplot() + 
      geom_point(data=back.maam, mapping=aes(x=bio12, y=bio15), col='red', size=0.5, alpha=0.5) + 
      geom_point(data=back.maca, mapping=aes(x=bio12, y=bio15), col='blue', size=0.5, alpha=0.5) + 
      geom_point(data=back.hyzo, mapping=aes(x=bio12, y=bio15), col='black', size=0.25, alpha=0.5) + 
      xlab('Bio12 - Annual Precip.') + ylab('Bio15 - Precip. Seasonality') + xlim(0,5400) + ylim(0,110) + 
      theme_classic()
ggsave('Bio12_Bio15.pdf', width=20, height=20, unit='cm')


summary(back.maam[,6:9]); summary(back.maca[,6:9]); summary(back.hyzo[,6:9])
#      bio4          bio10             bio12            bio15       
# Min.   :   0   Min.   :-0.2372   Min.   : 108.0   Min.   : 5.471  
# 1st Qu.:1002   1st Qu.:10.7552   1st Qu.: 338.0   1st Qu.:32.661  
# Median :1213   Median :12.9533   Median : 517.0   Median :48.384  
# Mean   :1189   Mean   :12.9409   Mean   : 661.9   Mean   :45.305  
# 3rd Qu.:1406   3rd Qu.:14.8837   3rd Qu.: 856.0   3rd Qu.:57.369  
# Max.   :1786   Max.   :23.1413   Max.   :5353.0   Max.   :96.510  
#      bio4            bio10            bio12            bio15       
# Min.   : 160.5   Min.   : 4.661   Min.   :  56.0   Min.   :  5.61  
# 1st Qu.: 709.6   1st Qu.:15.740   1st Qu.: 281.0   1st Qu.: 26.71  
# Median : 819.9   Median :18.662   Median : 386.0   Median : 44.53  
# Mean   : 793.0   Mean   :18.582   Mean   : 557.9   Mean   : 46.65  
# 3rd Qu.: 904.6   3rd Qu.:21.551   3rd Qu.: 545.0   3rd Qu.: 65.23  
# Max.   :1170.2   Max.   :36.032   Max.   :3426.0   Max.   :108.92  
#      bio4            bio10            bio12            bio15      
# Min.   : 657.4   Min.   : 5.977   Min.   : 191.0   Min.   :12.67  
# 1st Qu.: 766.0   1st Qu.:12.837   1st Qu.: 371.0   1st Qu.:25.88  
# Median : 819.1   Median :15.246   Median : 462.5   Median :37.53  
# Mean   : 831.3   Mean   :15.053   Mean   : 514.5   Mean   :39.22  
# 3rd Qu.: 879.9   3rd Qu.:17.387   3rd Qu.: 637.0   3rd Qu.:52.26  
# Max.   :1167.9   Max.   :22.608   Max.   :1230.0   Max.   :74.82 

sdata <- sdata[-zerodist(sdata )[,1],]

bivar.plot <- function(back1, back2, back3, occs1, occs2, occs3, xr, yr, l, yl, predic){
      require(alphahull)
      # standardizing all the backgrounds
      sd1x <- sd(back1[,predic[1]]); sd1y <- sd(back1[,predic[2]])
      sd2x <- sd(back2[,predic[1]]); sd2y <- sd(back2[,predic[2]])
      sd3x <- sd(back3[,predic[1]]); sd3y <- sd(back3[,predic[2]])
      back1 <- back1[,predic]; back1[,1] <- back1[,1]/sd1x; back1[,2] <- back1[,2]/sd1y
      back2 <- back2[,predic]; back2[,1] <- back2[,1]/sd2x; back2[,2] <- back2[,2]/sd2y
      back3 <- back3[,predic]; back3[,1] <- back3[,1]/sd3x; back3[,2] <- back3[,2]/sd3y
      # making alphahulls of each of the three backgrounds
      back1 <- back1[-zerodist(back1)[,1],]
      a1 <- alphahull::ahull(unique(back1), alpha=0.5)
      
}
back1 <- back.maam; back2 <- back.maca; back3 <- back.hyzo; predic=6:7
occs1  <- occ.maam_nohyzo; occs2 <- occ.maca_nohyzo; occs3 <- occ.mahy_hyzo.thin
xr <- c(0,1800); yr <- c(-2,36); xl <- 'Bio4 - Temp. Seasonality'; yl <- 'Bio10 - Mean Temp. Warmest Quarter'





