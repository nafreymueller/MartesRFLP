
##### Freymueller Maxent Pipeline -------------------------------------------------------------------------------------

# Custom functions made by Nick Freymueller for his M.S. work
# this pipeline has been generalized to work efficiently for paleo data - especially if . . .
# . . . running multiple species at a time

##### Packages -------------------------------------------------------------------------------------------------------------

# loading in the requisite packages

packs <- c('ade4', 'adehabitatMA', 'adehabitatHR', 'alphahull', 'dismo', 'dplyr', 'ecospat', 'ggplot2', 'jsonlite',
           'kuenm', 'matrixStats', 'raster', 'rgdal', 'rgeos', 'rJava',
           'sf', 'sp', 'splitstackshape', 'tidyr', 'utils', 'wesanderson', 'PerformanceAnalytics', 'SDMTools')
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

# finding where the dismo package is in your computer
# you also need to put a copy of the maxent.jar file there

# test if the maxent.jar file is in the dismo folder
list.files( system.file('java', package = 'dismo') )
# [1] "dismo.jar"  "maxent.jar"



##### Spatial Thinning --------------------------------------------------------------

### function that thins a species occurrences to only one occurrence per grid cell
## can be used in scenarios where a package like spThin may overly thin records

### ARGUMENTS
## taxa.df <- a data frame of occurrences
## rast <- a raster with the resolution you want to thin to
## x.col <- column name of the x-coordinate
## y.col <- column name of the y-coordinate
## time <- optional column name of time bins. Only needed if thinning over multiple time bins
## tax.names <- optional column name of the taxa of interest. Only needed if thinning >1 taxa at a time

thin_one_occ_per_cell <- function(taxa.df, rast, x.col='x', y.col='y', times=NULL, tax.names=NULL){
      require(terra)
      require(dplyr)
      require(DescTools)
      start.occs <- nrow(taxa.df)
      # making vectors of x and y coordinates to snap the occs to
      res.x <- res(rast)[1]
      res.y <- res(rast)[2]
      x.centers <- seq(from=as.numeric(ext(rast)[1]) + (res.x/2),
                       to=as.numeric(ext(rast)[2]) - (res.x/2),
                       by=res.x)
      y.centers <- seq(from=as.numeric(ext(rast)[3]) + (res.y/2),
                       to=as.numeric(ext(rast)[4]) - (res.y/2),
                       by=res.y)
      # making a modified data.frame of occs to be filtered
      if(is.null(times)){
            times <- rep(1, start.occs)
      } else {
            times <- taxa.df[,times]
      }
      if(is.null(tax.names)){
            tax.names <- rep(1, start.occs)
      } else {
            tax.names <- taxa.df[,tax.names]
      }
      occs1 <- data.frame(x.col=taxa.df[,x.col],
                          y.col=taxa.df[,y.col],
                          times,
                          tax.names)
      occs1$x.col <- sapply(occs1$x.col, function(x) Closest(x.centers, x) )
      occs1$y.col <- sapply(occs1$y.col, function(x) Closest(y.centers, x) )
      occs1 <- unique(occs1)
      taxa.df <- taxa.df[as.numeric(row.names(occs1)),]
      final.occs <- nrow(taxa.df)
      row.names(taxa.df) <- 1:final.occs
      writeLines(c(paste0('Unthinned occs: ', start.occs), paste0('Thinned occs: ', final.occs)) )
      return(taxa.df)
}


##### make k-folds --------------------------------------------------------------------------------------

### internal function that divides the data up into k-folds for cross validation
##  the user won't actually interact with this function personally, it's only needed for make.occurrence.df.terra

### ARGUMENTS
## x <- a data.frame of occurrences from a single taxon in a single area and time bin
## k <- number of folds for cross-validation
## seed <- a seed to be set

make.kfolds <- function(x, k=5, seed=NULL){
      if(class(x) != "data.frame"){
            class(x) <- "data.frame"
      }
      n <- nrow(x)
      rep.times <- n %/% k   # number of full reps for the rep function
      k.bins <- rep.times*k   # number of occs minus any remainder when dividing by k
      remainder <- n - k.bins   # find the remainder
      output <- rep(1:k, rep.times)   # make k bins of equal size
      if(!is.null(seed)){
            set.seed(seed)
      }
      if(remainder != 0){   # if there is a remainder, . . .
            extra <- sample(x=1:k, size=remainder)   # randomly sample it. . .
            output <- c(output, extra)   # . . . and add it to the output
      }
      output <- sample(output, size=n)   #  randomize the order of the output
      return(output)
}



##### Make Master Files ---------------------------------------------------------

### functions that make occurrence and background data.frames from a provided raster stack

### make.occurrence.df.terra extracts raster values from grid cells that have occurrences
##  this data frame will be in a format that will make it able to be merged with background data frames . . .
##  .  .  . generated with make.background.df.terra
##  the taxa must already be filtered to a single time bin and region, but multiple taxa can be ran at the same time

### ARGUMENTS
## r.stack <- a SpatRaster stack (using the terra package)
## taxa.df <- a data frame of species occurrences - must have species name, an x-coordinate, and a y-coordinate
# 
## x.col <- column name of the x-coordinate in the species data
## y.col <- column name of the y-coordinate in the species data
## time.bin  <- name of the time period the data is from - to be used for filtering later
## region <- name of the region the data is from - helpful for filtering later
## k.folds <- number of folds for cross-validation - to be passed to "make.kfolds" internally
## k.seed <- an optional seed to be sent to "make.kfolds"
## snap.dist <- maximum distance (in km) to snap points in NA-grid cells to the center of the nearest grid cell
## EPSG.in <- native crs of the input data (default is longitude/latitude in WS1984)
## EPSG.out <- crs used to project occ. points that are in NA-grid cells 

make.occurrence.df.terra <- function(r.stack, taxa.df, x.col='x', y.col='y', time.bin='time1', region=1, 
                                     k.folds=5, k.seed=NULL, snap.dist=50, EPSG.in=4326, EPSG.snap=NULL ){
      # requiring packages
      require(terra)
      require(sf)
      require(dplyr)
      # ensuring that k.folds is a positive integer between [1,26]
      if( !is.null(k.folds) ){
            if( k.folds < 1){
                  k.folds <- 1
                  print('Making occurrence data.frame without any k.folds.')
            } else if( k.folds > 26 ){
                  k.folds <- 26
                  warning('k.folds has been set to a maximum of 26.')
            } else if( k.folds%%1 != 0 ){
                  warning('k.folds has been rounded to the nearest whole number.')
                  k.folds <- round(k.folds,0)
            }
      } else {
            k.folds <- 1
            print('Making occurrence data.frame without any k.folds.')
      }
      # extracting the env-variables to a data.frame
      occs.df <- terra::extract(x=r.stack, y=taxa.df[,2:3])
      occs.df <- cbind(taxa.df, occs.df[,-1])
      colnames(occs.df)[2:3] <- c(x.col, y.col)
      # checking to see if not all points were matched
      if( nrow(occs.df) != nrow(occs.df[complete.cases(occs.df), ]) ){
            if(is.null(EPSG.snap)) stop('Points exist in NA grid cells, but EPSG.snap not provided ')
            # identifying points in need of matching
            n.snap <- nrow(occs.df) - nrow(occs.df[complete.cases(occs.df), ])
            warning(paste0('Snapping ', n.snap, ' points to nearest grid cell based on EPSG.snap and snap.dist'))
            # extracting and projecting raster cords
            r.cords <- crds(r.stack, df=TRUE)
            r.cords_sf <- st_as_sf(r.cords, coords=1:2, crs=EPSG.in)
            r.cords_sf <- st_transform(r.cords_sf, crs=EPSG.snap)
            # subsetting non-matched points and re-projecting them
            to.snap <- occs.df[!complete.cases(occs.df), 1:3]
            to.snap_sf <- st_as_sf(to.snap, coords=2:3, crs=EPSG.in)
            to.snap_sf <- st_transform(to.snap_sf, crs=EPSG.snap)
            # finding the nearest grid cell
            near.cell <- st_nearest_feature(to.snap_sf, r.cords_sf)
            r.cords_sf <- r.cords_sf[near.cell,]
            # distances
            dists <- round(as.numeric(st_distance(to.snap_sf, r.cords_sf, by_element = TRUE)/1000), 1)
            to.snap$dist <- dists
            # records to not snap and extract
            to.not.snap <- to.snap[which(to.snap$dist > snap.dist),]
            if(nrow(to.not.snap) > 0) {
                  warning(nrow(to.not.snap), ' points were beyond snap.dist and were not matched\nUse row.names(taxa.df) to get indices of matched points')
                  }
            # records to snap and extract
            to.snap <- to.snap[which(to.snap$dist <= snap.dist),]
            snap.names <- as.numeric(row.names(to.snap))
            snap.cells <- near.cell[which(to.snap$dist <= snap.dist)]
            # records already snapped and extracted
            id <- as.numeric(row.names(occs.df[complete.cases(occs.df), ]))
            occs.df <- occs.df[complete.cases(occs.df),]
            occs.df$ID <- id
            # extracting records to extract
            to.match <- to.snap
            to.match[,2:3] <- r.cords[snap.cells,]
            to.match <- terra::extract(x=r.stack, y=to.match[,2:3])[,-1]
            # assigning snapped records to to.snap
            snapped <- cbind(to.snap, to.match)
            snapped$dist <- NULL
            snapped$ID <- snap.names
            # merging snapped back with occs.df
            occs.df <- suppressMessages(full_join(occs.df, snapped))
            occs.df <- occs.df[order(occs.df$ID),]
            row.out <- occs.df$ID
            occs.df$ID <- NULL
      } else {row.out <- 1:nrow(occs.df)}# closing if-there-were-occs-in-NA-grid-cells

      # # extracting the env-variables to a data.frame
      # occs.df <- terra::extract(x=r.stack, y=taxa.df[,2:3])
      # occs.df <- cbind(taxa.df, occs.df[,-1])
      # colnames(occs.df)[2:3] <- c(x.col, y.col)
      # occs.df <- occs.df[complete.cases(occs.df), ]
      colnames(occs.df)[1] <-'Name'
      n.occs <- nrow(occs.df)
      # making the data.frame of the name, time.bin, region, and presence columns
      d.cols <- as.data.frame( matrix(0, nrow=n.occs, ncol=3) )
      d.cols[,1] <- time.bin
      d.cols[,2] <- region
      colnames(d.cols) <- c('time.bin', 'region', 'presence')
      # merging everything together
      occs.df <- do.call('cbind', list(occs.df[,1:3], d.cols[,1:2], occs.df[,4:ncol(occs.df)], d.cols[,3] ) )
      colnames(occs.df)[ ncol(occs.df) ] <- 'presence'
      
      
      # making the cross-validation columns
      if( k.folds > 1 ){
            ltrs <- letters[1:k.folds]
            # training columns
            train <- as.data.frame( matrix(1, nrow=n.occs, ncol=k.folds) )
            for(i in 1:k.folds){
                  colnames(train)[i] <- paste0('tr.', paste(ltrs[-(k.folds+1-i)], collapse='') )
            }
            # testing columns
            test <- as.data.frame( matrix(0, nrow=n.occs, ncol=k.folds) )
            for(i in 1:k.folds){
                  colnames(test)[i] <- paste0('test.', paste(ltrs[(k.folds+1-i)], collapse='') )
            }
            # splitting up the occurrence data.frame by taxon, then assigning the k.folds
            occs.list <- base::split(occs.df, occs.df[,1])
            n.taxa <- length(occs.list)
            for(i in 1:n.taxa){
                  
                  # setting all taxa with < k.folds occurrences to 0's. they will later be converted to exist in all the training and testing subsets
                  if( nrow( occs.list[[i]] ) < k.folds ){
                        occs.list[[i]]$presence <- 0
                  } else {
                        # filling in cross-validation folds for taxa who have > k.folds in terms of occurrences
                        occs.list[[i]]$presence <- make.kfolds(occs.list[[i]], k=k.folds, seed=k.seed)
                  }
            } # closing for(i in 1:n.taxa)
            # merging all the different taxa back together
            occs.df <- do.call('rbind', occs.list)
            colnames(occs.df)[ ncol(occs.df) ] <- 'presence'
            row.names(occs.df) <- row.out
            # assigning the k.fold parameters in occs.df
            for(i in 1:nrow(occs.df) ){
                  # for taxa with fewer occs than k.folds
                  if( occs.df$presence[i] == 0 ){
                        # train[i,] <- 1   # train is filled with 1's by default
                        test[i,] <- 1
                  } else {
                        # for taxa with greater occs than k.folds
                        cv <- occs.df$presence[i]
                        train[i, (k.folds+1-cv) ] <- 0
                        test[i, (k.folds+1-cv) ] <- 1
                  }
            } # closing for(i in 1:nrow(occs.df) )
            out <- do.call('cbind', list( occs.df, train, test ) )
      } # closing if( k.folds > 1 )
      # setting all the presence rows to  1
      out$presence <- 1
      # output
      return(out)
}


### make.background.df.terra extracts all non-NA grid cells into a data frame
##  this data frame will be in a format that will make it able to be merged with occurrence data frames . . .
##  .  .  . generated with make.occurrence.df.terra

### ARGUMENTS
## r.stack <- a SpatRaster stack (using the terra package)
## taxa.df <- a data frame of species occurrences - must have species name, an x-coordinate, and a y-coordinate
## x.col <- column name of the x-coordinate in the species data - choose the same name as x.col in make.occurrence.df.terra
## y.col <- column name of the y-coordinate in the species data - choose the same name as y.col in make.occurrence.df.terra
## time.bin  <- name of the time period the data is from - to be used for filtering later
## region <- name of the region the data is from - helpful for filtering later
## k.folds <- number of folds for cross-validation - to be passed to "make.kfolds" internally
## k.seed <- an optional seed to be sent to "make.kfolds"

make.background.df.terra <- function(r.stack, x.col='x', y.col='y', 
                                     time.bin='time1', region=1, k.folds=5){
      # ensuring that k.folds is a positive integer between [1,26]
      if( !is.null(k.folds) ){
            if( k.folds < 1){
                  k.folds <- 1
                  print('NOTE: Making background data.frame without any k.folds.')
            } else if( k.folds > 26 ){
                  k.folds <- 26
                  warning('k.folds has been set to 26.')
            } else if( k.folds%%1 != 0 ){
                  warning('k.folds has been rounded to the nearest whole number.')
                  k.folds <- round(k.folds)
            }
      } else {
            k.folds <- 1
            print('NOTE: Making background data.frame without any k.folds.')
      }
      # extracting the env-variables to a data.frame
      r <- xyFromCell(r.stack[[1]], cell=1:ncell(r.stack[[1]]))
      coords.df <- terra::extract(r.stack, r)
      coords.df <- cbind(r, coords.df)
      colnames(coords.df)[1:2] <- c(x.col, y.col)
      coords.df <- coords.df[complete.cases(coords.df), ]
      n.points <- nrow(coords.df)
      # making the data.frame of the name, time.bin, region, and presence columns
      d.cols <- as.data.frame( matrix(0, nrow=n.points, ncol=4) )
      d.cols[,1] <- 'Background'
      d.cols[,2] <- time.bin
      d.cols[,3] <- region
      colnames(d.cols) <- c('Name', 'time.bin', 'region', 'presence')
      # merging everything together
      out <- do.call('cbind', list( d.cols[,1],  coords.df[,1:2], d.cols[,2:3], coords.df[,3:ncol(coords.df)], d.cols[,4] ))
      colnames(out)[1] <- 'Name'
      colnames(out)[ ncol(out) ] <- 'presence'
      # making the cross-validation columns
      if( k.folds > 1 ){
            ltrs <- letters[1:k.folds]
            # training columns
            train <- as.data.frame( matrix(0, nrow=n.points, ncol=k.folds) )
            for(i in 1:k.folds){
                  colnames(train)[i] <- paste0('tr.', paste(ltrs[-(k.folds+1-i)], collapse='') )
            }
            # testing columns
            test <- as.data.frame( matrix(0, nrow=n.points, ncol=k.folds) )
            for(i in 1:k.folds){
                  colnames(test)[i] <- paste0('test.', paste(ltrs[(k.folds+1-i)], collapse='') )
            }
            out <- do.call('cbind', list( out, train, test ) )
      }
      # output
      return(out)
}



##### Prep Parameters ------------------------------------------------------------------------------------------------------

## internal function that preps parameters for Maxent models in R with the dismo package
# the user won't interact with this function personally. it is used in optimize.maxent.likelihood and maxent.crossval.error

# this function was made by Xiao Feng - all credit for this function is his alone
# sourced from https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix2_prepPara.R
#  https://github.com/shandongfx/workshop_maxent_R/blob/master/code

# should be used in concert with Appendix 3 from Feng et al. 2017 (PeerJ Preprints)
# https://peerj.com/preprints/3346.pdf (manuscript still unpublished as of August 2020)
# Appendix 3
# https://github.com/shandongfx/workshop_maxent_R/blob/master/code/Appendix3_maxentParameters_v2.pdf

# some arguments may change, or may/may not be used depending on if you're using raster data vs "samples-with-data" (SWD; column data) 

# A function that implements Maxent parameters using the general R manner
# leave "doclamp" as default - later in the code (internal fxns run), "doclamp" is set to FALSE
prepPara <- function(userfeatures=NULL,             # 41  NULL=autofeature, could be any combination of # c("L", "Q", "H", "P", "T")
                     #     MUST be specified as a single string (e.g., "LQ", "LQP", "LQHPT", etc.)
                     responsecurves=TRUE,           # 1
                     jackknife=TRUE,                # 3
                     outputformat="logistic",       # 4
                     outputfiletype="asc",          # 5
                     projectionlayers=NULL,         # 7
                     randomseed=FALSE,              # 10
                     removeduplicates=TRUE,         # 16
                     betamultiplier=NULL,           # 20, 53-56
                     biasfile=NULL,                 # 22
                     testsamplesfile=NULL,          # 23
                     replicates=1,                  # 24-25
                     replicatetype="crossvalidate", # 24-25
                     writeplotdata=TRUE,            # 37
                     extrapolate=TRUE,              # 39
                     doclamp=TRUE,                  # 42
                     beta_threshold=NULL,           # 20, 53-56
                     beta_categorical=NULL,         # 20, 53-56
                     beta_lqp=NULL,                 # 20, 53-56
                     beta_hinge=NULL,               # 20, 53-56
                     applythresholdrule=NULL        # 60
){
  #20, 29-33, & 41 features, default is autofeature
  if(is.null(userfeatures)){
    args_out <- c("autofeature")
  } else {
    args_out <- c("noautofeature")
    if(grepl("L",userfeatures)) args_out <- c(args_out,"linear") else args_out <- c(args_out,"nolinear")
    if(grepl("Q",userfeatures)) args_out <- c(args_out,"quadratic") else args_out <- c(args_out,"noquadratic")
    if(grepl("H",userfeatures)) args_out <- c(args_out,"hinge") else args_out <- c(args_out,"nohinge")
    if(grepl("P",userfeatures)) args_out <- c(args_out,"product") else args_out <- c(args_out,"noproduct")
    if(grepl("T",userfeatures)) args_out <- c(args_out,"threshold") else args_out <- c(args_out,"nothreshold")
  }
  
  # 1 - generate response curves for each variable
  if(responsecurves) args_out <- c(args_out,"responsecurves") else args_out <- c(args_out,"noresponsecurves")
  
  # 2
  #if(picture) args_out <- c(args_out,"pictures") else args_out <- c(args_out,"nopictures")
  
  # 3 - apply variable jackknife to see how the model changes if that variable is omitted, then if it's the ONLY variable used
  if(jackknife) args_out <- c(args_out,"jackknife") else args_out <- c(args_out,"nojackknife")
  
  # 4 - output format type. choose from c("logistic", "cumulative", "raw")
  args_out <- c(args_out,paste0("outputformat=",outputformat))
  
  # 5 - output file type. choose from c("asc", "mxe", "grd", "bil")
  args_out <- c(args_out,paste0("outputfiletype=",outputfiletype))
  
  # 7 - pathway to projection layers.
  # it seems that the projection layers should be the only files in that folder, just like with the MaxEnt .jar file
  if(!is.null(projectionlayers))    args_out <- c(args_out,paste0("projectionlayers=",projectionlayers))
  
  # 10 - will use different random number generators for selecting training vs testing data and background points (if applicable)
  if(randomseed) args_out <- c(args_out,"randomseed") else args_out <- c(args_out,"norandomseed")
  
  # 16 - remove duplicate coordinates that are in the same grid  - ONLY for raster data, not SWD
  if(removeduplicates) args_out <- c(args_out,"removeduplicates") else args_out <- c(args_out,"noremoveduplicates")
  
  # 20 & 53-56 - various beta (regularization) multipliers to be applied. default = 1
  # 20 applies all parameters by this regularization multiplier.
  # 53-56 can apply uniquely to different feature types
  # check if negative
  betas <- c( betamultiplier,beta_threshold,beta_categorical,beta_lqp,beta_hinge)
  if(! is.null(betas) ){
    for(i in 1:length(betas)){
      if(betas[i] <0) stop("betamultiplier has to be positive")
    }
  }
  if (  !is.null(betamultiplier)  ){
    args_out <- c(args_out,paste0("betamultiplier=",betamultiplier))
  } else {
    if(!is.null(beta_threshold)) args_out <- c(args_out,paste0("beta_threshold=",beta_threshold))
    if(!is.null(beta_categorical)) args_out <- c(args_out,paste0("beta_categorical=",beta_categorical))
    if(!is.null(beta_lqp)) args_out <- c(args_out,paste0("beta_lqp=",beta_lqp))
    if(!is.null(beta_hinge)) args_out <- c(args_out,paste0("beta_hinge=",beta_hinge))
  }
  
  # 22 - pathway to a bias file for selecting background points - ONLY for raster data, not SWD
  if(!is.null(biasfile))    args_out <- c(args_out,paste0("biasfile=",biasfile))
  
  # 23 - pathway to a test data file - can be in csv format
  if(!is.null(testsamplesfile))    args_out <- c(args_out,paste0("testsamplesfile=",testsamplesfile))
  
  # 24 - replicates = number of replicates to run (integer)
  # 25 - replicatetype = what type of replicates to run. choose from c('crossvalidate', 'bootstrap', 'subsample')
  replicates <- as.integer(replicates)
  if(replicates>1 ){
    args_out <- c(args_out,
                  paste0("replicates=",replicates),
                  paste0("replicatetype=",replicatetype) )
  }
  
  # 37 - write output files containing the data used to make response curves
  if(writeplotdata) args_out <- c(args_out,"writeplotdata") else args_out <- c(args_out,"nowriteplotdata")
  
  # 39 - allow extrapolation beyond the limits of the training data
  if(extrapolate) args_out <- c(args_out,"extrapolate") else args_out <- c(args_out,"noextrapolate")
  
  # 42 - apply clamping when projecting
  if(doclamp) args_out <- c(args_out,"doclamp") else args_out <- c(args_out,"nodoclamp")
  
  # 60 - threshold your model to binary 1/0
  # options are: c('Fixed cumulative value 1', 'Fixed cumulative value 5', 'Fixed cumulative value 10', 'Minimum training presence',
  # '10 percentile training presence', 'Equal training sensitivity and specificity', 'Maximum training sensitivity plus specificity').
  if(!is.null(applythresholdrule))    args_out <- c(args_out,paste0("applythresholdrule=",applythresholdrule))
  
  return(args_out)
}

# prepPara()
# [1] "autofeature"           "responsecurves"        "jackknife"             "outputformat=logistic" "outputfiletype=asc"   
# [6] "norandomseed"          "removeduplicates"      "writeplotdata"         "extrapolate"           "doclamp"



##### Create Null Object, Summary Object, Eval Object, and Nested File Structure -------------------------------------------

### functions that create the null data.frame, summary data.frame, and the eval data.frame

### ARGUMENTS 1 - arguments that exist to keep track of everything and do not change how the functions run
## taxon.name <- name of the entity that will get assigned to any null/summary/eval objects
## time.bin <-   name of the time bin that will get assigned to any null/summary/eval objects
## extent <-     name of the extent that will get assigned to any null/summary/eval objects

### ARGUMENTS 2 - arguments that set the model hyper-parameters and will change how the functions run
## cv.runs <-    name(s) of the cross-validation folds. folders will be generated with these names in create.folders.for.maxent.
#                recommended framework is to treat each cross-validation fold as a letter.
#                In the case of 5-fold cross validation, for running all the data (using the null.aic and optimize.maxent.likelihood . . .
#                . . . functions), you should specify 'abcde'. For running cross validation models with maxent.crossval.error, you . . .
#                . . . should specify c('abcd', 'abce', 'abde', 'acde', 'bcde').
## f.class <-    name(s) of the feature classes used in analysis.
## beta.values <- regularization multipliers used in analysis

# create.summary.df and create.eval.df will create every possible combination of cv.runs, f.class, and beta.values

# function to create the null data.frame
create.null.df <- function(taxon.name, time.bin, extent){
  null.object <- as.data.frame(matrix(NA, nrow=1, ncol=17))
  colnames(null.object) <- c( "Taxa", "Time.Bin", "Extent", "CrossVal", "cv.num", "Features", "Betas", "n", "k",
                              "ln.L", "AIC", "AICc", "delta.i", "delta.i.c", "w.i", "w.i.c", "lambdas" )
  null.object$Taxa <- taxon.name
  null.object$Time.Bin <- time.bin
  null.object$Extent <- extent
  null.object$CrossVal <- 'abcde'
  null.object$cv.num <- 1
  null.object$lambdas <- 'Incercept'
  
  return(null.object)
}


# function to create the summary data.frame - to be used as an input for optimize.maxent.likelihood
create.summary.df <- function(taxon.name,
                              time.bin,
                              extent,
                              cv.runs = 'abcde',
                              f.class = c('LQP', 'Q'),
                              beta.values  = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5) ){
                                
  summary.object <- expand.grid(CrossVal=cv.runs, Features = f.class, Betas = beta.values, stringsAsFactors = T)
  summary.object$cv.num <- as.numeric(summary.object$CrossVal)
  
  # ifelse functions defining what to do if taxon.name, time.bin, and extent are not specified
  if( !is.null(taxon.name) ){
    summary.object$Taxa <- taxon.name
  } else {
    summary.object$Taxa <- 'Taxon1'
  }
  if( !is.null(time.bin) ){
    summary.object$Time.Bin <- time.bin
  } else {
    summary.object$Time.Bin <- 'TimeBin1'
  }
  if( !is.null(extent) ){
    summary.object$Extent <- extent
  } else {
    summary.object$Extent <- 'Extent1'
  }
  
  # re-ordering columns
  summary.object <- summary.object[,c(5:7, 1, 4, 2:3)]
  # adding in all the other parameters
  summary.object$n <- NA   # sample size
  summary.object$k <- NA   # number of non-zero lambdas
  summary.object$ln.L <- NA   # log-likelihood
  summary.object$AIC <- NA   # AIC
  summary.object$AICc <- NA   # AICc corrected for small sample size
  summary.object$delta.i  <- NA   # delta.i for AIC - wont calculate everything until all models have been run
  summary.object$delta.i.c  <- NA   # delta.i for AICc - wont calculate everything until all models have been run
  summary.object$w.i  <- NA   # Akaike Weights for AIC - wont calculate everything until all models have been run
  summary.object$w.i.c  <- NA   # Akaike Weights for AIC - wont calculate everything until all models have been run
  summary.object$lambdas  <- NA   # list of all the non-zero lambdas
  summary.object[,4] <- as.character(summary.object[,4])
  summary.object[,6] <- as.character(summary.object[,6])
  
  return(summary.object)
}


# function to create the eval data.frame - to be used as an input for maxent.crossval.error
create.eval.df <- function(taxon.name,
                           time.bin,
                           extent,
                           cv.runs = c('abcd', 'abce', 'abde', 'acde', 'bcde'),
                           f.class = 'LQP',
                           beta.values = 1 ){
  
  summary.object <- expand.grid(CrossVal=cv.runs, Features = f.class, Betas = beta.values, stringsAsFactors = T)
  summary.object$cv.num <- 1:length(cv.runs)
  
  # ifelse functions defining what to do if taxon.name, time.bin, and extent are not specified
  if( !is.null(taxon.name) ){
    summary.object$Taxa <- taxon.name
  } else {
    summary.object$Taxa <- 'Taxon1'
  }
  if( !is.null(time.bin) ){
    summary.object$Time.Bin <- time.bin
  } else {
    summary.object$Time.Bin <- 'TimeBin1'
  }
  if( !is.null(extent) ){
    summary.object$Extent <- extent
  } else {
    summary.object$Extent <- 'Extent1'
  }
  
  # re-ordering columns
  summary.object <- summary.object[,c(5:7, 1, 4, 2:3)]
  # adding in all the other parameters
  summary.object$n <- NA   # sample size
  summary.object$k <- NA   # number of non-zero lambdas
  summary.object$ln.L <- NA   # log-likelihood
  summary.object$AIC <- NA   # AIC
  summary.object$AICc <- NA   # AICc corrected for small sample size
  summary.object$delta.i  <- NA   # delta.i for AIC - wont calculate everything until all models have been run
  summary.object$delta.i.c  <- NA   # delta.i for AICc - wont calculate everything until all models have been run
  summary.object$w.i  <- NA   # Akaike Weights for AIC - wont calculate everything until all models have been run
  summary.object$w.i.c  <- NA   # Akaike Weights for AIC - wont calculate everything until all models have been run
  summary.object$lambdas  <- NA   # list of all the non-zero lambdas
  summary.object[,4] <- as.character(summary.object[,4])
  summary.object[,6] <- as.character(summary.object[,6])
  
  # renaming the columns
  colnames(summary.object) <- c('Taxa', 'Time.Bin', 'Extent', 'CrossVal', 'cv.num', 'Features', 'Betas', 'n', 'thresh', 'test.sens',
                                'pROC_0.1', 'pval_0.1', 'pROC_1', 'pval_1', 'pROC_5', 'pval_5', 'lambdas')
  
  return(summary.object)
}


# function to create the nested fie structure that it will use to store the output (including html files) of maxent models 
# summary.eval.df <- the summary.df data.frame or the eval.df data.frame.
# the function will use information from the cv.runs, f.class, and beta.values columns to create this nested file structure
# wd <- working directory that the nested file structure is going to be generated in. if running multiple species, it is recommended . . .
# . . . that you make a folder for each species, then run this function in each species folder
create.folders.for.maxent <- function(summary.eval.df, wd = getwd() ){
  
  # setting the working directory
  setwd( getwd() )
  
  # for loop that goes through the summary/evaluation data.frame and creates the nested file structure
  for( i in 1:nrow(summary.eval.df) ){
    
    if( !dir.exists( paste(getwd(), summary.eval.df$CrossVal[i], sep='/' ) ) ){   # if the cross-validation folder exists
      writeLines( c('Creating folder:', paste(getwd(), summary.eval.df$CrossVal[i], sep='/'), sep='') )
      dir.create( paste(getwd(), summary.eval.df$CrossVal[i], sep='/') )
      setwd( paste(getwd(), summary.eval.df$CrossVal[i], sep='/') )
    } else {
      setwd( paste(getwd(), summary.eval.df$CrossVal[i], sep='/') )
    }
    
    if( !dir.exists( paste(getwd(), summary.eval.df$Features[i], sep='/' ) ) ){   # if the feature class folder exists
      writeLines( c('Creating folder:', paste(getwd(), summary.eval.df$Features[i], sep='/'), sep='') )
      dir.create( paste(getwd(), summary.eval.df$Features[i], sep='/') )
      setwd( paste(getwd(), summary.eval.df$Features[i], sep='/') )
    } else {
      setwd( paste(getwd(), summary.eval.df$Features[i], sep='/') )
    }
    
    if( !dir.exists( paste(getwd(), summary.eval.df$Features[i], sep='/' ) ) ){   # if the regularization folder exists
      writeLines( c('Creating folder:', paste(getwd(), paste('beta',  summary.eval.df$Betas[i], sep='_'),  sep='/'), sep='') )
      dir.create( paste(getwd(), paste('beta',  summary.eval.df$Betas[i], sep='_'), sep='/') )
      setwd( paste(getwd(), paste('beta',  summary.eval.df$Betas[i], sep='_'), sep='/') )
    } else {
      setwd( paste(getwd(), paste('beta',  summary.eval.df$Betas[i], sep='_'), sep='/') )
    }
    
    setwd('../../../')   # go three directories back
    
    }   # finishing main for loop
  
  

  
}


# setwd("~/Dropbox/SVP_Models/ModelOutput/Tyrano")
# 
# 
# examp <- create.summary.df('Tyrano', 'K', 'Laurimidia')
# create.folders.for.maxent(examp)
# 
# examp2 <- create.eval.df('Tyrano', 'K', 'Laurimidia')
# create.folders.for.maxent(examp2)




##### Null AICc ------------------------------------------------------------------------------------------------------------

## function that calculates AIC and AICc values for a null (i.e., intercept-only) model
#  to be compared with the "best" model generated via optimize.maxent.likelihood

### ARGUMENTS
# the arguments are the same as the optimize.maxent.likelihood model, except that null.df should be a data.frame created from the . . .
# . . . create.null.df object
# function will return 

null.aic <- function(null.df, occs, background, first.occ.col){
      
      # number of parameters
      null.df$k <- 1
      
      # assigning the column number to be sent through maxent
      col.number <- first.occ.col
      # filtering the species dataset by col.number and assigning to null.df
      sp <- occs %>% dplyr::filter(occs[,col.number] == 1)
      n <- nrow(sp)
      null.df[1,8] <- n
      
      # giving noting the variables/features/hyperparameters with non-zero lambdas
      null.df[1,17] <- 'Intercept'
      
      out.scale <- rep(1, times=n ) / nrow(background)
      
      log.like <- sum(log(out.scale))
      
      null.df[1,10] <- log.like
      
      ### calculating AIC and AICc
      # AIC
      AIC <- 2 - 2*log.like
      null.df[1,11] <- AIC
      # AICc for one term null model
      null.df[1,12] <- AIC + 4/(n - 2)
      
      return(null.df)
}



##### Optimize Maxent Likelihood -------------------------------------------------------------------------------------------

## function that runs a series of maxent models with varying parameter settings such that the user can . . .
#  .  .  . find one with  optimum settinggsoptimize their
# function will return a filled out summary model

optimize.maxent.likelihood <- function(summary.df,         # summary object that will keep all the model output
                                                           # the nested file structure created from create.folders.for.maxent MUST exist
                                       occs,               # species occurrence object
                                       background,         # sampled Background object (COLUMNS MUST BE IDENTICAL TO occs)
                                       predic,             # column numbers of the predictor variables
                                       first.occ.col,      # number of the first cross-validation column in occs/background
                                       home=getwd(),       # directory where all the models will be ran. 
                                       all.models = TRUE   # do you want to keep all versions of all the models?
                                       # helpful for quickly checking some models, but may consume loads (e.g., >1GB) . . .
                                       # . . . of hard disk space. Recommended to set to FALSE for exploratory analyses.
                                       # if FALSE, function will create a folder called "RunOver" and will . . .
                                       # . . . continuously write-over it for all models, and the only model you see . . .
                                       # . . . at the end will be the last model that was ran
){
  # Calculating the total number of models
  nmodels <- nrow(summary.df)
  
  # prompting the user if they want to store models in the RAM
  print(paste0('The time is ', Sys.time(), '. You are running ', nmodels, ' total MaxEnt models.'))
  
  # setting up the progress bar
  prog <- txtProgressBar(min=0, max=nrow(summary.df), style=3,  char='+')
  
  for(i in 1:nrow(summary.df)){
    
    # setting the working directory for each folder
    if(all.models == TRUE){
      setwd( paste(home, summary.df[i,4], summary.df[i,6], paste('beta', summary.df[i,7], sep='_'), sep='/') )
    } else {
      # create the RunOver folder
      if( !dir.exists( paste(home, 'RunOver', sep='/') ) ){
        dir.create( paste(home, 'RunOver', sep='/') )
        setwd( paste(home, 'RunOver', sep='/') )
      } else {
        setwd( paste(home, 'RunOver', sep='/') )
      }
    }
    
    ###########################################################################
    
    ### MaxEnt things happen here
    
    ### preparing the data for the maxent model
    # filtering the occ object by it's respective cross-validation identity
    cv.number <- summary.df$cv.num[i]
    # assigning the column number to be sent through maxent
    col.number <- first.occ.col + cv.number - 1
    # filtering the species dataset by col.number and assigning to summary.df
    sp <- occs %>% dplyr::filter(occs[,col.number] == 1)
    n <- nrow(sp)
    summary.df[i,8] <- n
    # bending the occ and background data together
    mx.data <- rbind(sp, background)
    
    ### running the actual maxent model
    mx.model <- dismo::maxent(p = mx.data[,col.number],
                              x = mx.data[,predic],
                              path =  paste0(getwd()),
                              args = prepPara(userfeatures = summary.df[i,6], betamultiplier = summary.df[i,7], doclamp = FALSE)
    )
    
    ### calculating k and the names of the lambdas
    lambda.file <- as.data.frame(mx.model@lambdas) %>% `colnames<-`('lambdas')
    # lambdas data frame
    lambdas.df <- splitstackshape::cSplit(lambda.file, sep=',', splitCols='lambdas')
    colnames(lambdas.df) <- c('feature', 'lambda', 'min', 'max')
    # finding the non-zero lambdas
    non.zero.lambdas <- lambdas.df %>% dplyr::filter(!is.na(max)) %>% dplyr::filter(lambda != 0)
    class(non.zero.lambdas) <- 'data.frame'
    # assigning the number of parameters
    k <- nrow(non.zero.lambdas)
    # if beta is too high, the model gets over-regularized to the point that all the lambda coefficients get set to 0
    #  this effectively becomes an intercept only model, which is effectively the global mean
    # in this sense, k should get set to 0
    if(k == 0){
      k <- 1
    }
    summary.df[i,9] <- k
    # giving noting the variables/features/hyperparameters with non-zero lambdas
    summary.df[i,17] <- toString(non.zero.lambdas[,1])
    
    ### calculating the log-likelihood, then AIC and AICc
    # if statement calculating if there is an appropriate AIC value
    # e.g., can't fit 4 observations (occurrence points) with 5 variables
    if(n - k < 2){
      summary.df[i,10] <- NA
      summary.df[i,11] <- NA
      summary.df[i,12] <- NA
    } else {   # if  it is possible to calculate AIC and AICc
      # logistic model output
      mx.back <- dismo::predict(mx.model, background[,predic])
      # sum of all background point values - mx.back / back.sum should = 1
      back.sum <- sum(mx.back)
      # logistic values of the (k-1)/k occurrences
      mx.occs <- dismo::predict(mx.model, sp[,predic])
      # scaling to make compatible for calculating AIC
      occs.raw <- mx.occs / back.sum
      # log(likelihood)
      log.like  <- sum(log(occs.raw))
      summary.df[i,10]<- log.like
      
      ### calculating AIC and AICc
      # AIC
      AIC <- 2*k - 2*log.like
      summary.df[i,11] <- AIC
      # AICc
      summary.df[i,12] <- AIC + 2*((k^2 + k) / (n - k - 1))
      
    }
    # updating the prograss bar for each run to get an idea of how long things will take
    setTxtProgressBar(prog, i)
    
    
    ###########################################################################
    # returning to the home directory
    setwd(home)
    
  }   # closes the for loop
  
  return(summary.df)
}





##### MaxEnt Cross-Validation Error -----------------------------------------------------------------------------------------------

# function that runs five-fold cross-validation once you've found the optimum hyper-parameters for your maxent model(s)
# functions returns a list containing 3 objects:
# 1.) a filled out eval object
# 2.) an object containing all the maxent models; and . . .
# 3.) a data.frame with dimensions [ 1:nrow(background), 1:nrow(eval.df) ] containing the projections of all the models in

maxent.crossval.error <- function(eval.df,            # object generated from create.eval.df that will keep all the model output
                                  occs,               # species occurrence object
                                  background,         # sampled Background object (columns MUST be identical to columns in occs)
                                  predic,             # column numbers of the predictor variables
                                  first.occ.col,      # number of the first training in occs/background (assumes k = 5)
                                  # NOTE: this is not the presence column! it's the column to the right of it
                                  first.test.col,     # number of the first testing column in occs/background (assumes k = 5)
                                  home=getwd(),       # directory where all the models will be ran
                                  all.models = TRUE,  # do you want to keep all versions of all the models?
                                  # setting to TRUE is helpful for quickly checking some models, but may consume >1GB . . .
                                  # . . . of hard disk space. Recommended to set to FALSE for exploratory analyses.
                                  # if FALSE, function will create a folder called "RunOver" and will . . .
                                  # . . . continuously write-over it for all models, and the only model you see . . .
                                  # . . . at the end will be the last model that was ran
                                  omission.rate=0,    # user specified omission rate to be calculated for the threshold
                                  # express as a proportion from 0-1
                                  all.background,     # background data from the extents the species exists in
                                  # if using all the potential background points, all.background should == background
                                  pROC.reps = 500     # number of iterations each partialROC test will go through
){
  # Calculating the total number of models
  nmodels <- nrow(eval.df)
  
  # stop if an incompatible omission rate is specified
  if(omission.rate > 1  || omission.rate < 0){
    stop('Specify an omission rate between 0-1.')
  }
  
  # prompting the user if they want to store models in the RAM
  print(paste0('The time is ', Sys.time(), '. You are running ', nmodels, ' total MaxEnt models.'))
  
  # setting up the progress bar
  prog <- txtProgressBar(min=0, max=nrow(eval.df), style=3,  char='+')
  
  # setting up various list objects to send output  to
  # list for the maxent models
  maxent.list <- list()
  # list for the testing background data
  test.back.list <- list()
  # list for all the occs
  all.occs.list <- list()
  
  for(i in 1:nrow(eval.df)){
    
    # setting the working directory for each folder
    if(all.models == TRUE){
      setwd( paste(home, eval.df[i,4], eval.df[i,6], paste('beta', eval.df[i,7], sep='_'), sep='/') )
    } else {
      # create the RunOver folder
      if( !dir.exists( paste(home, 'RunOver', sep='/') ) ){
        dir.create( paste(home, 'RunOver', sep='/') )
        setwd( paste(home, 'RunOver', sep='/') )
      } else {
        setwd( paste(home, 'RunOver', sep='/') )
      }
    }
    
    ###########################################################################
    
    ### MaxEnt things happen here
    
    ### preparing the data for the maxent model
    # filtering the occ object by it's respective cross-validation identity
    cv.number <- eval.df$cv.num[i]
    # assigning the column number to be sent through maxent
    col.number <- first.occ.col + cv.number - 1
    # filtering the species dataset by col.number and assigning to eval.df
    sp <- occs %>% dplyr::filter(occs[,col.number] == 1)
    n <- nrow(sp)
    eval.df[i,8] <- n
    
    # generating the testing data
    
    test.col.number <- first.test.col + cv.number - 1
    sp.test <- occs %>% dplyr::filter(occs[,test.col.number] == 1)
    
    
    # bending the occ and background data together
    mx.data <- rbind(sp, background)
    
    ### running the actual maxent model
    mx.model <- dismo::maxent(p = mx.data[,col.number],
                              x = mx.data[,predic],
                              path =  paste0(getwd()),
                              args = prepPara(userfeatures = eval.df[i,6], betamultiplier = eval.df[i,7], doclamp = FALSE)
    )
    
    ### calculating k and the names of the lambdas
    lambda.file <- as.data.frame(mx.model@lambdas) %>% `colnames<-`('lambdas')
    # lambdas data frame
    lambdas.df <- splitstackshape::cSplit(lambda.file, sep=',', splitCols='lambdas')
    colnames(lambdas.df) <- c('feature', 'lambda', 'min', 'max')
    # finding the non-zero lambdas
    non.zero.lambdas <- lambdas.df %>% dplyr::filter(!is.na(max)) %>% dplyr::filter(lambda != 0)
    class(non.zero.lambdas) <- 'data.frame'
    # giving noting the variables/features/hyperparameters with non-zero lambdas
    eval.df[i,17] <- toString( non.zero.lambdas[,1] )
    
    
    ## evaluating the models
    # predicting the training data
    mx.train.occ <-  dismo::predict(mx.model, sp[,predic])
    # predicting the training background data
    mx.train.back <- dismo::predict(mx.model, background[,predic])
    # predicting the testing data
    mx.test.occ <- dismo::predict(mx.model, sp.test[,predic])
    # projecting to all extents
    mx.test.back <- dismo::predict(mx.model, all.background[,predic])
    # projecting to all occ points
    mx.all.occs <- dismo::predict(mx.model, occs[,predic])
    
    # calculate the threshold
    if(omission.rate == 0){   # if using the LTP threshold, dismo calculates slightly too high of a threshold
      thresh <- min(mx.train.occ)
    } else {
      # make a model evaluation object
      mx.eval <- dismo::evaluate(p=mx.train.occ, a=mx.train.back)
      thresh <- dismo::threshold(mx.eval, stat='sensitivity', sensitivity= (1 - omission.rate) )
    }
    # assigning the threshold value to the output file
    eval.df[i,9] <- thresh
    # calculating the test sensitivity
    sens <- length(which(mx.test.occ >= thresh)) / length(mx.test.occ)
    eval.df[i,10] <-  sens
    
    # making a raster of the testing background data
    # for some reason, kuenm calculates wonky AUC_ratio values (i.e., > 2, which is impossible) unless you specify a raster
    r <- raster(nrows=1, ncols=length(mx.test.back) )
    r[r] <- mx.test.back
    
    ## calculating AUC_ratios from a partialROC test
    # error = 0.1%
    pROC_0.1 <- kuenm::kuenm_proc(occ.test = mx.test.occ,   # numeric vector of the predicted suitability values on the testing data
                                  model = r,                # raster model of the predicted suitability values for the background
                                  threshold = 0.1,          # potential error threshold (expressed as a percent)
                                  rand.percent = 50,        # percentage of data to be used in each bootstrap rep
                                  iterations = pROC.reps    # number of repititions
    )
    # assigning the average AUC_ratio from pROC.reps iterations to eval.df
    eval.df[i,11] <- as.numeric(pROC_0.1$pROC_summary[1])
    # assigning the partialROC p-value to eval.df
    eval.df[i,12] <- as.numeric(pROC_0.1$pROC_summary[2])
    
    # error = 1%
    pROC_1 <- kuenm::kuenm_proc(occ.test = mx.test.occ,   # numeric vector of the predicted suitability values on the testing data
                                model = r,                # raster model of the predicted suitability values for the background
                                threshold = 1,          # potential error threshold (expressed as a percent)
                                rand.percent = 50,        # percentage of data to be used in each bootstrap rep
                                iterations = pROC.reps    # number of repititions
    )
    # assigning the average AUC_ratio from pROC.reps iterations to eval.df
    eval.df[i,13] <- as.numeric(pROC_1$pROC_summary[1])
    # assigning the partialROC p-value to eval.df
    eval.df[i,14] <- as.numeric(pROC_1$pROC_summary[2])
    
    # error = 5%
    pROC_5 <- kuenm::kuenm_proc(occ.test = mx.test.occ,   # numeric vector of the predicted suitability values on the testing data
                                model = r,                # raster model of the predicted suitability values for the background
                                threshold = 5,          # potential error threshold (expressed as a percent)
                                rand.percent = 50,        # percentage of data to be used in each bootstrap rep
                                iterations = pROC.reps    # number of repititions
    )
    # assigning the average AUC_ratio from pROC.reps iterations to eval.df
    eval.df[i,15] <- as.numeric(pROC_5$pROC_summary[1])
    # assigning the partialROC p-value to eval.df
    eval.df[i,16] <- as.numeric(pROC_5$pROC_summary[2])
    
    
    mx.test.back.df <- as.data.frame( as.matrix(mx.test.back, ncol=1) )
    
    # assigning the objects to the various lists
    maxent.list[[i]] <- mx.model
    test.back.list[[i]] <- mx.test.back
    all.occs.list[[i]] <- mx.all.occs
    
    # updating the prograss bar for each run to get an idea of how long things will take
    setTxtProgressBar(prog, i)
    
    
    ###########################################################################
    # returning to the home directory
    setwd(home)
    
  }   # closes the for loop
  
  # bind the testing background data into a single data
  back.projections <- as.data.frame( do.call('cbind', test.back.list ) )
  colnames(back.projections) <- eval.df$CrossVal
  
  # bind all occs together
  occ.projections <- as.data.frame(do.call('cbind', all.occs.list))
  colnames(occ.projections) <- eval.df$CrossVal
  
  # make a list of the output
  out.list <- list()
  out.list$maxent.models <- maxent.list
  out.list$back.projection <- back.projections
  out.list$occ.projection <- occ.projections
  out.list$summary <- eval.df
  
  return(out.list)
}



##### MaxEnt Evaluation ----------------------------------------------------------------------------------------------------

# function that takes the eval object generated from maxent.crossval.error and calculates the weighted mean and standard deviation
# the weighted mean is calculating my testing sensitivity (1 - omission rate) multiplied by the partial ROC/AUC value.
# this ensures that if a model does not discriminate between presences/non-presences well, it will receive comparatively lower weight
# later package versions will include Boyce index as an additional calibration technique and will offer the user the ability to . . .
# . . . choose which metrics to use for assessing model reliability
# the weighted mean and standard deviation can be plotted to infer model variability/uncertainty

### ARGUMENTS ###
## eval <- eval object generated from maxent.crossval.error that will project the model to every grid cell used in the training region
## pROC.error <- the user-specified omission rate.
# (choose from 0.1, 1, 5 - however they almost always end up super correlated with each other)
maxent.eval <- function(eval, pROC.error=1){
  
  # extracting model sensitivity
  sens <- eval$summary$test.sens
  
  # which pROC error amount to use?
  if(pROC.error == 0.1){
    AUC_ratio <- eval$summary$pROC_0.1
  } else if(pROC.error == 1){
    AUC_ratio <- eval$summary$pROC_1
  } else if(pROC.error == 5){
    AUC_ratio <- eval$summary$pROC_5
  } else {
    stop('Select an appropriate partialROC error amount.')
  }
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * AUC_ratio
  weights[is.nan(weights)] <- 0
  
  # making a matrix of the background points
  models.mat <- as.matrix(eval$back.projection)
  # weighted means and standard deviations
  w.mean <- matrixStats::rowWeightedMeans(x=models.mat, w=weights)
  w.sd <- matrixStats::rowWeightedSds(x=models.mat, w=weights)
  # merging the weighted means/sd's of all background points
  back.df <- data.frame(w.mean, w.sd)
  
  # making a matrix of the background points
  occs.mat <- as.matrix(eval$occ.projection)
  # weighted means and standard deviations
  w.mean <- matrixStats::rowWeightedMeans(x=occs.mat, w=weights)
  w.sd <- matrixStats::rowWeightedSds(x=occs.mat, w=weights)
  # merging the weighted means/sd's of all occ points
  occ.df <- data.frame(w.mean, w.sd)
  
  
  # making the output list
  out.list <- list()
  out.list$back <- back.df
  out.list$occ <- occ.df
  out.list$weights <- weights
  
  return(out.list)
  
}



##### continuous.boyce ----------------------------


continuous.boyce <- function(model.pred,        # predicted suitability values for the entire model
                             occs.pred,         # predicted suitability values for the occurrences
                             window.size=0.1,   # size of the moving window (as a proportion of the suitability range)
                             reso=1000          # Resolution (in focals) of the moving window
){       # if TRUE, plot the predicted:expected ratio vs the suitability values
      
      # makes intervals of the data
      interval <- range(fit)
      # definin window.size
      window.size <- diff(interval)*window.size
      # uses a moving window
      vec.mov <- seq(from = interval[1],
                     to = interval[2] - window.size,
                     by = (interval[2] - interval[1] - window.size)/reso)  # 101 long
      vec.mov[reso + 1] <- vec.mov[reso + 1] + 1  #Trick to avoid error with closed interval in R
      interval <- cbind(vec.mov, vec.mov + window.size)
      
      # for loop calculating the continuous boyce window
      f.i <- c()
      
      numbs <- nrow(interval)
      
      for(i in 1:numbs){
            # make dummy variables that will get written over every loop
            model.pred.bin <- model.pred
            occs.pred.bin <- occs.pred
            # predicted amount
            model.pred.bin[model.pred[] >= interval[i,1] & model.pred[] <= interval[i,2]] <- "l"
            model.pred.bin[model.pred.bin != "l"] <- 0
            # expected amount
            occs.pred.bin[occs.pred[] >= interval[i,1] & occs.pred[] <= interval[i,2]] <- "l"
            occs.pred.bin[occs.pred.bin != "l"] <- 0
            # predicted and expected ratio (F)
            p.i <- length(which(occs.pred.bin == "l"))/length(occs.pred)
            e.i <- length(which(model.pred.bin == "l"))/length(model.pred.bin)
            f.i[i] <- p.i/e.i
      }  # close for(i in 1:numbs) loop
      
      # vector that keeps the non-NaN data
      no.nas <- which(f.i != "NaN")
      # removes non-NA values from F vector
      f.i <- f.i[no.nas]
      
      # Spearman's rho correlation
      if (length(f.i) < 2) {
            rho <- NA
      } else {
            # remove duplicate values
            no.dups <- c(1:length(f.i))[f.i != c(f.i[-1], FALSE)]
            # calculates spearman correlation (i.e. Boyce index) after removing duplicated values
            rho <- cor(f.i[no.dups], vec.mov[no.nas][no.dups], method = "spearman")
      }
      
      # habitat suitability
      HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
      HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
      HS <- HS[no.nas]  #exlude the NaN
      # results
      results <- list(F.ratio = f.i, Spearman.cor = round(rho, 3), HS = HS)
      return(results)
} # close function


##### maxent.eval.summarize -------------------------------------------------------------------------------------------------

# function that summarizes information across multiple cross-validation reps and find the "best" one using the 
# one-standard error rule

maxent.eval.summarize <- function(eval, pROC.error=0.1){
      require(dplyr)
      # eval summary
      eval.sum <- eval$summary
      ### getting model weights
      # extracting model sensitivity
      sens <- eval.sum$test.sens
      # which pROC error amount to use?
      if(pROC.error == 0.1){
            AUC_ratio <- eval.sum$pROC_0.1
      } else if(pROC.error == 1){
            AUC_ratio <- eval.sum$pROC_1
      } else if(pROC.error == 5){
            AUC_ratio <- eval.sum$pROC_5
      } else {
            stop('Select an appropriate partialROC error amount.')
      }
      # Boyce Index
      
      # weights = sensitivity*AUC_ratio
      AUC_ratio[is.nan(AUC_ratio)] <- 0
      AUC_ratio[is.na(AUC_ratio)] <- 0
      weights <- sens * AUC_ratio
      eval.sum$weights <- weights
      eval.sum$AUC_ratio <- AUC_ratio
      # grouping by Features and Betas and calculating the performance weight
      eval.group <- eval.sum %>% group_by(Features, Betas) %>% summarize(across(c(test.sens, AUC_ratio, weights), mean))
      # calculating the one-standard error rule
      nmodels <- nrow(eval.group)
      SE <- sd(eval.group$weights)/(nmodels^0.5)
      oneMinusSE <- max(eval.group$weights) - SE
      eval.group$standard1mSE <- sapply(eval.group$weights, function(x) abs((x-oneMinusSE)/SE) )
      # returning eval.group
      return(eval.group %>% arrange(standard1mSE))
}



##### Calculating the MaxEnt Threshold -------------------------------------------------------------------------------------

# function that calculates the threshold value for all the training data

### ARGUMENTS ###
## eval <- eval object generated from maxent.crossval.error that will project the model to every grid cell used in the training region
## occs <- the occurrence data.frame
## predic <- column numbers of the predictor variables 
## pROC.error <- the user-specified omission rate.
# (choose from 0.1, 1, 5 - however they almost always end up super correlated with each other)

# function will return the lowest training preference threshold.
# future package versions will give the opportunity to select different thresholds

maxent.thresh <- function(eval, occs, predic, pROC.error=1){
  
  n.rows <- nrow(occs)
  n.cols <- nrow(eval$summary)
  
  # making a data.frame of the occurrences
  occ.values <- as.data.frame( matrix(NA, nrow=n.rows, ncol=n.cols ) )
  colnames(occ.values) <- eval$summary$CrossVal
  
  # for loop predicting all the cross validation maxent models
  for(i in 1:nrow(eval$summary) ){
    mx.occs <- dismo::predict(object=eval$maxent.models[[i]], x=occs[,predic])
    occ.values[,i] <- mx.occs
  }
  
  # 
  # extracting model sensitivity
  sens <- eval$summary$test.sens
  
  # which pROC error amount to use?
  if(pROC.error == 0.1){
    AUC_ratio <- eval$summary$pROC_0.1
  } else if(pROC.error == 1){
    AUC_ratio <- eval$summary$pROC_1
  } else if(pROC.error == 5){
    AUC_ratio <- eval$summary$pROC_5
  } else {
    stop('Select an appropriate partialROC error amount.')
  }
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * AUC_ratio
  weights[is.nan(weights)]<- 0
  
  
  # making a matrix of the background points
  models.mat <- as.matrix(occ.values)
  # weighted means
  w.mean <- matrixStats::rowWeightedMeans(x=models.mat, w=weights)
  
  # setting the threshold as the lowest training presence
  ltp <- min(w.mean)
  
  return(ltp)
}



##### Projecting the Best Models to All the Extents ------------------------------------------------------------------------

# function to project models to any and all extents you want to
# effectively making new master occurrence and master background files that can easily be saved and re-loaded again

maxent.everything <- function(eval, means, everything, thresh=NULL, pROC.error=1, predic, name='pc'){
  
      
  # making a data.frame of the occurrences
  every.value <- as.data.frame( matrix(NA, nrow=nrow(everything), ncol=nrow(eval$summary)) )
  colnames(every.value) <- eval$summary$CrossVal
  
  # for loop predicting all the cross validation maxent models
  for(i in 1:nrow(eval$summary) ){
    mx.everything <- dismo::predict(object=eval$maxent.models[[i]], x=everything[,predic])
    every.value[,i] <- mx.everything
  }
  
  # extracting model sensitivity
  sens <- eval$summary$test.sens
  
  # which pROC error amount to use?
  if(pROC.error == 0.1){
    AUC_ratio <- eval$summary$pROC_0.1
  } else if(pROC.error == 1){
    AUC_ratio <- eval$summary$pROC_1
  } else if(pROC.error == 5){
    AUC_ratio <- eval$summary$pROC_5
  } else {
    stop('Select an appropriate partialROC error amount.')
  }
  
  # weights = sensitivity*AUC_ratio
  weights <- sens * AUC_ratio
  weights[is.nan(weights)] <- 0
  
  # making a matrix of the background points
  models.mat <- as.matrix(every.value)
  # weighted means
  w.mean <- matrixStats::rowWeightedMeans(x=models.mat, w=weights)
  
  # making the output data frame
  out.df <- as.data.frame( matrix(NA, nrow=length(w.mean), ncol=2) )
  colnames(out.df) <- c(paste0(eval$summary[1,1], '_contin'), paste0(eval$summary[1,1], '_thresh') )
  # assigning the weighted means to the output data.frame
  out.df[,1] <- w.mean
  # replicating the weighted means so they can be thresolded
  out.df[,2] <- w.mean
  # calculating the threshold - if no threshold is provided, use the lowest training presence
  if(is.null(thresh)){
        thresh <- min(means$occ$w.mean)
        warning('No threshold provided.\nUsing lowest weighted training presence as the threshold.\nThis value is likely arbitrarily low.')
  }
  # thresolding the values
  out.df[,2][out.df[,2] >= thresh] <- 1
  out.df[,2][out.df[,2] < thresh] <- 0
  
  
  
  return(out.df)
  
}



##### Informed MESS Analysis -----------------------------------------------------------------------------------------------

# function that calculates Multivariate Environmental Suitability Surface (MESS), but allows for some extrapolation . . .
# . . . based on how response curves look.

informed.mess <- function(ref.extent,        # reference data.frame - must contain same column names as 'data'
                          mess.extent,        # data.frame containing the extent to be projected to - ideally this should be all extents
                          coord.cols=1:2,    # column numbers of the coordinates of the data/ref.data (e.g., long/lat)
                          predic=3:5,    # column numbers of ONLY the final predictor variables
                          tolerance=NULL,    # tolerance vector to extrapolate beyond the limits of the data. see example below
                          # if tolerance is not specified (default), function will calculate basic MESS
                          keep.layers=FALSE  # decision to retain the MESS values for all variables instead of just the final MESS data
                          # changing to TRUE can help easily identify which variables are causing the MESS value . . .
                          # . . . to be so low (i.e., it may be only 1 variable that is responsible)
                          # can easily re-check this later
){
  # isolating the coordinates and predictor variables
  new.vars <- mess.extent[, predic]
  ref.vars <- ref.extent[, predic]
  # 
  if( !is.null(tolerance) ){
    if(length(predic) != length(tolerance)/2 ){
      stop('Ensure that tolerance is 2x as long as predic')
    }
    # changing tolerance to a data.frame
    tolerance <- as.data.frame( matrix(tolerance, ncol=2, byrow=TRUE) )
    nvars <- length(predic)
    # pre-changing values outside the range of the input data
    # extra rows
    extra.rows <- matrix(NA, nrow=2, ncol=nvars)
    colnames(extra.rows) <- colnames(ref.vars)
    ref.vars <- rbind(ref.vars, extra.rows )
    for(i in 1:nvars){
      # if we can extrapolate beyond the lower end of variable i, allow mess to not 
      if(tolerance[i,1] == 1){
        min.var <- min(ref.vars[,i], na.rm=TRUE) - 0.0001
        new.vars[ new.vars[,i] < min.var, i ] <- min.var
        ref.vars[nrow(ref.vars)-1, i] <- min.var
      } else {
        ref.vars[nrow(ref.vars)-1, i] <- min(ref.vars[,i], na.rm=TRUE)
      }
      
      if(tolerance[i,2] == 1){
        max.var <- max(ref.vars[,i], na.rm=TRUE) + 0.0001
        new.vars[ new.vars[,i] > max.var, i ] <- max.var
        ref.vars[nrow(ref.vars), i] <- max.var
      } else {
        ref.vars[nrow(ref.vars), i] <- max(ref.vars[,i], na.rm=TRUE)
      }
      
    } # end for(i in 1:nvars) loop
    
  }
  
  # running the mess analysis
  mess.vars <- as.data.frame( sapply(1:ncol(new.vars), function(i) .messi3(new.vars[, i], ref.vars[, i])) )
  # re-asigning of the extrapolating points to have a mess value of 0. 
  nref <- nrow(ref.vars)
  fix.vars <- as.data.frame( sapply(1:ncol(new.vars), function(i) .messi3(ref.vars[ (nref-1):nref , i], ref.vars[, i])) )
  for(i in 1:nvars){
    mess.vars[mess.vars[,i] == fix.vars[1,i],i] <- 0
    mess.vars[mess.vars[,i] == fix.vars[2,i],i] <- 0
    
  }
  
  # making a simple thresholded version of the mess analysis for easy plotting
  final.mess <- as.data.frame( apply(mess.vars, 1, min) )
  colnames(final.mess) <- 'mess.raw'
  mess.thresh <- final.mess$mess.raw
  mess.thresh[mess.thresh > 0] <- 1
  mess.thresh[mess.thresh < 0] <- -1
  final.mess <- cbind(final.mess, mess.thresh)
  
  # deciding whether to keep individual mess layers
  if(keep.layers==TRUE){
    colnames(ref.vars) <- paste0('mess_', colnames(ref.vars))
    final.mess <- cbind(final.mess, ref.vars)
    return(final.mess)
  } else {
    return(final.mess)
  }
  
}



# internal function originally from dismo that runs the actual mess analysis in data.frame format
.messi3 <- function(p,v) { # p=new.vars   v=ref.vars
  # seems 2-3 times faster than messi2
  v <- stats::na.omit(v)
  f <- 100*findInterval(p, sort(v)) / length(v)
  minv <- min(v)
  maxv <- max(v)
  res <- 2*f 
  f[is.na(f)] <- -99
  i <- f>50 & f<100
  res[i] <- 200-res[i]
  
  i <- f==0 
  res[i] <- 100*(p[i]-minv)/(maxv-minv)
  i <- f==100
  res[i] <- 100*(maxv-p[i])/(maxv-minv)
  res
}




##### Informed MESS conceptually -------------------------------------------------------------------------------------------

# principles of informed mess analysis


#            |                                                      #     # we haven't captured the entire response curve for this . . .
#      S   1-|          |                      |                 __ #     # . . . variable. this is a problem when extrapolating . . .
#      u     |                                           ?__-----   #     # . . . to non-analog environments
#      i     |          |                      |    ?___/           #
#      t     |                                  ___/                #     # we're very confident that projecting to environments with . . .
#      i     |          |                      |                    #     # . . . very low temperatures would likely mean low suitability
#      b     |                             ___  ____?____?____?____ #
#      i     |          |              ___/    |                    #     # however, we aren't confident about this at high temperatures.
#      l     |                      __/           \__?              #     # the suitability could tail off at any point, but the model . . .
#      i     |          |        __/           |      \__?          #     # . . . will likely assume ↑temp = ↑suitability forever
#      t     |                 _/                         \_?       #
#      y   0-|  ?___?___|_____/                |              \__?_ #     # to note this in the informed.mess function, you need to . . .
#            |_____________________________________________________ #     # . . . visualize the response curves in the maxent output file
#                       |       Variable       |                    #
#                                  1                                #     # note this in the 'tolerance' argument
#                       |  (e.g., Temperature) |                    #



### the tolerance argument is a vector telling the function what to do with the non-analog environments it encounters
## tolerance is a vector containing 1's and 0's based on what the response curves look like

### specify 1 or 0 for both ends of the response curves for all the final variables in the model
## check the lambdas column in the output summary file for the final variables
# for example, if the lambdas for a madel are: coldmon^2, warmmon^2, & coldmon*warmmon, . . . 
# . . . there are 2 final variables in the model (coldmon and warmmon), so tolerance should be 4 numbers long
# for another example, if the lambdas for a madel are: coldmon^2, warmmon^2, coldmon*warmmon, & seasonality, . . . 
# . . . there are 3 final variables in the model (coldmon, warmmon, and seasonality), so tolerance should be 6 numbers long
# for a final example, if the lambdas for a madel are: coldmon & coldmon^2,, . . . 
# . . . there is 1 final variable in the model (coldmon), so tolerance should be 2 numbers long


### the following example contains three final variables (Temperature, Precip, and Seasonality)
# thus, tolerance should be a vector with length 6
# tolerance <- c(_,_,
#                _,_,
#                _,_)

#            |                                                      #     # because we ARE confident in extrapolating to . . .
#      S   1-|          |                      |                    #     # . . .  low temps, put a 1 for the 1st spot in tolerance
#      u     |                                                      #
#      i     |          |                      |                    #     # tolerance <- c(1,_,
#      t     |    1                                       0         #     #                _,_,
#      i     |          |                      |                    #     #                _,_)
#      b     |                             ___                      #
#      i     |          |              ___/    |                    #     # however, because we're NOT confident in extrapolating . . .
#      l     |                      __/                             #     # . . . to high temps, put a 0 in the 2nd spot in tolerance
#      i     |          |        __/           |                    #
#      t     |                 _/                                   #     # tolerance <- c(1,0,
#      y   0-|          |_____/                |                    #     #                _,_,
#            |_____________________________________________________ #     #                _,_)
#                       |       Variable       |                    #
#                                  1                                #
#                       |  (e.g., Temperature) |                    #



#            |                                                      #     # because we ARE confident in extrapolating to . . .
#      S   1-|          |                      |                    #     # . . . low precip, put a 1 for the 3rd spot in tolerance
#      u     |                                                      #
#      i     |          |                      |                    #     # tolerance <- c(1,0,
#      t     |    1                                       1         #     #                1,_,
#      i     |          |                      |                    #     #                _,_)
#      b     |                                                      #
#      i     |          |          ___         |                    #     # because we ARE ALSO confident about extrapolating to . . .
#      l     |                  __/   \__                           #     # . . . high precip, put a 1 in the 4th spot in tolerance
#      i     |          |    __/         \__   |                    #
#      t     |             _/               \__                     #     # tolerance <- c(1,0,
#      y   0-|          |_/                    |                    #     #                1,1,
#            |_____________________________________________________ #     #                _,_)
#                       |       Variable       |                    #
#                                  2                                #
#                       |    (e.g., Precip)    |                    #



#            |                                                      #     # given that this response curve does NOT turn downwards on . . .
#      S   1-|          |                      |                    #     # . . . the low end, you COULD put a 0 for the 5th spot in . . .
#      u     |                                                      #     # . . . seasonality. HOWEVER, low seasonality is unlikely to . . .
#      i     |          |                      |                    #     # . . . be non-suitable for any species, we can assume that . . .
#      t     |    1                                      1          #     # . . . this extrapolation will NOT be problematic.
#      i     |          |                      |                    #
#      b     |           ___                                        #     # tolerance <- c(1,0,
#      i     |          |   \___               |                    #     #                1,1,
#      l     |                  \___                                #     #                1,_)
#      i     |          |           \__        |                    #
#      t     |                         \_                           #     # because we ARE confident about extrapolating to at high . . .
#      y   0-|          |                \_____|                    #     # . . . seasonality, put a 1 in the 6th spot in tolerance
#            |_____________________________________________________ #
#                       |       Variable       |                    #     # tolerance <- c(1,0,
#                                  3                                #     #                1,1,
#                       | (e.g., Seasonality)  |                    #     #                1,1)


### basically, a 0 tells the code to believe the MESS value that it calculates, whereas a 1 tells the code to . . .
#   . . . change the MESS value because we have confidence that the projection will still be reasonable

### thus, for the previous three response curves, specify the following for mess.tol
# mess.tol <- c(1,0,
#               1,1,
#               1,1)





##### calculate.mop --------------------------------------------------------------------------------------------------------


## calculate.mop calculates Mobility Oriented Parity (MOP; sensu Owens et al. 2013) from a set of reference points

# function identifies env-space that represents novel multivariate combinations relative to a set of reference points
# use this in your analyses if you want to be conservative about projecting ENMs to new areas/time periods


### ARGUMENTS
## ref.extent <- a data.frame or matrix containing the variable values for the reference/training extent, such as . . .
#                . . . any object generated with the create.background.df function
## predic <- the column numbers of the predictor variables in ref.extent
## new.extent <- a data.frame or matrix containing the variable values for the new extent to be projected to, such as . . .
#                . . . any object generated with the create.background.df function
## if new.extent is not provided, function will calculate MOP based on distance from each point in ref.extent to . . .
#  . . . the nearest set of points in ref.extent
## if new.extent is provided, function will calculate MOP based on distance from each point in new.extent to . . .
#  . . . the nearest set of points in ref.extent


## the function returns a data.frame of multivariate distances between a given set of points and other points in the . . .
#  . . . training extent with the same number of rows as either new.extent (if provided) or ref.extent
## distances are multivariate euclidean distances based on a set of reference points that have been standardized (μ=0, σ=1)
## it is recommended that the user run MOP on the ref.extent first to set a threshold for how much . . .
#  . . . environmental difference from ref.extent they are willing to tolerate
## suggested thresholds are near the 
#  note that this threshold is a trade-off. if ref.extent is small (i.e., few grid cells), the user may want to select . . .
#  . . . a higher MOP threshold as there may be greater chance for having holes/gaps in env-space.
#  thus, selecting an artificially low MOP threshold is likely to filter out interpret-able parts of env-space
#  this will also be an issue with the greater number of dimensions ref.extent has.
#  more dimensions = more chance for multivariate holes in env-space. these are often super hard/impossible to detect
## in the output data.frame, the column names are:
# pt1 = distance from a given point to the nearest point in env-space
# pt5 = distance from a given point to the nearest 5 points in env-space
# pt10 = distance from a given point to the nearest 10 points in env-space
# pt1perc = distance from a given point to the nearest 1% of points in env-space
# pt5perc = distance from a given point to the nearest 5% of points in env-space


## worth noting that this will easily blow out R's memory allotment if the number of points is high (>50000) or . . .
#  . . . if the dimensionality is high
# this is because the number of calculations in the distance matrix that is calculated is equal to . . .
# . . . (nrow(ref.extent) + nrow(new.extent))^2
# for example, having 50000 total points results in 2.5B calculations, which takes 3 minutes per dimension and . . .
# . . . will take a stupidly long time to complete
# an example with ~34000 background points took ~13 minutes
# you may see an error like this:
# Error: vector memory exhausted (limit reached?)
# e.g., this failed when there were ~150000 background points over 4 dimensions, even when setting . . .
# . . . R_MAX_VSIZE=100Gb as an environment
# example from https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached 

calculate.mop <- function(ref.extent, predic, new.extent=NULL){
      # initial condition checking to see if the number of background points is >50000
      if(!is.null(new.extent)){
            total.pts <- nrow(ref.extent) + nrow(new.extent)
      } else {
            total.pts <- nrow(ref.extent)
      }
      if(total.pts >= 50000){
            stop('Calculating ', total.pts, ' points will blow out the memory. Please subset total points to < 50000.')
      }
      # setting start time
      start.time <- Sys.time()
      message( paste0( 'The time is ', start.time ) )
      # number of variables ran
      n.vars <- length(predic)
      # number of background points
      n.points <- nrow(ref.extent)
      # reducing the data to only the variables of interest for the ref extent
      ref.extent <- ref.extent[,predic]
      # defining mop.extent if new.extent is provided
      # if new.extent is not provided, define mop.extent as ref.extent
      if( !is.null(new.extent) ){
            # number of points in new extent
            n.new.points <- nrow(new.extent)
            # reducing the data to only the variables of interest for the new extent
            new.extent <- new.extent[,predic]
            # ensuring that the columns of ref.extent and new.extent match
            # if so, bind them together (new.extent 1st!). if not, return an error
            if( identical( colnames(ref.extent), colnames(new.extent) ) ){
                  mop.extent <- rbind(new.extent, ref.extent)
            } else {
                  stop(paste0( ref.extent,' and ', new.extent,' have non-identical columns'  ) )
            }
      } else {
            mop.extent <- ref.extent
            n.new.points <- n.points
      }
      # calculating the mean and sd of each variable
      rescale.df <- as.data.frame( matrix(NA, nrow=2, ncol=n.vars) )
      for( i in 1:n.vars ){
            rescale.df[1,i] <- mean(ref.extent[,i])   # mean of each variable
            rescale.df[2,i] <- sd(ref.extent[,i])   # sd of each variable
      }
      # for loop standardizing each of the variables in both extents
      for (i in 1:ncol(ref.extent) ){
            # center each variable on 0
            mop.extent[,i] <- mop.extent[,i] - rescale.df[1,i]
            # divide each variable by the standard deviation
            mop.extent[,i] <- mop.extent[,i]/rescale.df[2,i]
      }
      # defining the output data.frame that is n.new. points long and measures distances from any point to nearby points
      output.df <- as.data.frame(matrix(NA, nrow=n.new.points, ncol=5))
      colnames(output.df) <- c('pt1',   # distance from each point to the nearest 1 point
                               'pt5',   # average distance from each point to the nearest 5 points
                               'pt10',   # average distance from each point to the nearest 10 points
                               'pt1perc',   # average distance from each point to the nearest 1% of points
                               'pt5perc')   # average distance from each point to the nearest 5% of points
      # defining a distance matrix between every background point and every other background point (for mop.extent)
      start.d <- Sys.time()
      message('Calculating distance matrix')
      d.matrix <- fields::rdist(mop.extent)
      end.d <- Sys.time()
      # if statement contingent on new.extent being provided
      if( is.null(new.extent) ){
            # defining the diagonal and lower half of the matrix as NA values
            # doing this because the distance between points a and b in the matrix is noted in cells [a,b] and [b,a] and . . .
            # . . . we don't want to use that distance twice 
            diag(d.matrix) <- NA
            d.matrix[upper.tri(d.matrix)] <- NA
            # defining the start time
            mop.start <- Sys.time()
            # setting up the progress bar
            message( paste0('Calculating MOP for the reference extent') )
            prog <- txtProgressBar(min=0, max=nrow(output.df), style=3,  char='+')
            # massive for loop calculating the distance from a given point to the nearest X points
            for( i in 1:n.new.points ){
                  # the distance from point i to every other point is noted in d.matrix[i,] and d.matrix[,i]
                  # filters out distance from point i to any other point and sorts if
                  d.point <- sort( c( d.matrix[i,], d.matrix[,i] ) )
                  # takes the ranked distance(s) from point i to every other point and averages the closest points
                  output.df[i,1] <- d.point[1]   # distance from each point to the nearest 1 point
                  output.df[i,2] <- mean(d.point[1:5])   # average distance from each point to the nearest 5 points
                  output.df[i,3] <- mean(d.point[1:10])   # average distance from each point to the nearest 10 points
                  output.df[i,4] <- mean(d.point[1:floor(n.points/100)])   # average distance from each point to the nearest 1% of points
                  output.df[i,5] <- mean(d.point[1:floor(n.points/20)])   # average distance from each point to the nearest 5% of points
                  # updating the progress bar for each run to get an idea of how long things will take
                  setTxtProgressBar(prog, i)
            }
      } else {
            # filter d.matrix to only distances between new.extent and points in the old extent
            d.matrix <- d.matrix[ -(1:n.new.points) , 1:n.new.points ]
            # defining the start time
            mop.start <- Sys.time()
            # setting up the progress bar
            message( paste0('Calculating MOP for the new extent') )
            prog <- txtProgressBar(min=0, max=nrow(output.df), style=3,  char='+')
            # massive for loop calculating the distance from a given point to the nearest X points
            for( i in 1:n.new.points ){
                  # the distance from point i to every other point is noted in d.matrix[i,] and d.matrix[,i]
                  # filters out distance from point i to any other point and sorts if
                  d.point <- sort(d.matrix[,i])
                  # takes the ranked distance(s) from point i to every other point and averages the closest points
                  output.df[i,1] <- d.point[1]   # distance from each point to the nearest 1 point
                  output.df[i,2] <- mean(d.point[1:5])   # average distance from each point to the nearest 5 points
                  output.df[i,3] <- mean(d.point[1:10])   # average distance from each point to the nearest 10 points
                  output.df[i,4] <- mean(d.point[1:floor(n.points/100)])   # average distance from each point to the nearest 1% of points
                  output.df[i,5] <- mean(d.point[1:floor(n.points/20)])   # average distance from each point to the nearest 5% of points
                  # updating the progress bar for each run to get an idea of how long things will take
                  setTxtProgressBar(prog, i)
            }
      } # closes ifelse is.null(new.extent) statement
      # returning how long it took to run aspects of the function
      # end time
      end.time <- Sys.time()
      # calculating distance matrix
      message(paste0('Time to calculate distance matrix: ', round(end.d - start.d, 5), ' ', units(end.d - start.d)) )
      # running MOP
      message( paste0('Time to run MOP: ', round(end.time - mop.start, 5), ' ', units(end.time - mop.start)) )
      # total function time
      message( paste0('Time to run function: ', round(end.time - start.time, 5), ' ', units(end.time - start.time)) )
      # return output
      return(output.df)
}


##### Variable Ranges ---------------------------------------------------

# function that returns the minimum and maximum ranges of each variable for a given set of cross-validation models
get.var.ranges <- function(eval){ # eval object generated from the maxent.crossval.error function
      r <- apply(X=eval$maxent.models[[1]]@absence, 2, 'range')
      colnames(r) <- colnames(eval$maxent.models[[1]]@absence)
      rownames(r) <- c('min', 'max')
      return(r)
}


##### Variable Importance ---------------------------------------------------

# function that returns the weighted variable contribution and importance for a given set of cross-validation models
# per Phillips 2006, 
# contribution depends on the process of how the model is fit. intrepret values of highly-correlated variables with caution
# importance is calculated based on randomizing a single variable at a time AFTER the model is made and calculating how much the AUC drops
get.var.contrib.importance <- function(eval, pROC.error = 5){ # eval object generated from the maxent.crossval.error function
      ### getting model weights
      # extracting model sensitivity
      sens <- eval$summary$test.sens
      # which pROC error amount to use?
      if(pROC.error == 0.1){
            AUC_ratio <- eval$summary$pROC_0.1
      } else if(pROC.error == 1){
            AUC_ratio <- eval$summary$pROC_1
      } else if(pROC.error == 5){
            AUC_ratio <- eval$summary$pROC_5
      } else {
            stop('Select an appropriate partialROC error amount.')
      }
      # weights = sensitivity*AUC_ratio
      weights <- sens * AUC_ratio
      weights[is.nan(weights)] <- 0
      weights[is.na(weights)] <- 0
      # number of cross-validation reps
      nreps <- length(eval$maxent.models)
      # number of variables
      nvars <- ncol(eval$maxent.models[[1]]@presence)
      # making data.frames for variables and importance
      contrib.df <- matrix(NA, nrow=nvars, ncol=nreps)
      colnames(contrib.df) <- paste0('cv_', 1:nreps)
      rownames(contrib.df) <- sort(colnames(eval$maxent.models[[1]]@presence)) # maxent sorts variables as strigs, so we have to here also
      impor.df <- contrib.df
      # for loop extracting the variable contribution and importance
      for(i in 1:nreps){
            contrib.df[,i] <- eval$maxent.models[[i]]@results[grep("contribution", rownames(eval$maxent.models[[i]]@results)), ]
            impor.df[,i] <- eval$maxent.models[[i]]@results[grep("importance", rownames(eval$maxent.models[[i]]@results)), ]
      }
      
      # weighted means and sd
      w.contrib <- data.frame(w.mean = matrixStats::rowWeightedMeans(x=contrib.df, w=weights),
                              w.sd = matrixStats::rowWeightedSds(x=contrib.df, w=weights),
                              row.names = rownames(contrib.df)
      )
      w.impor <- data.frame(w.mean = matrixStats::rowWeightedMeans(x=impor.df, w=weights),
                            w.sd = matrixStats::rowWeightedSds(x=impor.df, w=weights),
                            row.names = rownames(contrib.df)
      )
      
      out.list <- list(w.contrib, w.impor)
      names(out.list) <- c('VarContribution', 'VarImportance')
      return(out.list)
      
}



##### Extract Response Curves ---------------------------------------------------

# function will give the response curves for each variable while keeping all other variables at their median values for the training region
get.maxent.response.curves <- function(eval, # eval object generated from the maxent.crossval.error function
                                       expand=5, # extrapolation percentage - how far beyond the range of input values to project response curves to
                                       ngradient=100, # number of points used to calculate the response curve gradient
                                       pROC.error = 5 #  user-specified omission rate. choose from 0.1, 1, 5
){
      ### getting model weights
      # extracting model sensitivity
      sens <- eval$summary$test.sens
      # which pROC error amount to use?
      if(pROC.error == 0.1){
            AUC_ratio <- eval$summary$pROC_0.1
      } else if(pROC.error == 1){
            AUC_ratio <- eval$summary$pROC_1
      } else if(pROC.error == 5){
            AUC_ratio <- eval$summary$pROC_5
      } else {
            stop('Select an appropriate partialROC error amount.')
      }
      # weights = sensitivity*AUC_ratio
      weights <- sens * AUC_ratio
      weights[is.nan(weights)] <- 0
      weights[is.na(weights)] <- 0
      ### preparing arguments
      # number of cross-validation reps
      nreps <- length(eval$maxent.models)
      # make a data.frame of the output
      response.df <- as.data.frame(matrix(NA, nrow=ngradient, ncol=nreps))
      colnames(response.df) <- paste0('cv', seq(1:ncol(response.df)))
      # make background points for the gradient
      background <- eval$maxent.models[[1]]@absence
      cn <- colnames(background)
      # number of vars
      nvars <- ncol(background)
      # structure of the data.frame that will be used to get values for each response curve
      m <- matrix(nrow=1, ncol=ncol(background))
      # defining medians
      m <- as.numeric(apply(background, 2, median)) # assigns median point of response curve
      medians <- as.data.frame(matrix(m, nrow=ngradient, ncol=length(m), byrow=TRUE))
      colnames(medians) <- cn
      input.range <- as.data.frame(matrix(NA, nrow=ngradient, ncol=length(m), byrow=TRUE)) # input ranges of response curves
      colnames(input.range) <- cn
      # ranges of input variables
      ranges <- apply(X=background, 2, 'range')
      # creating the output.list
      out.list <- list()
      # filling out input.range
      for(i in 1:length(m)){
            # min and max of the variable
            minmax <- ranges[,i]
            r <- minmax[2]-minmax[1]
            expand1 <- r*round(expand)/100
            # assigning ranges to each variable number
            input.range[,i] <- minmax[1] - expand1 + 0:(ngradient-1) * (r + 2*expand1)/(ngradient-1)
      } # close length(cn) loop
      # nested for loop extracting each variable from each cross-valication rep
      for(var in 1:nvars){ # going variable by variable, then rep by rep
            # making a list of cross-val reps
            list.reps <- list()
            # for loop analyzing the cross-val reps
            for(rep in 1:nreps){
                  # making a data.frame of medians, but changing the var in question to the values of input.range
                  mx <- medians
                  mx[,var] <- input.range[,var]
                  # projecting maxent model to mx and assigning it to list.reps
                  list.reps[[rep]] <- dismo::predict(eval$maxent.models[[rep]], mx)
            }   # closing for(rep in nreps)
            
            # binding everything together
            reps <- do.call('cbind', list.reps)
            colnames(reps) <- paste0(cn[var], '_cv', 1:nreps)
            w.mean <- matrixStats::rowWeightedMeans(x=reps, w=weights)
            w.sd <- matrixStats::rowWeightedSds(x=reps, w=weights)
            reps <- as.data.frame(reps)
            reps$w.mean <- w.mean
            reps$w.sd <- w.sd
            colnames(reps) <-  c(paste0(cn[var], '_cv', 1:nreps), paste0(cn[var], '_mean'), paste0(cn[var], '_sd'))
            reps$input <- input.range[,var]
            
            out.list[[var]] <- reps
      }   # closing for(var in nvars)
      
      names(out.list) <- cn
      
      return(out.list)
}

##### Process Maxent Response Curves ---------------------------------------------------

# function that generates a list of data.frames that represent mean ±1 sd bounding boxes around uncertatinties
# input is an object generated by get.maxent.response.curves
process.maxent.response.curves <- function(rc){
      # number of variables
      nvars <- length(rc)
      # making output.list
      output.list <- list()
      # for loop making data.frames of the uncertainties of the variable contributions
      for(i in 1:nvars){
            
            lower <- data.frame(x=rc[[i]][,grep('input', colnames( rc[[i]] )) ],
                                y=rc[[i]][,grep('mean', colnames( rc[[i]] )) ] - rc[[i]][,grep('sd', colnames( rc[[i]] )) ]
            )
            upper <- data.frame(x=rc[[i]][,grep('input', colnames( rc[[i]] )) ],
                                y=rc[[i]][,grep('mean', colnames( rc[[i]] )) ] + rc[[i]][,grep('sd', colnames( rc[[i]] )) ]
            )
            upper <- upper[rev(rownames(upper)), ]
            rownames(upper) <- NULL
            both <- rbind(lower , upper)
            both[nrow(both) + 1,] <- both[1,]
            rownames(both) <- NULL
            # fixing any values below or above (0,1) to (0,1)
            both$y <- replace(both$y, both$y>1, 1)
            both$y <- replace(both$y, both$y<0, 0)
            
            output.list[[i]] <-  both
      }
      names(output.list) <- names(rc)
      return(output.list)
}


##### suit.uncert.plot -----------------------------------------------------------------------------------------------------

### makes a bivariate plot of weighted suitability and weighted uncertainty (standard deviation)

suit.uncert.plot <- function(means){
  require(ggplot2)
      
  back <- means$back
  occs <- means$occ
  
  ltp <- min(occs$w.mean)
  
  output <- ggplot(data=back, aes(x=w.mean, y=w.sd)) + geom_point(colour='black', size=0.75) +
    geom_point(data=occs, aes(x=w.mean, y=w.sd), size=2, colour='red', shape=18 ) +
    ylim(0, 0.55) + xlim(0,1) + theme_classic() + coord_fixed(1/0.55) +
    annotate(geom='text', x=0.05, y=0.45, label=round(ltp, 4), hjust=0 )
  
  return(output)
  
}



##### END ------------------------------------------------------------------------------------------------------------------




