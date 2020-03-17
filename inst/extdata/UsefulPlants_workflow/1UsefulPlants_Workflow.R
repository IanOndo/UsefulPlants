## If the OS is Windows, set mclapply to the hackish version. Otherwise, leave the definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {UsefulPlants::mclapply.hack},
                    Linux   = {parallel::mclapply},
                    Darwin  = {parallel::mclapply})


#====================================================================
#== 0. Get parameters form command-line
#====================================================================

args = commandArgs(trailingOnly=TRUE)

for(arg in args) {

  arg_name = stringr::str_extract(arg,".+?(?=\\=)")

  switch(arg_name,

         'base_directory' 			= {BASE_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'species_directory' 		= {SPECIES_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'env_directory' 			= {ENV_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'which_uses'    = {WHICH_USE = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'points_proj' 		= {POINTS_PROJ = paste0("+",gsub(".+?(?<=\\=\\+)","",arg,perl=TRUE))},

         'run_name'			= {RUN_NAME = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'mc_cores' 			= {MC_CORES = as.integer(gsub(".+?(?<=\\=)","",arg,perl=TRUE))}
  )
}
cat(paste0('#',paste(rep('-',times=100),collapse="")))
cat('\n')
cat(paste0("> Run ","'",RUN_NAME,"'",' started on: ',Sys.time()));started.at=proc.time()
cat('\n')
cat('#====================================================================\n')
cat('#== 1. Setup directories\n')
cat('#====================================================================\n\n')
xyOnly = FALSE
if(exists(ENV_DIR)){
  init_main_directories(base_directory=BASE_DIR, species_directory=SPECIES_DIR, env_directory=ENV_DIR, ncores=MC_CORES)
}else{
  xyOnly = TRUE
  init_main_directories(base_directory=BASE_DIR, species_directory=SPECIES_DIR, ncores=MC_CORES)
}
cat('\n\n')

if(!xyOnly){
  cat('#====================================================================\n')
  cat('#== 2. Prepare environmental layers\n')
  cat('#====================================================================\n')
  # TODO: write script to prepare envrionmental predictors
}

cat('#====================================================================\n')
cat('#== 3. Prepare species occurrence data\n')
cat('#====================================================================\n\n')

cat('#--------------------------------------------------------------------\n')
cat('### 3.a sort species by number of occurrence records\n')
cat('#--------------------------------------------------------------------\n')
# Sort species according to algorithms to be used. The function check the remaining number of sites/cells after cleaning (i.e. number of sites/cells where the species occur after spatial filtering)
# The spatial filtering distance is 20km (i.e. minimum distance between records is 20km)
# < 2 cells => points directory
# < 3 cells => bounding box directory
# < 10 cells => convex hull directory
# 20 > # of cells > 10  => SDM10_20
# > 20 cells => SDM

cat('\n\n')
if(!xyOnly){
  cat('#--------------------------------------------------------------------\n')
  cat('### 3.b refine the sorting but using grid background \n')
  cat('#--------------------------------------------------------------------\n')
  # TODO: write script to assign species whose # presence cells < 10 to the geographic distance-based folder
}

cat('#--------------------------------------------------------------------\n')
cat('### 3.b specify which species to run\n')
cat('#--------------------------------------------------------------------\n')

# List of species to run with MaxEnt algorithm
if(!xyOnly)
  speciesListMaxEnt=list.files(file.path(BASE_DIR,"envModels","maxent","occurrences"), full.names=T)
else
  speciesListMaxEnt=character(0)

# List of species to run with geographic distance-based algorithm
speciesListGeoDist=list.files(file.path(BASE_DIR,"geoModels","geo_dist","occurrences"), full.names=T)

# List of species to run as single points
speciesListPt=list.files(file.path(BASE_DIR,"pointModels","points","occurrences"), full.names=T)

# List of species to run
speciesList = c(speciesListPt, speciesListGeoDist, speciesListMaxEnt)

# List of algorithms to run
if(length(speciesList)==0L)
  stop("Oops !! Something probably went wrong with the species sorting...")

cat(paste(length(speciesList),'species will be run:\n'))
cat(paste('*',length(speciesListMaxEnt),"species with 'MaxEnt'\n"))
cat(paste('*',length(speciesListRB),"species with a 'Geographic distance-based' model\n"))
cat(paste('*',length(speciesListPt),"species will be modelled as 'Points'\n"))

cat('\n\n')

if(!xyOnly){
  cat('#====================================================================\n')
  cat('#== 4. Model settings\n')
  cat('#====================================================================\n')

  cat('...Setting up parameters for Maxent models...\n')

  maxent_settings <- list(
    path_to_maxent="C:/Users/io03kg/Desktop/MAXENT/maxent/maxent.jar",
    visible=FALSE,
    writemess=FALSE,
    writebackgroundpredictions=TRUE,
    maximumbackground=20000,
    betamultiplier=c(1,),
    threshold=FALSE,
    outputformat='raw',
    outputgrids=FALSE
  )
  cat(paste('Run settings saved in:',dirname(allBaseDirs$envDir)))
  cat('\n\n')
}else{
  cat('#====================================================================\n')
  cat('#== 4. Model settings\n')
  cat('#====================================================================\n')

  geographic_distance_algorithm = "dist"
  geographic_grid_resolution = 0.1666667

  cat('\n\n')
}

cat('#====================================================================\n')
cat('#== 5. Run Models\n')
cat('#====================================================================\n')

cat('...Starting computations...\n\n')
started.at = Sys.time()
mclapply(speciesList, function(j){

  # Get data for just that species
  if(file.exists(j)){
    species.data <- tryCatch(data.table::fread(j, showProgress=FALSE),error=function(err) return(NULL))
    k = gsub("\\.csv|_","",basename(j))
  }else{
    stop("Unable to locate the file :", j)
  }
  # if an error occurred during the reading returns an NULL
  if (is.null(species.data) || ncol(species.data) < 2){
    if(verbose) warning(paste0("Cannot read file from species",k))
    dn <- file.path(BASE_DIR,"problems")
    if(!dir.exists(dn))
      dir.create(dn, recursive =TRUE)
    fn <-file.path(dn,'Cannot_read_file.csv')
    if(!file.exists(fn))
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    else if(!k %in% read.csv(fn)$Species)
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    return(NULL)
  }
  # if the species has no points
  if(nrow(species.data)< 1){
    if(verbose) warning(paste0("Species '",k,"' has no points"))
    dn <- file.path(BASE_DIR,"problems")
    if(!dir.exists(dn))
      dir.create(dn, recursive =TRUE)
    fn <-file.path(dn,'Has_no_points.csv')
    if(!file.exists(fn))
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    else if(!k %in% read.csv(fn)$Species)
      write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
    return(NULL)
  }
  out = switch(speciesAlgo[j],

               "maxent" = UsefulPlants::run_maxentModel(),

               "rangebag" = UsefulPlants::run_geoModel(loc_dat = j,
                                                       species_name = k,
                                                       algorithm = geographic_distance_algorithm,
                                                       outputdir = file.path(BASE_DIR,"geoModels","geo_dist"),
                                                       output_res = geographic_grid_resolution,
                                                       verbose = FALSE),

               "point" = UsefulPlants::
  )
  if(class(out)=='try-error') cat(out)
  gc()
  raster::removeTMPFiles(h=1)
}, mc.cores=MC_CORES)
finished.at = Sys.time()
time.elapsed <- finished.at - started.at
cat('...End of computations...\n\n')
cat(paste0('> End of Run ',"'",RUN_NAME,"'",' on: ',Sys.time()),'\n\n');finished.at=proc.time()
cat(paste0('> Running time: ',as.numeric(time.elapsed),attr(time.elapsed,'units'),'.\n'))
cat(paste0('#',paste(rep('-',times=100),collapse="")))

