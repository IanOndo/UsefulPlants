#' Sort species by number of occurrence records retrieve from occurrence data files.
#'
#' Assign species occurrence data to different folders based on minimum number of occurrence records
#'
#' @param base_directory A character string specifying the directory where to store the outputs of the workflow
#' @param species_directory A character string specifying the directory where are stored the species occurrence records.
#' @param env_directory A character string specifying the directory where are stored the environmental layers.
#' @param min_occ_envmodel A numeric integer specifying the minimum number of occurrence records required to assign a species to the environmental model directory
#' @param min_occ_geomodel A numeric integer specifying the minimum number of occurrence records required to assign a species to the geographic model directory
#' @param nchunk A numeric integer specifying the number of chunks in which the vector of species names must be divided. Must be > 1 to enable parallel processing.
#' @param ncores A numeric integer specifying the number of cores (or CPU's) to use in case of parallel processing. Default is `number of cores available` - 1.
#' @export
sortByInitNumRecords <- function(base_directory, species_directory, env_directory, algorithm, min_occ_envmodel=10, min_occ_geomodel=3, nchunk=20, ncores=NULL, verbose=FALSE){

  cat("\n")
  if(verbose) cat(">...checking input directories...")
  # ensure all arguments are provided
  if(missing(base_directory)){
    cat("FAILED\n")
    stop("The base directory is missing.")
  }
  if(missing(species_directory)){
    cat("FAILED\n")
    stop("The directory with species occurrence data is missing.")
  }

  # ensure that the directories provided exist
  if(!dir.exists(base_directory)){
    cat("FAILED\n")
    stop("Base directory:", base_directory,"does not exist or could not be found.");flush.console()
  }

  if(verbose) cat("OK\n")
  if(verbose) cat("> ...ensure species directory has data...")
  # ensure species directory has data
  list.species <- list.files(species_directory, pattern="\\.csv$", full.names=TRUE, recursive=TRUE) # select species
  if(length(list.species)==0L)
    stop(paste0("\nUnable to find csv files in directory:",species_directory,". Please provide a directory with csv files."))
  empty.string <- nchar(list.species)==0L
  if(any(empty.string)){
    if(verbose)
      warning(paste("\nRemoving",sum(empty.string),"file(s) with no species names."))
    list.species <- list.species[!empty.string]
  }
  if(verbose) cat("OK\n")

  if(!is.numeric(nchunk)){
    stop("nchunk must be numeric")
  }
  if(nchunk<=0){
    stop("nchunk must be > 0")
  }

  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min(ncores, parallel::detectCores()) # number of cores to use

  if(verbose) cat("> ...split species files list in ",nchunk,"chunks and sort species by number of occurrence records...")

  species_split_list <- chunk2(list.species, nchunk)

  # get the number of records by species
  if(nchunk>1)
    doParallel::registerDoParallel(ncores)
  else foreach::registerDoSEQ()

  num_records_by_species <- foreach::foreach(f=1:nchunk, .combine='c') %dopar% {

    sapply(species_split_list[[f]], function(x){

      dat <- read.csv(x)
      #---------------------------------------------------------------------------
      #= 1. get coords (with a regular expression to look for longitude, latitude)
      #---------------------------------------------------------------------------
      x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
      y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
      coordHeaders <- c(x_lon,y_lat)

      if(all(!coordHeaders %in% names(dat))){
        return(0)
      }
      #---------------------------------------------------------------------------
      #= 2. remove duplicated coords
      #---------------------------------------------------------------------------
      xtmp  <- dat[, coordHeaders]
      nrow(xtmp[!duplicated(xtmp),])

      }, USE.NAMES=FALSE)
  }

  if(nchunk>1)
    doParallel::stopImplicitCluster()

  # assign species to model directories
  has_bad_headers <- num_records_by_species==0 # if the data file has bad headers
  if(any(has_bad_headers)){
    dn <- file.path(base_directory,"problems/has_bad_headers")
    if(!dir.exists(dn))
      dir.create(dn, recursive =TRUE)
    file.copy(from=list.species[has_bad_headers],
              to=file.path(dn, basename(list.species[has_bad_headers])))
  }
  are_environmental_model <- num_records_by_species >= min_occ_envmodel
  are_point_model <- num_records_by_species < min_occ_geomodel
  are_geo_model <- !are_environmental_model & !are_point_model

  if(any(are_environmental_model)){
      dirs_created <- init_model_directory(base_directory, type="environmental", algorithm="maxent")
      if(!dirs_created){
        file.copy(from=list.species[are_environmental_model],
                  to=file.path(base_directory,"envModels","maxent","occurrences",basename(list.species[are_environmental_model])))

        # copy environmental variables
        if(!missing(env_directory) && !inherits(try(read_layers(env_directory),silent=TRUE),'try-error')){
          files.copied <- file.copy(from=list.files(env_directory, pattern="\\.tif$"), to=file.path(directory,"envModels","modelsVariables", basename(list.files(env_directory, pattern="\\.tif$"))))
          if(any(!files.copied))
            stop("The following files could not be copied: ", list.files(env_directory, pattern="\\.tif$")[!files.copied])
        }
      }
  }

  if(any(are_geo_model)){
    dirs_created <- init_model_directory(base_directory, type="geographic", algorithm="geo_dist")
    if(!dirs_created)
      file.copy(from=list.species[are_geo_model],
                to=file.path(base_directory,"geoModels","geo_dist","occurrences",basename(list.species[are_geo_model])))
  }

  if(any(are_point_model)){
    dirs_created <- init_model_directory(base_directory, type="point", algorithm="point")
    if(!dirs_created)
      file.copy(from=list.species[are_point_model],
                to=file.path(base_directory,"pointModels","point","occurrences",basename(list.species[are_point_model])))
  }

}




#' helper function that split a vector x in n chunks of approximatively equal size
#' @export
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
