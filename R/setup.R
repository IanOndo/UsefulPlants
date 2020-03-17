#' Initialize SDM output directories
#'
#' Creates directories that will contain the outputs of the workflow
#'
#' @param base_directory A character string specifying the directory where to store the outputs of the workflow
#' @param species_directory A character string specifying the directory where are stored the species occurrence records.
#' @param env_directory A character string specifying the directory where are stored the environmental layers.
#' @param verbose A logical. Default is `TRUE`.
#' @return None
#' @export
init_main_directories <- function(base_directory, species_directory, env_directory, verbose=TRUE, force_sorting=FALSE, ...){

  if(verbose) cat("> ...checking input directories...")
  # ensure all arguments are provided
  if(missing(base_directory)){
    cat("FAILED\n")
    stop("The base directory is missing.")
  }
  if(missing(species_directory)){
    cat("FAILED\n")
    stop("The directory with species occurrence data is missing.")
  }
  if(missing(env_directory)){
    cat("\n")
    warning("The directory with environmental layers is missing. You will not be able to run environmnental niche models, but only geographic models for your species.")
  }

  # ensure that the directories provided exist
  if(!dir.exists(base_directory)){
    cat("FAILED\n")
    stop("Base directory:", base_directory,"does not exist or could not be found.");flush.console()
  }
  if(!dir.exists(species_directory)){
    cat("FAILED\n")
    stop("The directory:", species_directory,"does not exist or could not be found.");flush.console()
  }
  if(!missing(env_directory) && !dir.exists(env_directory)){
    cat("FAILED\n")
    stop("The directory:", env_directory,"does not exist or could not be found.");flush.console()
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

  if(!missing(env_directory)){
    if(verbose) cat("> ...ensure environment directory has data...")
    # ensure environment directory has data
    env_layers <- try(read_layers(env_directory),silent=TRUE)
    if(!missing(env_directory) && inherits(env_layers,'try-error')){
      stop(cat(env_layers))
    }
    if(verbose) cat("OK\n")
  }

  # check up the base directory
  has_data <- FALSE
  files_to_remove <- list.files(base_directory, pattern="\\.csv$", recursive=TRUE, full.names=TRUE)
  has_data <-length(files_to_remove)>0L
  if(has_data){
    warning("The base directory already contain data.") #
  }

  list_args <- list(base_directory, species_directory)
  if(!missing(env_directory))
    list_args <- append(list_args, env_directory)

  dot.args <- list(...)
  if(length(dot.args>0L))
    list_args <- append(list_args, dot.args[names(dot.args) %in% c("min_occ_envmodel","min_occ_geomodel", "nchunk", "ncores")])

  if(!has_data | force_sorting){

    if(force_sorting & has_data){
      if(verbose) cat("> ...Removing previous files...")
      files.removed <- file.remove(files_to_remove, recursive=TRUE)
      if(any(!files.removed))
        cat("FAILED\n")
      else
        cat("OK\n")
    }

    if(verbose) cat("> ...Sorting species by number of occurrence records found in the file...")
    do.call(sortByInitNumRecords, list_args)
    if(verbose) cat("OK\n")
  }
  else
    if(verbose) cat("> ...No sorting required...\n")

}

#' Initialize model directory
#'
#' Creates sub-directories that will contain the outputs of the type of model specified
#'
#' @param directory A character string specifying the directory where to store the outputs of the model
#' @param type A character string specifying the type of model directory to create
#' @param algorithm A character string specifying the name of the algorithm used as model
#' @export
init_model_directory <- function(directory, type="point", algorithm="point", verbose=FALSE){

  if(verbose) cat("> ...checking input directory...")
  # ensure all arguments are provided
  if(missing(directory)){
    cat("FAILED\n")
    stop("The directory is missing.")
  }

  # ensure that the directory provided exist
  if(!dir.exists(directory)){
    cat("FAILED\n")
    cat("Directory:", directory,"does not exist or could not be found.");flush.console()
    return(-1)
  }
  if(verbose) cat("OK\n")

  # list of directories to create
  dirs_to_create <- switch(type,

   "environmental" = file.path(directory,"envModels",algorithm,
                    c("occurrences",
                      "modelsVariables",
                      "modelsStats",
                      "modelsFigures",
                      "modelsPredictions",
                      "modelsUncertainties",
                      "modelsMaps")),

   "geographic" = file.path(directory,"geoModels",algorithm,
                            c("occurrences",
                              "modelsStats",
                              "modelsFigures",
                              "modelsPredictions",
                              "modelsUncertainties",
                              "modelsMaps")),

   "point" = file.path(directory,"pointModels",algorithm,
                       c("occurrences",
                         "modelsFigures",
                         "modelsPredictions",
                         "modelsMaps"))

  )

  dirs_created <- sapply(dirs_to_create, function(dn) if(!dir.exists(dn)) dir.create(dn,recursive=TRUE))
  if(any(!dirs_created)){
    cat("The following directories could not be created: ", dirs_to_create[!dirs_created])
    return(-1)
  }

  return(0)
}


