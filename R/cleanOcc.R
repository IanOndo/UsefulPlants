#'                                 Clean Species Occurrences
#'
#' Cleans species occurrences using CoordinateCleaner and TDWG region
#'
#' @param species_occ Either a directory or a csv file providing occurrence records to clean.
#' @param cutoff_year A numeric specifying the cutoff year beyond which records should be removed.
#' @param cutoff_coord_uncertainty A numeric specifying the cutoff degree of coordinates uncertainty (in m) beyond which records should be removed.
#' @param buffer_cen A numeric specifying the size of the buffer around each province or country centroids, where records should be ﬂagged as problematic. Default = 1 kilometre
#' @param buffer_cap A numeric specifying the size of the buffer around each capital coordinate (the centre of the city), where records should be ﬂagged as problematic. Default = 10 kilometres
#' @param buffer_inst A numeric specifying the size of the buffer around each institution, where records should be ﬂagged as problematic. Default = 1 kilometre
#' @param buffer_gbif_HQ A numeric specifying the size of the buffer around the gbif headquarters, where records should be ﬂagged as problematic. Default = 1 kilometre
#' @param output_dir A directory where to save species occurrence records data after cleaning.
#' @param use_TDWG A logical. Should the occurrences be cleaned by the TDWG regions. Default is \code{TRUE}.
#' @param species_id [OPTIONAL] A character string specifying the id of the species in ipni or Kew world checklist database. Only used if \code{use_TDWG=TRUE}.
#' @param mc_core A numeric integer specifying the number of cores to be used for parallel computing. Should be >1 for parallel processing.
#' @return None
#' @export
cleanOcc <- function(species_occ = system.file("extdata/UsefulPlants_workflow/Occ_dir/data_output", package='UsefulPlants'),
                     cutoff_year = 1945,
                     cutoff_coord_uncertainty = 20000,
                     buffer_cent = 1000,
                     buffer_cap = 10000,
                     buffer_inst = 1000,
                     buffer_gbif_HQ = 1000,
                     output_dir = system.file("extdata/Occ_dir/data_output", package='UsefulPlants'),
                     use_TDWG = TRUE,
                     species_id = NULL,
                     mc_cores	= NULL)
{
  #----------------------------------------
  #= 0.Load packages
  #----------------------------------------
  pkgs.to.load  <- c('data.table','CoordinateCleaner','parallel','stringr', 'dplyr', 'raster')
  pkgs.loaded   <- sapply(pkgs.to.load,require,character.only=TRUE)
  if(!any(pkgs.loaded)){
    warning(paste('Packages',paste(pkgs.to.load[!pkgs.loaded],collapse=', '),'failed to load'))
    stop("Try to re-install packages that failed to load")
  }
  #--------------------
  #= 1. Check inputs
  #--------------------
  if(is.null(species_occ))
    stop("Please provide valid occurrence data")

  file_flag= tryCatch(file.exists(species_occ), error=function(err) FALSE) && !tryCatch(dir.exists(species_occ), error=function(err) FALSE)
  dir_flag = tryCatch(dir.exists(species_occ), error=function(err) FALSE)
  data_flag = !file_flag & inherits(species_occ,c("data.frame","data.table"))

  if(!any(c(file_flag,dir_flag,data_flag)))
    stop('Please provide a valid directory, csv file or data.frame/data.object object with occurrence records')

  if(file_flag && length(readLines(species_occ))==0L)
    stop(paste("The file provided:",species_occ," is empty!"))

  if(!any(sapply(c(cutoff_year,cutoff_coord_uncertainty,buffer_cent,buffer_cap,buffer_inst,buffer_gbif_HQ),is.numeric)))
    stop('Arguments: cutoff_year, cutoff_coord_uncertainty, buffer_cent, buffer_cap, buffer_inst, buffer_gbif_HQ must be numeric ! \n')

  # ensure that the output directory provided exists
  if(!dir.exists(output_dir)){
    warning(paste("Output directory:", output_dir,"does not exist or could not be found, and will be set to the current working directory ./cleanedOcc"));flush.console()
    output_dir=file.path(getwd(), 'cleanedOcc')
    dir.create(output_dir, recursive=TRUE)
  }

  # check if cleaning by TDWG regions is possible
  if(use_TDWG & !requireNamespace("TDWG", quietly = TRUE)){
    warning("Cleaning by TDWG region is not possible because the package 'TDWG' is not installed.")
    use_TDWG = FALSE
  }

  if(use_TDWG & !is.null(species_id)){
    if(dir_flag){
      if(length(list.files(species_occ,recursive=TRUE)) != length(species_id))
        stop("Argument 'species id' must have the same length as the number of species to be cleaned.")
    }
    if(length(names(species_id))==0L)
      stop("Argument 'species_id' must be a named vector")

    if(data_flag)
      stop("Cleaning using species ids on input R object dataset is not available yet. Please provide another type of data (i.e. path to a csv file or a directory), or set species_id argument to NULL.")

    file_tmp <- tempfile("speciesIDs",fileext=".csv")
    write.csv(data.frame(SpeciesID = species_id, Species=names(species_id)), file_tmp, row.names=FALSE)
    species_id <- file_tmp
  }

  if(!is.null(mc_cores) & !is.numeric(mc_cores))
    stop("Argument 'mc_cores' must be a numeric integer")

  if(is.null(mc_cores)) mc_cores = detectCores()-1

  #---------------------
  #= 2. Set parameters
  #---------------------
  if(data_flag){
    file_tmp <- tempfile("cleanOcc",fileext=".csv")
    write.csv(species_occ, file_tmp, row.names=FALSE)
    species_occ <- file_tmp
  }

  params = c(

    'species_occ'	= species_occ,		# path to occurrence records data

    'data_type' = (1:3)[c(file_flag,dir_flag,data_flag)],

    'cutoff_year' = cutoff_year,

    'cutoff_coord_uncertainty' = cutoff_coord_uncertainty,

    'buffer_cent' = buffer_cent,

    'buffer_cap' = buffer_cap,

    'buffer_inst' = buffer_inst,

    'buffer_gbif_HQ'= buffer_gbif_HQ,

    'output_dir' 	  = output_dir,

    'species_id' = species_id,

    'use_TDWG'  = use_TDWG,

    'mc_cores' 			= round(mc_cores)  # number of cores to use
  )

  #====================================================================
  #== 3. Run species search
  #====================================================================
  commandArgs <- function(...) paste0(names(params),'=',params)
  assign('commandArgs',commandArgs,envir=.GlobalEnv)
  source(system.file("extdata/UsefulPlants_workflow/occCleaning.R", package='UsefulPlants'))
}
