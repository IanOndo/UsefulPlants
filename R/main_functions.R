#'                                 Run geographic distance-based model workflow
#'
#' Partition both occurrence records and background into evaluation bins based on some spatial rules.
#'
#'
#' @param loc_dat A character string specifying the path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param algorithm A character string specifying the name of the algorithm to use.
#' @param outputdir A character string specifying the directory where to save the raster of interpolated values
#' @return An integer number. `0` in case of success otherwise return `-1`
#' @export
run_maxentModel <- function(loc_dat, outputdir, verbose=TRUE, maxent_settings=list(), ...){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  if(verbose) cat(sprintf(">>[[Starting the modelling of %s]]<<\n\n", k));flush.console()

  if(verbose) cat(sprintf("...[Setting up the parameters]...", k));flush.console()
  list_maxent_args <- list(loc_dat=loc_dat, outputdir=outputdir)

  extra_maxent_args <- list(...)
  if(length(extra_maxent_args)>0L){
    extra_maxent_args <- extra_maxent_args[names(extra_maxent_args) %in% c("k","coordHeaders","bg_masks","species_name","eval_metrics","varying_parameter_name")]
    if(length(extra_maxent_args)>0L)
      list_maxent_args <- append(list_maxent_args, extra_maxent_args)
  }
  env_dn <- file.path(outputdir,"modelsVariables")
  if(!dir.exists(env_dn))
    stop("Unable to locate directory: ", env_dn)

  list_maxent_args[['env_dat']] <- env_dn
  list_maxent_args <- append(list_maxent_args, maxent_settings)

  out_stat_dn <- file.path(outputdir,"modelsStats",k)
  if(!dir.exists(out_stat_dn)) dir.create(out_stat_dn, recursive=TRUE)
  if(verbose) cat("> ...[Training with block cross-validation]...");flush.console()
  eval_trained_model <- try(do.call(block_cv_maxent, list_maxent_args), silent=TRUE)


  if(!inherits(eval_trained_model,"try-error")){
    out_model_fn <- file.path(out_stat_dn, paste0(k,".rds"))
    if(!file.exists(out_model_fn))
      stop("unable to locate file: ", out_model_fn)
    if(verbose) cat("OK\n");flush.console()
    trained_model <- readRDS(out_model_fn)$model.object

    if(verbose) cat("> ...[Projecting]...");flush.console()
    out_pred_dn <- file.path(outputdir,"/modelsPredictions")
    if(!dir.exists(out_pred_dn)) dir.create(out_pred_dn, recursive=TRUE)
    projected_model <- try(project_geoModel(domain=, model=trained_model, output_directory=out_pred_dn, output_name=k, save_output=TRUE), silent=TRUE)
    if(!inherits(projected_model, "try-error")){
      if(verbose) cat("OK\n");flush.console()
      if(verbose) cat(sprintf(">>[[COMPLETED]]<<\n\n", k));flush.console()
      return(0)
    }else{
      if(verbose) cat("FAILED\n");flush.console()
      if(verbose) cat(sprintf(">>[[FAILURE]]<<\n\n", k));flush.console()
      if(verbose) cat(projected_model);flush.console()
      return(-1)
    }
  }

  if(verbose) cat(eval_trained_model)
  return(-1)
}

#'                                 Run geographic distance-based model workflow
#'
#' Uses only known occurrence locations to predict the likelihood of finding species in an area surrounding its known presence.
#' The likelihood of finding the species is thus inversely proportional to the distance to the nearest known presence record.
#'
#' @param loc_dat A character string specifying the path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param algorithm A character string specifying the name of the algorithm to use.
#' @param outputdir A character string specifying the directory where to save the raster of interpolated values
#' @param output_res A numeric value specifying the resolution of the output raster in decimal degree. Default is 0.1666667 i.e ~ 20km
#' @return An integer number. `0` in case of success otherwise return `-1`
#' @export
run_geoModel <- function(loc_dat, species_name=NULL, algorithm="idw", outputdir, output_res= 0.1666667, verbose=TRUE, ...){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  if(!is.numeric(output_res)|| is.na(output_res))
    stop("The output resolution must be numeric")

  if(is.null(species_name)){
    warning("species_name is null. Will assign absolute data-time as species name.")
    species_name=paste0("Unknown_", gsub(" ","_",Sys.time()))
  }

  if(verbose) cat(sprintf(">>[[Starting the modelling of %s]]<<\n\n", species_name));flush.console()

  if(verbose) cat("> ...[Training]...");flush.console()
  out_stat_dn <- file.path(outputdir,"modelsStats")
  if(!file.exists(out_stat_dn))
    dir.create(out_stat_dn, recursive=TRUE)
  eval_trained_model <- try(train_geoModel(loc_dat, algorithm=algorithm, outputdir=out_stat_dn, species_name=species_name),silent=TRUE)

  if(!inherits(eval_trained_model,"try-error")){
    out_model_fn <- file.path(out_stat_dn, paste0(species_name,".rds"))
    if(!file.exists(out_model_fn))
      stop("unable to locate file: ", out_model_fn)
    if(verbose) cat("OK\n");flush.console()
    trained_model <- readRDS(out_model_fn)$model.object

    if(verbose) cat("> ...[Preparing data for projection]...");flush.console()
    out_bg_fn <- file.path(outdir,"modelsBackground",paste0(species_name,".rds"))
    if(!file.exists(out_bg_fn)){
      cat("No region available for projecting the model")
      return(-1)
    }
    bg_domain <- readRDS(out_bg_fn)
    r <- raster::raster(bg_domain, res = output_res)
    projected_domain <- fasterize::fasterize(sf::st_cast(bg_domain,"MULTIPOLYGON"), r, fun="count")

    if(verbose) cat("> ...[Projecting]...");flush.console()
    out_pred_dn <- file.path(outputdir,"modelsPredictions")
    if(!dir.exists(out_pred_dn))
      dir.create(out_pred_dn, recursive=TRUE)

    projected_model <- try(project_geoModel(domain=projected_domain, model=trained_model, output_directory=out_pred_dn, output_name=species_name, save_output=TRUE), silent=TRUE)

    if(!inherits(projected_model, "try-error")){
      if(verbose){
        cat("OK\n");flush.console()
        cat(sprintf(">>[[COMPLETED]]<<\n\n"));flush.console()
      }
      return(0)
    }else{
      if(verbose){
        cat("FAILED\n");flush.console()
        cat(sprintf(">>[[FAILURE]]<<\n\n"));flush.console()
        cat(projected_model);flush.console()
      }
      return(-1)
    }
  }
  if(verbose){
    cat("FAILED\n");flush.console()
    cat(sprintf(">>[[FAILURE]]<<\n\n"));flush.console()
    cat(eval_trained_model);flush.console()
  }

  return(-1)
}

#'                                 Run geographic point-based model
#'
#' Uses known occurrence locations as its own.
#'
#'
#' @param loc_dat A character string specifying the path to a csv file with longitude and latitude coordinates of occurrence records.
#' @param algorithm A character string specifying the name of the algorithm to use.
#' @param outputdir A character string specifying the directory where to save the raster of interpolated values
#' @param output_res A numeric value specifying the resolution of the output raster in decimal degree. default is 0.1666667 i.e ~ 20km
#' @return An integer number. `0` in case of success otherwise return `-1`
#' @export
run_pointModel <- function(loc_dat, species_name=NULL, outputdir, output_res= 0.1666667, verbose=TRUE, ...){

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(missing(outputdir))
    stop("Output directory is missing.")

  if(!is.numeric(output_res)|| is.na(output_res))
    stop("The output resolution must be numeric")

  if(is.null(species_name)){
    warning("species_name is null. Will assign absolute data-time as species name.")
    species_name=paste0("Unknown_", gsub(" ","_",Sys.time()))
  }
  if(verbose) cat(sprintf(">>[[Starting the modelling of %s]]<<\n\n", species_name));flush.console()

  dat <- read.csv(loc_dat)
  #---------------------------------------------------------------------------
  #= 1. get coords (with a regular expression to look for longitude, latitude)
  #---------------------------------------------------------------------------
  x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
  y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]|[Bn][Ii][Nn]",x = names(dat),value=TRUE)[1]
  coordHeaders <- c(x_lon,y_lat)

  if(all(!coordHeaders %in% names(dat))){
    cat("Unable to extract coordinates headers from the data")
    return(-1)
  }

  sf_loc_data <- loc_dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))

  #---------------------------------------------------------------------------
  #= 2. fasterize
  #---------------------------------------------------------------------------
  bbox_domain <- make_bbox_domain(sf_loc_data)
  r <- raster::raster(bbox_domain, res = output_res)
  projected_domain <- fasterize::fasterize(sf_loc_data, r, fun="count")

}

