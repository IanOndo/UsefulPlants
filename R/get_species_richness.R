#'                                 Compute the species richness of geographic units
#'
#' Calculates the species richness at a geographic level
#'
#' @param species_name A character vector of species names
#' @param species_id A character vector of ids either from IPNI or the World Check List [Recommended if available].
#' @param status Optional. A character string specifying the type of geographic region where to search for the species. Can be : 'native', 'introduced' or 'both'. Default is 'native'.
#' @param level A numeric integer between 1 and 4. The unit level of the TWDG at which the species richness must be calculated. Default is level 2.
#' @param use_name_matching A logical. Should the species names be matched with the Kew checklist (i.e. try to guess actual accepted name) ? Default is FALSE.
#' @param do.parallel A logical. Should the computation be parallelized ? Default is TRUE
#' @param ncores A numeric integer specifying the number of cores to use for parallel computation. Ignored if \code{do.parallel} is set to FALSE.
#' @param verbose A logical. Should the name of the species processed be printed on the screen ? Default is TRUE.
#' @return A SpatialPolygonsDataFrame with species richness values stored in the data.
#' @export
get_species_richness <- function(species_name=NULL, species_id=NULL, status='native', level=2, use_name_matching = FALSE, do.parallel=TRUE, ncores=NULL, verbose=TRUE){

  # check if TDWG package is available
  if(!requireNamespace("TDWG", quietly = TRUE)){
    stop("Computing species richness is not possible because the package 'TDWG' is not installed.")
  }

  if(is.null(species_name) & is.null(species_id))
    stop("Please provide at least a species name or a species id vector")

  if(!is.null(species_name) & !is.null(species_id)){
    if(length(species_name)!=length(species_id))
      warning("Arguments 'species_name' and 'species_id' have different lengths.")
  }

  if(!is.null(species_name) && (!all(is.character(species_name)) | all(species_name=="")) ){
    stop("Please provide a valid species name vector")
  }

  if(!is.null(species_id) && (!all(is.character(species_id)) | all(species_id=="")) ){
    stop("Please provide a valid species id vector")
  }

  if(!inherits(level,"numeric") || (level<1 | level>4)){
    cat('failed\n')
    stop("Argument 'initial level' must be an integer number between 1 and 4")
  }

  if(do.parallel & (length(species_name)>1 | length(species_id)>1)){
    # Parallel processing
    ncores = if(is.null(ncores)) parallel::detectCores()-1 else min(ncores,parallel::detectCores()) # number of cores to use
    doParallel::registerDoParallel(ncores)
  }else{
    # Sequential processing
    foreach::registerDoSEQ()
  }

  if(!is.null(species_id))
    species_vector = species_id
  else
    species_vector = species_name

  output_geo_codes <- foreach::foreach(k=1:length(species_vector), .combine='c', .packages=c("TDWG","data.table","sp")) %dopar% {

    # if(use_name_matching && data.table::key(TDWG:::kew_checklist)!="full_name_without_family") data.table::setkey(TDWG:::kew_checklist,"full_name_without_family")
    # if(!use_name_matching && data.table::key(TDWG:::kew_checklist)!="acc_full_name_without_family") data.table::setkey(TDWG:::kew_checklist,"acc_full_name_without_family")

    if(verbose){
      if(!is.null(species_name)){
         cat(sprintf("Species name: %s\n\n", species_vector[k]));flush.console()
      }
    }

    if(!is.null(species_id))
      TDWG::get_geocode(species_id=species_vector[k], status=status, level=level, use_name_matching=use_name_matching)
    else
      TDWG::get_geocode(species_name=species_vector[k], status=status, level=level, use_name_matching=use_name_matching)
  }

  if(do.parallel & (length(species_name)>1 | length(species_id)>1))
    doParallel::stopImplicitCluster()

  table_output <- table(output_geo_codes)

  switch(level,

         {
           output_spdf <- TDWG:::tdwg_level1
           idx <- match(TDWG:::tdwg_level1[['LEVEL1_COD']], as.numeric(names(table_output)), nomatch=0)
         },
         {
           output_spdf <- TDWG:::tdwg_level2
           idx <- match(TDWG:::tdwg_level2[['LEVEL2_COD']], as.numeric(names(table_output)), nomatch=0)
         },
         {
           output_spdf <- TDWG:::tdwg_level3
           idx <- match(TDWG:::tdwg_level3[['LEVEL3_COD']], names(table_output), nomatch=0)
         },
         {
           output_spdf <- TDWG:::tdwg_level4
           idx <- match(TDWG:::tdwg_level4[['area_code_l3']], as.numeric(names(table_output)), nomatch=0)
         }

  )

  #assign species richness values to geographic units
  output_spdf@data <- data.frame(species_richness=rep(NA,length(output_spdf)))
  output_spdf@data[['species_richness']]<- replace(idx, idx>0, table_output[idx])

  #return SpatialPolygonsDataFrame
  return(output_spdf)
}

