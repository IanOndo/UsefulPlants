#'                                 Query Species Occurrences from BIEN
#'
#' Queries species occurrences from BIEN
#'
#' @param species_name A single species name
#' @return A data.table object with occurrence data
#' @export
queryBIEN <- function(species_name, verbose=TRUE){

  if(missing(species_name) | !is.character(species_name) | species_name==""){
    stop("Please provide a valid species name")
  }

  #-------------------------------------
  #=1. obtain species occurrence records
  #-------------------------------------
  bienDATA0 <- try(BIEN::BIEN_occurrence_species(species = species_name,
                                                 cultivated = T,
                                                 only.new.world = F,
                                                 all.taxonomy = T,
                                                 native.status = T,
                                                 observation.type = T,
                                                 political.boundaries = T,
                                                 natives.only = F,
                                                 collection.info = T),silent=TRUE)

  if(inherits(bienDATA0,'try-error') || nrow(bienDATA0)==0L){
    cat('\n\n')
    cat(paste('<< Unable to download occurrence records for', species_name,'. Aborting...>>'))
    return(NULL)
  }

  # set to data.table object
  data.table::setDT(bienDATA0)

  #--------------------
  #=2. format the data
  #--------------------
  bienDATA <- subset(bienDATA0, select=c("name_matched", "scrubbed_species_binomial", "longitude", "latitude", "country", "date_collected",
                                        "observation_type", "catalog_number", "is_cultivated_observation", "is_introduced"))
  # change the column names
  data.table::setnames(bienDATA, c("fullname", "species", "decimalLongitude", "decimalLatitude", "countryCode", "year",
                                   "basisOfRecord", "gbifID", "is_cultivated_observation", "establishmentMeans"))
  # format the data:
  if(verbose)
    cat('...Converting date into year...Assigning country code...Converting is_cultivated_observation code...Converting introduced (establishmentMeans) code...\n')

  bienDATA[, `:=`(year = as.numeric(format(as.Date(year,"%Y%m%d"),"%Y")),
                  countryCode = countrycode::countrycode(countryCode, origin =  'country.name', destination = 'iso3c', nomatch = NA),
                  is_cultivated_observation = ifelse(is_cultivated_observation==1,"Yes",ifelse(is_cultivated_observation== 0, "No", is_cultivated_observation)),
                  establishmentMeans = ifelse(establishmentMeans==1,"Introduced",ifelse(establishmentMeans== 0, "Native", establishmentMeans)))]
  if(verbose) cat('...Adding matched columns and sourceID...\n')
  bienDATA[,`:=`(coordinateUncertaintyInMeters = NA,
                 individualCount = NA,
                 institutionCode = NA,
                 sourceID = 'BIEN')]
  # set column order
  if(!exists('colNames'))
    colNames = c("species","fullname","decimalLongitude","decimalLatitude","countryCode","coordinateUncertaintyInMeters","year","individualCount","gbifID","basisOfRecord","institutionCode","establishmentMeans",
                 "is_cultivated_observation","sourceID")
  data.table::setcolorder(bienDATA, colNames)
  # set the key to the species column to enable fast binary search
  data.table::setkey(bienDATA, 'species')

  if(verbose) cat('...DONE...\n\n')
  return(data.table::copy(bienDATA))

}
