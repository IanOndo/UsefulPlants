#'                                 Query Species Occurrences from GBIF
#'
#' Queries species occurrences from GBIF
#'
#' @param species_name A single species name
#' @param gbif_login An object of class GBIFLogin to log in to GBIF to begin the download
#' @return A data.table object with occurrence data
#' @export
queryGBIF <- function(species_name, gbif_login, gbif_download_dir=NULL, rank='species', kingdom='plants', status_ping=5, verbose=TRUE){

  if(missing(species_name) | !is.character(species_name) | species_name==""){
    stop("Please provide a valid species name")
  }

  if(missing(gbif_login)){
    if(all(Sys.getenv(c('GBIF_USER','GBIF_PWD','GBIF_EMAIL'))==""))
      stop("Please provide a GBIFLogin object with your GBIF credentials")
    else
      gbif_login <-occCite::GBIFLoginManager(user = Sys.getenv('GBIF_USER'),
                                             email = Sys.getenv('GBIF_EMAIL'),
                                             pwd = Sys.getenv('GBIF_PWD'))
  }

  if(is.null(gbif_download_dir)){
    warning('Occurrence data will be downloaded in your current working directory.')
    gbif_download_dir = getwd()
  }

  #----------------------------------------------
  #=1. obtain taxon key ID for the target species
  #----------------------------------------------
  res <- rgbif::name_backbone(name = species_name, rank=rank, kingdom=kingdom)

  if(is.null(res$speciesKey)){
    res <- rgbif::name_suggest(q = species_name, rank=rank)
    if(!is.null(res$key[1])){
      if(res$canonicalName[1]==species_name)
        taxonkey <- as.numeric(res$key[1])
      else
        return(NULL)
    }else
      return(NULL)
  }else{
    taxonkey <- as.numeric(res$speciesKey)
  }


  #-------------------------------------------------------
  #=2. obtain species occurrence records from taxon key ID
  #-------------------------------------------------------
  occ_query <- try(rgbif::occ_download(paste0("taxonKey =", taxonkey),
                    "hasCoordinate = true",
                    format= "SIMPLE_CSV",
                    user = gbif_login@'username',
                    pwd = gbif_login@'pwd',
                    email = gbif_login@'email'),silent=TRUE)

  if(inherits(occ_query,'error')){
    cat('\n\n')
    cat(paste('<< Unable to download occurrence records for', species_name,'. Aborting...>>'))
    return(NULL)
  }

  #---------------------
  #=2a. obtain metadata
  #---------------------
  # # https://stackoverflow.com/questions/55847360/wait-for-rgbif-download-to-complete-before-proceeding
  is_running <- TRUE
  while(is_running){ # check until ready i.e. status is "SUCCEEDED"
    meta <- rgbif::occ_download_meta(occ_query)
    status <- meta$status
    is_running <- !status %in% c("SUCCEEDED", "KILLED")
    Sys.sleep(status_ping) # sleep between pings
  }
  # occ_file <- rgbif::occ_download_queue(.list = occ_query, status_ping=status_ping)


  #---------------------------------------------
  #=2b. download data from GBIF to local machine
  #---------------------------------------------
  occ_data <- rgbif::occ_download_get(occ_query, path = gbif_download_dir, overwrite=TRUE)  # if too large to download into R- download manually through GBIF or split into parts above, if
  if(attr(occ_data,'size')==0L) return(NULL)                                        # if manual download, manual import also needed, see rgbifFORMAT downloads for manual import and formatting
  #---------------------
  #=2c. import data
  #---------------------
  gbifDATA<- rgbif::occ_download_import(occ_data, header=TRUE, showProgress=FALSE, na.strings=c("",NA),
                                            select=c("species",
                                                     "taxonRank",
                                                     "infraspecificEpithet",
                                                     "decimalLongitude",
                                                     "decimalLatitude",
                                                     "countryCode",
                                                     "coordinateUncertaintyInMeters",
                                                     "year",
                                                     "gbifID",
                                                     "basisOfRecord",
                                                     "institutionCode",
                                                     "establishmentMeans",
                                                     "individualCount"))
  if(!inherits(gbifDATA,'error')){

    # set to data.table object
    data.table::setDT(gbifDATA)

    #--------------------
    #=3. format the data
    #--------------------
    gbifDATA<-na.omit(gbifDATA, cols= c("decimalLatitude", "decimalLongitude"))
    # format the data:
    gbifDATA[, `:=`(taxonRank = ifelse(taxonRank %in% c("SPECIES","GENUS"), NA, ifelse(taxonRank=="FORM", "f.", ifelse(taxonRank=="SUBSPECIES", "subsp.", ifelse(taxonRank=="VARIETY","var.", taxonRank)))),
                    countryCode = countrycode::countrycode(countryCode, origin =  'iso2c', destination = 'iso3c', nomatch = NA),
                    establishmentMeans = ifelse(establishmentMeans=="INTRODUCED","Introduced",ifelse(establishmentMeans== "NATIVE", "Native", establishmentMeans)))]
    # create 'fullname', 'is_cultivated_observation' and 'sourceID' columns
    gbifDATA[,`:=`(fullname = ifelse(is.na(infraspecificEpithet), paste(species), paste(species, taxonRank, infraspecificEpithet)),
                                                   is_cultivated_observation = NA,
                                                   sourceID = 'GBIF')]
    # delete taxonRank and infraspecificEpithet columns
    gbifDATA[, `:=`(taxonRank = NULL,
                    infraspecificEpithet = NULL)]

    # remove fossil records
    gbifDATA <- gbifDATA[basisOfRecord!="FOSSIL_SPECIMEN"]

    # set column order
    if(!exists('colNames'))
      colNames = c("species","fullname","decimalLongitude","decimalLatitude","countryCode","coordinateUncertaintyInMeters","year","individualCount","gbifID","basisOfRecord","institutionCode","establishmentMeans",
                   "is_cultivated_observation","sourceID")
    data.table::setcolorder(gbifDATA, colNames)
    # set the key to the species column to enable fast binary search
    data.table::setkey(gbifDATA, 'species')

  }else{
    cat('\n\n')
    cat('...<< Unable to import data from gbif database >>...')
    return(NULL)
  }

  if(verbose) cat('...DONE...\n\n')
  return(gbifDATA)
}

