#'                                 Query Species Occurrences from GBIF
#'
#' Queries species occurrences from GBIF
#'
#' @param species_name A vector of species names
#' @param user A vector of type character specifying a GBIF username.
#' @param email A vector of type character specifying the email associated with a GBIF username.
#' @param pwd A vector of type character containing the user's password for logging in to GBIF.
#' @param rank the given rank at which to search. Default is 'species'.
#' @param kingdom the kingdom in which to search. Default is 'Plantae'.
#' @return A data.table object with occurrence data
#' @export
queryGBIF <- function(species_name, user=NULL, email=NULL, pwd=NULL, gbif_download_dir=NULL, rank='species', kingdom='Plantae', status_ping=3, time_out=300, index=NULL, check_output="", verbose=TRUE){

  if(missing(species_name) | any(!is.character(species_name)) | any(species_name=="")){
    stop("Please provide a valid species name")
  }

  if(is.null(user) | is.null(email) |is.null(pwd)){
    if(all(Sys.getenv(c('GBIF_USER','GBIF_PWD','GBIF_EMAIL'))==""))
      stop("Please provide all your GBIF credentials")
    else{
      user = Sys.getenv('GBIF_USER')
      email = Sys.getenv('GBIF_EMAIL')
      pwd = Sys.getenv('GBIF_PWD')
    }
  }

  if(any(is.invalid <- !sapply(c(user, email, pwd),is.character))){
    cat('Input', c('user', 'email', 'pwd')[is.invalid],'invalid. Please specify a character string',c('user', 'email', 'pwd')[is.invalid],'\n')
    return(NULL)
  }

  if(is.null(gbif_download_dir)){
    warning('Occurrence data will be downloaded in your current working directory.')
    gbif_download_dir = getwd()
  }

  #----------------------------------------------
  #=0. check processing step
  #----------------------------------------------
  do.query <- TRUE

  call.fun 	  <- match.call(expand.dots=FALSE)
  tmp.args    <- c("",'species_name', 'user','email', 'pwd', 'gbif_download_dir', 'rank', 'kingdom')
  call.tmp    <- call.fun[match(tmp.args, names(call.fun),nomatch=0)]

  gbif_query_fn <- file.path(dirname(raster::rasterTmpFile()),"gbif_query.rds")
  if(file.exists(gbif_query_fn)){
    gbif_query  	<- readRDS(gbif_query_fn)
    if(identical(call.tmp,gbif_query[[1]])){
      occ_query <- gbif_query[[2]]
      do.query <- FALSE
    }

  }

  if(do.query){

    #----------------------------------------------
    #=1. obtain taxon key ID for the target species
    #----------------------------------------------
    taxon_keys <- suppressMessages(get_gbiftaxonkey(species_name=species_name, rank=rank, kingdom=kingdom, check.output=check_output, index=index))
    taxon_keys <- na.omit(taxon_keys)

    if(length(taxon_keys)==0L)
      return(NULL)

    #-------------------------------------------------------
    #=2. obtain species occurrence records from taxon key ID
    #-------------------------------------------------------
    is_pending = TRUE
    time_elapsed = 0
    while(is_pending){
      num_downloads <- rgbif::occ_download_list()[['meta']][['count']]
      pending_keys <- subset(rgbif::occ_download_list(limit=num_downloads)[['results']], status %in% c("PREPARING","RUNNING"), select='key',drop=T)
      Sys.sleep(status_ping)
      time_elapsed = time_elapsed + status_ping
      is_pending = length(pending_keys)>3 & time_elapsed < time_out
    }
    if(length(pending_keys)>3L)
      suppressMessages(job_cancelled <- sapply(pending_keys[1:(length(pending_keys)-3)], rgbif::occ_download_cancel))

    occ_query <- try(rgbif::occ_download(
                      rgbif::pred_in("taxonKey", taxon_keys),
                      rgbif::pred("hasCoordinate",TRUE),
                      format= "SIMPLE_CSV",
                      user = user,
                      pwd = pwd,
                      email = email),silent=TRUE)

    if(inherits(occ_query,'try-error')){
      cat('\n')
      cat(paste('<< Unable to obtain occurrence records for:', paste(head(species_name),collapse=", "),'.... [Aborting]...>>\n'))
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
    if(status =='KILLED'){
      cat('\n')
      cat(paste('<< Unable to obtain metadata for:', paste(head(species_name),collapse=", "),'... [Aborting]...>>\n'))
      return(NULL)
    }

    # save query
    saveRDS(list(call.tmp, occ_query), file=gbif_query_fn)
  }

  #---------------------------------------------
  #=2b. download data from GBIF to local machine
  #---------------------------------------------
  occ_data <- try(rgbif::occ_download_get(occ_query, path = gbif_download_dir, overwrite=TRUE), silent=TRUE)  # if too large to download into R- download manually through GBIF or split into parts above, if

  if(inherits(occ_data,'try-error') ){#|| attr(occ_data,'size')==0L
    cat('\n')
    if(verbose) cat('...<< Unable to download data from GBIF to local machine >>...\n')
    return(NULL)
  }                                      # if manual download, manual import also needed, see rgbifFORMAT downloads for manual import and formatting

  #---------------------
  #=2c. import data
  #---------------------
  gbifDATA<- try(rgbif::occ_download_import(occ_data, header=TRUE, showProgress=FALSE, na.strings=c("",NA), fill=FALSE, quote = "",
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
                                                     "individualCount")), silent=TRUE)
  if(!inherits(gbifDATA,'try-error')){

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
    gbifDATA[, `:=`(fullname = ifelse(is.na(infraspecificEpithet), paste(species), paste(species, taxonRank, infraspecificEpithet)),
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
    cat('\n')
    cat('...<< Unable to import data from gbif database >>...\n')
    return(NULL)
  }

  if(verbose) cat('...DONE...\n\n')
  return(gbifDATA)
}

#'                                 Extract Taxon key from GBIF
#'
#' Retrieves the taxon key value of a species from GBIF database
#' Modified from https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/
#' @param species_name A vector of species names
#' @param rank the given rank at which to search. Default is 'species'.
#' @param kingdom the kingdom in which to search. Default is 'Plantae'.
#' @return A numeric value corresponding to the taxon key in the GBIF database
#' @export
get_gbiftaxonkey <- function(species_name, rank='species', kingdom='Plantae', save.outputs=TRUE, index=NULL, check.output=""){

  result <- taxize::get_gbifid_(sciname = species_name,  method="backbone")

  if(all(lengths(result)==0)) return(NULL)

  #-------------------------------
  #= 0. Keep original species name
  #     Combine data into one
  #-------------------------------
  taxon_key <- result %>%

    purrr::imap(~ .x %>% dplyr::mutate(original_sciname = .y)) %>%

    dplyr::bind_rows()

  #-------------------------------
  #= 1. Filter by kingdom and rank
  #-------------------------------
  taxon_key <- taxon_key %>%

    dplyr::filter(rank == rank & kingdom == kingdom)

  if(nrow(taxon_key)==0) return(NULL)

  #-------------------------------
  #= 2. Keep exact matches only
  #-------------------------------
  taxon_key <- taxon_key %>%

    dplyr::filter(matchtype == "EXACT")

  if(nrow(taxon_key)==0) return(NULL)

  #-------------------------------
  #= 3. Save outputs if required
  #-------------------------------
  if(save.outputs & dir.exists(check.output)){
    index = ifelse(!is.null(index), index, round(runif(1,min=1,max=10000)))
    write.csv(taxon_key, file=file.path(check.output,paste0('check_matched_names_',index,'.csv')), row.names=FALSE)
  }

  #-------------------------------
  #= 4. Extract the taxon key
  #-------------------------------
  taxon_key <- taxon_key %>%

    dplyr::pull("specieskey")

  # if(is.null(res$speciesKey)){
  #   res <- rgbif::name_suggest(q = species_name, rank=rank)
  #   if(!is.null(res$key[1])){
  #     if(res$canonicalName[1]==species_name)
  #       taxonkey <- as.numeric(res$key[1])
  #     else
  #       return(NULL)
  #   }else
  #     return(NULL)
  # }else{
  #   taxonkey <- as.numeric(res$speciesKey)
  # }

  return(taxon_key)
}
