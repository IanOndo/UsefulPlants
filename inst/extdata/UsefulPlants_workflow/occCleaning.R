## If the OS is Windows, set mclapply to the hackish version. Otherwise, leave the definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {UsefulPlants::mclapply.hack},
                    Linux   = {parallel::mclapply},
                    Darwin  = {parallel::mclapply})

#====================================================================
#== 0. Get parameters from command-line
#====================================================================
SPECIES_ID = ""

args = commandArgs(trailingOnly=TRUE)

for(arg in args) {

  arg_name = stringr::str_extract(arg,".+?(?=\\=)")

  switch(arg_name,
         'species_occ'  = {SPECIES_OCC = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'data_type'  = {DATA_TYPE = as.integer(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'cutoff_year' = {CUTOFF_YEAR = as.integer(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'cutoff_coord_uncertainty' = {CUTOFF_COORD_UNCERTAINTY = as.numeric(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'buffer_cent' = {BUFFER_CENT = as.numeric(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'buffer_cap' = {BUFFER_CAP = as.numeric(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'buffer_inst' = {BUFFER_INST = as.numeric(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'buffer_gbif_HQ'= {BUFFER_GBIF_HQ = as.numeric(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'output_dir'   = {OUTPUT_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'use_TDWG'   = {USE_TDWG = as.logical(gsub(".+?(?<=\\=)","",arg,perl=TRUE))},

         'species_id'   = {SPECIES_ID = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'mc_cores' 	= {MC_CORES = as.integer(gsub(".+?(?<=\\=)","",arg,perl=TRUE))}
  )
}

cat(paste0('#',paste(rep('-',times=100),collapse="")))
cat('\n')
cat(paste0("> Run started on: ",Sys.time()));
cat('\n')
cat('#====================================================================\n')
cat('#== 1. Check input parameters\n')
cat('#====================================================================\n\n')
started.at <- Sys.time()

switch(DATA_TYPE,

       # file
       {
        if(!grepl("\\.csv$", SPECIES_OCC))
          stop(paste0(basename(point_data),"is not a csv file. Please provide a csv file."))
        mySpecies_occurrences = try(read.csv(SPECIES_OCC, header=TRUE), silent=TRUE)
        if(!inherits(mySpecies_occurrences,'try-error'))
          list.species <- gsub("\\.csv","",basename(SPECIES_OCC))
        else
          stop('...<< Unable to read species data >>...')
       },

       # directory
      {
        dn <- dir(OUTPUT_DIR, pattern="^Cleaned_data_by_", full.names=TRUE)
        if(length(dn)>0L)
          unlink(dn, recursive=TRUE)

        list.species <- list.files(SPECIES_OCC, pattern="\\.csv$", full.names=TRUE, recursive=TRUE) # select species
        if(length(list.species)==0L)
          stop(paste0("Unable to find csv files in directory:",SPECIES_OCC,". Please provide a directory with csv files."))
        empty.string <- nchar(list.species)==0L
        if(any(empty.string)){
          if(verbose)
            warning(paste("Removing",sum(empty.string),"file(s) with no species names."))
          list.species <- list.species[!empty.string]
        }
        mySpecies_occurrences <- SPECIES_OCC
      },

      # data
      {
        mySpecies_occurrences = try(read.csv(SPECIES_OCC, header=TRUE),silent=TRUE)
        if(!inherits(mySpecies_occurrences,'try-error')){
          file.removed <- file.remove(SPECIES_OCC)
          if(!file.removed){
            cat('\n\n')
            warning(paste('Unable to remove species name file:', SPECIES_OCC))
          }
          list.species <- gsub("\\.csv","",basename(SPECIES_OCC))
        }else{
          stop('...<< Unable to read species data >>...')
        }
      }
)

if(length(mySpecies_occurrences)==0L){
  stop("Oops !! Something probably went wrong when reading the species names...")
}else{
  cat('> Parameters check: OK...\n')
}

if(file.exists(SPECIES_ID)){
  mySpecies_ids = try(read.csv(SPECIES_ID, header=TRUE)[,1],silent=TRUE)
  mySpecies_ids_names <- try(read.csv(SPECIES_ID, header=TRUE)[,2],silent=TRUE)
  if(!inherits(mySpecies_ids,'try-error') & !inherits(mySpecies_ids_names,'try-error')){
    names(mySpecies_ids) <- mySpecies_ids_names
    file.removed <- file.remove(SPECIES_ID)
    if(!file.removed){
      cat('\n\n')
      warning(paste('Unable to remove species id file:', SPECIES_ID))
    }
  }
  else
    mySpecies_ids = NULL
}else{
  mySpecies_ids = NULL
}

cat('\n\n')

cat('#====================================================================\n')
cat('#== 2. Run TDWG cleaning tool\n')
cat('#====================================================================\n')
if(USE_TDWG){
  status ='both'
  cat(paste0('...cleaning occurrences by ',ifelse(status=='both','known',status),' range...'))
  cat('\n\n')
  clean.out <- TDWG::rangeCleaner(mySpecies_occurrences,
                     what='occurrences',
                     working_dir=OUTPUT_DIR,
                     species_id = NULL,
                     status='both',
                     initial_level=2,
                     do.parallel=MC_CORES>1 & length(list.species)>1,
                     use_name_matching=FALSE,
                     ncores = MC_CORES,
                     full_data = TRUE,
                     verbose=FALSE)
  if(is.null(clean.out))
    stop(paste('No available data after cleaning by',ifelse(status=='both','known',status),'range'))

  if(all(lengths(clean.out)==0))
    stop(paste('No available data after cleaning by',ifelse(status=='both','known',status),'range'))

  clean.out = clean.out[!lengths(clean.out)==0]

  if(mySpecies_occurrences==OUTPUT_DIR)
    clean_dir <- file.path(OUTPUT_DIR,paste0("Cleaned_data_by_",ifelse(status=='both','known',status),"_range"))
  else
    clean_dir <- OUTPUT_DIR

  list.species <- list.files(clean_dir, pattern="\\.csv$", full.names=TRUE, recursive=TRUE)
  is.empty <- lengths(sapply(list.species,readLines,USE.NAMES=FALSE))==0

  if(all(is.empty))
    stop('Oops ! something when wrong with the rangeCleaner function: all generated files are empty !! ')

  list.species <- list.species[!is.empty]
  OUTPUT_DIR <- clean_dir

}else{
  cat('...No cleaning by TDWG regions...')
}
cat('\n\n')

cat('#====================================================================\n')
cat('#== 3. Run standard cleaning tool\n')
cat('#====================================================================\n')
rm(mySpecies_occurrences)

cat('...Cleaning occurrence records...\n\n')

run.started.at <- Sys.time()
mclapply(list.species, function(j){

    verbose=FALSE

    out = try({

            # Get data for just that species
            if(file.exists(j)){
              species.data <- tryCatch(data.table::fread(j, showProgress=FALSE),error=function(err) return(NULL))
              k = gsub("\\.csv|_","",basename(j))
            }else return(NULL)

            # if an error occurred during the reading returns an NULL
            if(is.null(species.data) || ncol(species.data) < 2){
              if(verbose) warning(paste0("Cannot read file from species",k))
                dn <- file.path(OUTPUT_DIR,"problems")
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
                dn <- file.path(OUTPUT_DIR,"problems")
                if(!dir.exists(dn))
                  dir.create(dn, recursive =TRUE)
                fn <-file.path(dn,'Has_no_points.csv')
                if(!file.exists(fn))
                  write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
                else if(!k %in% read.csv(fn)$Species)
                  write.table(data.frame(Species=k), file=fn, row.names=FALSE, col.names = !file.exists(fn), sep=",",  append=TRUE)
              return(NULL)
            }
            # transform back to data.frame
            data.table::setDF(species.data)

            #---------------------------------------------------------
            ##= 1. Remove records beyond a certain age
            #---------------------------------------------------------
            clean <- species.data %>%

              dplyr::filter(is.na(year) | year >= CUTOFF_YEAR)

            if(nrow(clean)==0) return(NULL)

            #---------------------------------------------------------
            ##= 2. Remove cultivated specimens
            #---------------------------------------------------------
            # clean <- clean %>%
            #
            #   dplyr::filter(is.na(is_cultivated_observation) | is_cultivated_observation == "Yes"
            #                 & is.na(basisOfRecord) | basisOfRecord!= "Cultivated habitat")
            #
            # if(nrow(clean)==0) return(NULL)

            #---------------------------------------------------------
            ##= 3. Remove un-misgeoreferenced specimens
            #---------------------------------------------------------
            clean <- clean %>%

              dplyr::filter(!is.na(decimalLongitude)) %>%

              dplyr::filter(!is.na(decimalLatitude)) %>%

              # Occurences with impossible lat-long coordinates
              dplyr::filter(decimalLatitude < 90 | decimalLatitude > -90 | decimalLongitude < 180 | decimalLongitude > -180)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 4. Remove coordinates with large uncertainty in m
            #---------------------------------------------------------
            clean <- clean %>%

              dplyr::filter(coordinateUncertaintyInMeters <= CUTOFF_COORD_UNCERTAINTY | is.na(coordinateUncertaintyInMeters))

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 5. Remove coordinates with no decimal place (rounded)
            #---------------------------------------------------------
            clean <- clean %>%

              dplyr::filter((decimalLongitude*10) %% 10 != 0) %>%

              dplyr::filter((decimalLatitude*10) %% 10 != 0)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 6. Clean records whose individual counts are 0s
            #---------------------------------------------------------
            clean <- clean %>%

              dplyr::filter(individualCount > 0 | is.na(individualCount)) %>%

              dplyr::filter(individualCount < 99 | is.na(individualCount))

            if(nrow(clean)==0) stop('No more available information in the dataset')

            #=========================================================
            ##
            ##=  Start CoordinateCleaner cleaning process
            ##
            #=========================================================
            # Rename lon/lat fields
            clean <- clean %>%

              dplyr::rename(decimallongitude = 'decimalLongitude', decimallatitude = 'decimalLatitude')

            #---------------------------------------------------------
            ##= 7. Clean records with 0 long and lat
            #---------------------------------------------------------
            clean <- clean %>%

              CoordinateCleaner::cc_zero(buffer = 0.5, verbose=verbose)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 8. Clean records with equal latitude and longitude
            #      coordinates, either exact or absolute
            #---------------------------------------------------------
            clean <- clean %>%

              CoordinateCleaner::cc_equ(test = "absolute", verbose=verbose)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 9. Clean records outside reported country
            #---------------------------------------------------------
            # clean <- clean %>%
            #
            #   CoordinateCleaner::cc_coun(iso3='countryCode', ref=NULL, verbose=verbose)
            #
            # if(nrow(clean)==0) return(NULL)

            #---------------------------------------------------------
            ##= 10. Clean records with country and province centroids
            #---------------------------------------------------------
            clean <- clean %>%

              CoordinateCleaner::cc_cen(buffer = BUFFER_CENT, geod=TRUE, test='both', species='species', ref=NULL, verify=FALSE, verbose=verbose)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 11. Clean records in country capitals
            #---------------------------------------------------------
            clean <- clean %>%

              CoordinateCleaner::cc_cap(buffer = BUFFER_CAP, geod=TRUE, species='species', ref=NULL, verify=FALSE, verbose=verbose)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 12. Clean records with institutional coordinates
            #---------------------------------------------------------
            clean <- clean %>%

              CoordinateCleaner::cc_inst(buffer = BUFFER_INST, geod=TRUE, species='species', verify=FALSE, verify_mltpl = 10, verbose=verbose)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #---------------------------------------------------------
            ##= 13. Clean records with GBIF HQ coordinates
            #---------------------------------------------------------
            clean <- clean %>%

              CoordinateCleaner::cc_gbif(buffer = BUFFER_GBIF_HQ, geod=TRUE, species='species', verify=FALSE, verbose=verbose)

            if(nrow(clean)==0) stop('No more available information in the dataset') #return(NULL)

            #=========================================================
            ##
            ##=  End CoordinateCleaner cleaning process
            ##
            #=========================================================
            # Rename lon/lat fields
            clean <- clean %>%

              dplyr::rename(decimalLongitude = 'decimallongitude', decimalLatitude = 'decimallatitude')

            #---------------------------------------------------------
            ##= Save cleaned records to output folder
            #---------------------------------------------------------
            NEW_OUTPUT_DIR <- file.path(OUTPUT_DIR, "Cleaned_data_by_standardized_criteria")
            if(!dir.exists(NEW_OUTPUT_DIR))
              dir.create(NEW_OUTPUT_DIR, recursive=TRUE)

            write.csv(clean, file=file.path(NEW_OUTPUT_DIR, paste0(k,'.csv')), row.names=FALSE)

            clean

    }, silent=TRUE)

    if(inherits(out,'try-error')) cat(out)
    gc()

},mc.cores=MC_CORES, mc.preschedule = FALSE)
finished.at <- Sys.time()
runtime.elapsed <- finished.at - run.started.at
total.time.elapsed <- finished.at - started.at

cat('...Cleaning completed...\n\n')
cat(paste0('> End of Run on: ',Sys.time()),'\n\n');finished.at=proc.time()
cat(paste0('> Running time: ',as.numeric(runtime.elapsed),' ', attr(runtime.elapsed,'units'),'\n'))
cat(paste0('> Total time: ',as.numeric(total.time.elapsed),' ', attr(total.time.elapsed,'units'),'\n'))
cat(paste0('#',paste(rep('-',times=100),collapse="")))
cat('\n\n')
