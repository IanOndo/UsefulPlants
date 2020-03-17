## If the OS is Windows, set mclapply to the hackish version. Otherwise, leave the definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {mclapply.hack},
                    Linux   = {mclapply},
                    Darwin  = {mclapply})

#====================================================================
#== 0. Get parameters from command-line
#====================================================================

args = commandArgs(trailingOnly=TRUE)

for(arg in args) {

  arg_name = stringr::str_extract(arg,".+?(?=\\=)")

  switch(arg_name,

         'species_name' = {SPECIES_NAME = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'data_source'  = {DATA_SOURCE = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'download_dir' = {DOWNLOAD_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'output_dir'   = {OUTPUT_DIR = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'run_name'     = {RUN_NAME = gsub(".+?(?<=\\=)","",arg,perl=TRUE)},

         'mc_cores' 		= {MC_CORES = as.integer(gsub(".+?(?<=\\=)","",arg,perl=TRUE))}
  )
}
cat(paste0('#',paste(rep('-',times=100),collapse="")))
cat('\n')
cat(paste0("> Run ","'",RUN_NAME,"'",' started on: ',Sys.time()));
cat('\n')
cat('#====================================================================\n')
cat('#== 1. Check input parameters\n')
cat('#====================================================================\n\n')

# get the name of the run
run_name = RUN_NAME

# get the list of species from file
mySpecies_names = try(as.character(read.csv(SPECIES_NAME, header=TRUE)[,1]),silent=TRUE)

cat('> Species name(s):')
if(!inherits(mySpecies_names,'try-error')){
  if(grepl("queryOcc",SPECIES_NAME)){
      file.removed <- file.remove(SPECIES_NAME)
      if(!file.removed){
        cat('\n\n')
        warning(paste('Unable to remove species name file:', SPECIES_NAME))
      }
      cat('...OK\n')
  }
}else{
  stop('...<< Unable to read species names data >>...')
}
if(length(mySpecies_names)==0L)
  stop("Oops !! Something probably went wrong when reading the species names...")

# get the list of data sources to use
myData_sources  = unlist(stringr::str_split(DATA_SOURCE, pattern=","))

cat('> Data source(s):')
if(length(myData_sources)!=0L && sum(nchar(myData_sources))>0L){
  cat('...OK\n')
}else{
  stop('...<< Please provide at least one valid data source >>...')
}

# ensure that the directories provided exists
if(exists('DOWNLOAD_DIR'))
  if(!dir.exists(DOWNLOAD_DIR)) dir.create(DOWNLOAD_DIR, recursive=TRUE)
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)

# Is there any offline databases required ?
myOfflineData_sources <- intersect(myData_sources, c('biotime','rainbio', 'genesys', 'genesys','spLink','cwr_gbif'))
if(length(myOfflineData_sources)>0L){
  # Should some data sources be preformatted ?
  need_formatting = sapply(myOfflineData_sources, function(src) !file.exists(system.file("extdata/Occ_dir/data_formatted",paste0(src,"DATA.rds"), package='UsefulPlants')))
}else{
  need_formatting = FALSE
}
doFormatting = any(need_formatting)

cat('\n\n')

cat('#====================================================================\n')
cat('#== 2. Prepare databases\n')
cat('#====================================================================\n')
if(doFormatting){
  cat('...Formatting required databases...')
  cat('\n\n')
  commandArgs <- function(...) paste0(intersect(names(need_formatting)[need_formatting], c('biotime','rainbio', 'genesys','spLink','cwr_gbif')))
  assign('commandArgs',commandArgs,envir=.GlobalEnv)
  source(system.file("extdata/UsefulPlants_workflow/dbFormat.R",package='UsefulPlants'))
}else{
  cat('...No need to format required databases now...')
}

cat('\n\n')


cat('#====================================================================\n')
cat('#== 3. Pre-schedule database query\n')
cat('#====================================================================\n\n')
if(length(myOfflineData_sources)>0L){

  cat('#--------------------------------------------------------------------\n')
  cat('### 3.a search in which offline database to query your species\n')
  cat('#--------------------------------------------------------------------\n')
  list_of_formatted_databases <- list.files(system.file("extdata/Occ_dir/data_formatted", package='UsefulPlants'), pattern='\\.rds$', full.names=TRUE)
  if(length(list_of_formatted_databases)==0L)
    stop('<< Unable to find formatted database >>')

  myOfflineDataBases <- lapply(list_of_formatted_databases, readRDS)
  names(myOfflineDataBases) <- gsub('DATA\\.rds','', basename(list_of_formatted_databases))

  mydatabaseQueries <- sapply(myOfflineDataBases, is.keyval.exists, keyval=mySpecies_names)

  # compute number of species available by databases and total number of species available
  NbSpecies_by_offline_data_sources <- colSums(mydatabaseQueries)
  NbSpecies_available <- sum(apply(mydatabaseQueries,1,any))

  if(NbSpecies_available>0){
    cat(paste(NbSpecies_available,'species available in the local databases specified:\n'))

    for(k in 1:length(myOfflineDataBases))
      cat(paste('*',NbSpecies_by_offline_data_sources[k],paste0("species from ", toupper(names(myOfflineDataBases)[k])," database\n")))
  }else{
    cat('...<< None of your species are available in offline databases required >>...')
  }

}else{
  cat('...No need to pre-schedule required databases...')
}

cat('\n\n')

cat('#====================================================================\n')
cat('#== 4. Run queries\n')
cat('#====================================================================\n')
has_gbif = 'gbif' %in% myData_sources

if(has_gbif){
  #Creating a GBIF login
  myGBIFLogin <- list(user = "spironon",
                      email = "s.pironon@kew.org",
                      pwd = "Useful54");
  # if(is.null(myGBIFLogin)){
  #   if(length(myData_sources)==1)
  #     stop("Impossible to login to GBIF")
  #   else{
  #     warning("Impossible to login to GBIF. Gbif database will not be queried.")
  #     myData_sources<-myData_sources[-c('gbif')]
  #     has_gbif = 'gbif' %in% myData_sources
  #   }
  # }
  # set gbif user account credentials
  if(all(Sys.getenv(c('GBIF_USER','GBIF_PWD','GBIF_EMAIL'))==""))
    Sys.setenv(GBIF_USER=myGBIFLogin$'user', GBIF_PWD=myGBIFLogin$'pwd', GBIF_EMAIL=myGBIFLogin$'email')
}

cat('...Extracting occurrence records from', paste(myData_sources[myData_sources!='gbif'],collapse=", "),'...\n\n')

started.at <- Sys.time()
mclapply(c(1:length(mySpecies_names)), function(j){

  # get database to query
  offline_db  <- NULL
  if(exists('mydatabaseQueries'))
    offline_db  <- colnames(mydatabaseQueries)[mydatabaseQueries[j,]]
  online_db   <- setdiff(myData_sources, myOfflineData_sources)
  db_to_query <- c(offline_db, online_db[!online_db=="gbif"])

  # loop over databases
  for(db in db_to_query){

    out = switch(db,

                 'bien' = {
                   bienData <- queryBIEN(species_name = mySpecies_names[j])

                   if(!is.null(bienData)){

                     if(!dir.exists(file.path(OUTPUT_DIR, mySpecies_names[j])))
                       dir.create(file.path(OUTPUT_DIR, mySpecies_names[j]),recursive=TRUE)

                     write.table(bienData, file=file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv')),
                                 row.names=FALSE,
                                 col.names = !file.exists(file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv'))),
                                 sep=",",
                                 append=TRUE)
                   }
                 },

                 'genesys' = {

                   if(!exists('genesysDATA'))
                    genesysDATA <- readRDS(system.file("extdata/Occ_dir/data_formatted","genesysDATA.rds", package='UsefulPlants'))

                   if(!dir.exists(file.path(OUTPUT_DIR, mySpecies_names[j])))
                     dir.create(file.path(OUTPUT_DIR, mySpecies_names[j]),recursive=TRUE)

                    write.table(genesysDATA[mySpecies_names[j]], file=file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv')),
                                row.names=FALSE,
                                col.names = !file.exists(file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv'))),
                                sep=",",
                                append=TRUE)
                 },

                 'cwr_gbif' = {

                   if(!exists('cwr_gbifDATA'))
                     cwr_gbifDATA <- readRDS(system.file("extdata/Occ_dir/data_formatted","cwr_gbifDATA.rds", package='UsefulPlants'))

                   if(!dir.exists(file.path(OUTPUT_DIR, mySpecies_names[j])))
                     dir.create(file.path(OUTPUT_DIR, mySpecies_names[j]),recursive=TRUE)

                   write.table(cwr_gbifDATA[mySpecies_names[j]], file=file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv')),
                               row.names=FALSE,
                               col.names = !file.exists(file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv'))),
                               sep=",",
                               append=TRUE)
                 },

                 'spLink' = {

                   if(!exists('spLinkDATA'))
                     spLinkDATA <- readRDS(system.file("extdata/Occ_dir/data_formatted","spLinkDATA.rds", package='UsefulPlants'))

                   if(!dir.exists(file.path(OUTPUT_DIR, mySpecies_names[j])))
                     dir.create(file.path(OUTPUT_DIR, mySpecies_names[j]),recursive=TRUE)

                   write.table(spLinkDATA[mySpecies_names[j]], file=file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv')),
                               row.names=FALSE,
                               col.names = !file.exists(file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv'))),
                               sep=",",
                               append=TRUE)
                 },

                 'rainbio' = {

                   if(!exists('rainbioDATA'))
                     rainbioDATA <- readRDS(system.file("extdata/Occ_dir/data_formatted","rainbioDATA.rds", package='UsefulPlants'))

                   if(!dir.exists(file.path(OUTPUT_DIR, mySpecies_names[j])))
                     dir.create(file.path(OUTPUT_DIR, mySpecies_names[j]),recursive=TRUE)

                   write.table(rainbioDATA[mySpecies_names[j]], file=file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv')),
                               row.names=FALSE,
                               col.names = !file.exists(file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv'))),
                               sep=",",
                               append=TRUE)
                 },

                 'biotime' = {

                   if(!exists('biotimeDATA'))
                     biotimeDATA <- readRDS(system.file("extdata/Occ_dir/data_formatted","biotimeDATA.rds", package='UsefulPlant'))

                   if(!dir.exists(file.path(OUTPUT_DIR, mySpecies_names[j])))
                     dir.create(file.path(OUTPUT_DIR, mySpecies_names[j]),recursive=TRUE)

                   write.table(biotimeDATA[mySpecies_names[j]], file=file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv')),
                               row.names=FALSE,
                               col.names = !file.exists(file.path(OUTPUT_DIR, mySpecies_names[j], paste0(mySpecies_names[j],'.csv'))),
                               sep=",",
                               append=TRUE)
                 }

    )
    if(inherits(out,'try-error')) cat(out)
    gc()
  }
},mc.cores=MC_CORES, mc.preschedule = FALSE)

if(has_gbif){

  cat('\n\n...Extracting occurrence records from gbif...\n\n')

  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  nchunk = 3
  chunked.list <- chunk(mySpecies_names, nchunk)

  mclapply(c(1:length(chunked.list)), function(j){

    out = try({

      # set gbif user account credentials
      if(all(Sys.getenv(c('GBIF_USER','GBIF_PWD','GBIF_EMAIL'))==""))
        Sys.setenv(GBIF_USER=myGBIFLogin$'user', GBIF_PWD=myGBIFLogin$'pwd', GBIF_EMAIL=myGBIFLogin$'email')

      if(!dir.exists(DOWNLOAD_DIR))
        dir.create(DOWNLOAD_DIR,recursive=TRUE)

      gbifData <- queryGBIF(species_name = chunked.list[[j]], gbif_download_dir = DOWNLOAD_DIR, time_out=10800, check_output=DOWNLOAD_DIR, index=j)

      if(!is.null(gbifData)){

        for(sp in 1:length(chunked.list[[j]])){

          if(!dir.exists(file.path(OUTPUT_DIR, chunked.list[[j]][sp])))
            dir.create(file.path(OUTPUT_DIR, chunked.list[[j]][sp]),recursive=TRUE)

          write.table(gbifData[chunked.list[[j]][sp]], file=file.path(OUTPUT_DIR, chunked.list[[j]][sp], paste0(chunked.list[[j]][sp],'.csv')),
                      row.names=FALSE,
                      col.names = !file.exists(file.path(OUTPUT_DIR, chunked.list[[j]][sp], paste0(chunked.list[[j]][sp],'.csv'))),
                      sep=",",
                      append=TRUE)
        }

      }

   }, silent=TRUE)

   if(inherits(out,'try-error')) cat(out)
   gc()

  },mc.cores=MC_CORES, mc.preschedule = FALSE)

}

finished.at <- Sys.time()
time.elapsed <- finished.at - started.at

cat('...Extraction completed...\n\n')
cat(paste0('> End of Run ',"'",run_name,"'",' on: ',Sys.time()),'\n\n');finished.at=proc.time()
cat(paste0('> Running time: ',as.numeric(time.elapsed),' ', attr(time.elapsed,'units'),'\n'))
cat(paste0('#',paste(rep('-',times=100),collapse="")))
cat('\n\n')
