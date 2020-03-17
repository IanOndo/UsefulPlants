db_sources = commandArgs(trailingOnly=TRUE)

#=====================================================================
#=1. Reading and formatting offline databases (downloaded from online)
#=====================================================================

# define column names
colNames = c("species",
             "fullname",
             "decimalLongitude",
             "decimalLatitude",
             "countryCode",
             "coordinateUncertaintyInMeters",
             "year",
             "individualCount",
             "gbifID",
             "basisOfRecord",
             "institutionCode",
             "establishmentMeans",
             "is_cultivated_observation",
             "sourceID")

# flags
genesys_flag = spLink_flag = rainbio_flag = biotime_flag = cwr_db_flag = FALSE
overwrite = FALSE

cat('> reading and formatting offline databases...\n')

# start of loop
for(src in db_sources) {

  switch(src,

       'genesys' = {

         #---------------
         #= reading
         #---------------
         genesysDATA <- try(data.table::fread(paste("unzip -p",system.file("extdata/Occ_dir/data_sources/GENESYSdata.zip", package='UsefulPlants')), header=TRUE, showProgress=FALSE, na.strings=c("",NA),
                                              select= c("GENUS",
                                              "SPECIES",
                                              "SUBTAXA",
                                              "DECLONGITUDE",
                                              "DECLATITUDE",
                                              "ORIGCTY",
                                              "COLLSRC",
                                              "INSTCODE",
                                              "UUID",
                                              "COLLDATE",
                                              "COORDUNCERT",
                                              "SAMPSTAT")), silent=TRUE)

         if(!inherits(genesysDATA,'try-error')){

           #---------------
           #= formatting
           #---------------
           # change the column names
           data.table::setnames(genesysDATA, c("genus", "species", "subtaxa", "decimalLongitude", "decimalLatitude", "countryCode", "basisOfRecord", "institutionCode", "gbifID", "year", "coordinateUncertaintyInMeters", "is_cultivated_observation"))
           # format the data: concatenate genus+species, create full species name, extract the year of the collection, transform information about cultivated observation and basis of records in categorical variable
           genesysDATA[, `:=`(species = paste(genus, species),
                              fullname = ifelse(is.na(subtaxa), paste(genus, species), paste(genus, species, subtaxa)),
                              year = as.numeric(stringr::str_extract(year,"[[:digit:]]{4}")),
                              is_cultivated_observation = ifelse(is_cultivated_observation==999, NA, ifelse(is_cultivated_observation > 300, "Yes", "No")),
                              basisOfRecord = cut(basisOfRecord, breaks=c(10, 16, 29, 59, 63),  labels=c("Wild habitat", NA, "Cultivated habitat","Wild habitat")))]
           # create additional columns: 'establishmentMeans', 'individualCount' and 'sourceID' and delete 'genus' and 'subtaxa'
           genesysDATA[ ,`:=`(genus = NULL,
                              subtaxa = NULL,
                              establishmentMeans = NA,
                              individualCount = NA,
                              sourceID = 'GENESYS')]
           # set column order
           data.table::setcolorder(genesysDATA, colNames)
           # set the key to the species column to enable fast binary search
           data.table::setkey(genesysDATA, 'species')
           genesys_flag = TRUE
         }else{
           cat('\n\n')
           cat('...<< Unable to read data from genesys database >>...')
         }


       },

       'cwr_gbif' = {

         #---------------
         #= reading
         #---------------
         cwr_gbifDATA <- try(data.table::fread(paste("unzip -p",system.file("extdata/Occ_dir/data_sources/CWRdatabaseGBIF.zip", package='UsefulPlants')), header=TRUE, showProgress=FALSE, quote="", na.strings=c("",NA),
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
                                      "establishmentMeans")), silent=TRUE)

         if(!inherits(cwr_gbifDATA,'try-error')){

           #---------------
           #= formatting
           #---------------
           # remove NA coordinates
           cwr_gbifDATA  <- na.omit(cwr_gbifDATA, cols= c("decimalLatitude", "decimalLongitude"))
           # format the data: concatenate genus+species, create full species name, extract the year of the collection, transform information about cultivated observation and basis of records in categorical variable
           cwr_gbifDATA[, `:=`(taxonRank = ifelse(taxonRank %in% c("SPECIES","GENUS"), NA, ifelse(taxonRank=="FORM", "f.", ifelse(taxonRank=="SUBSPECIES", "subsp.", ifelse(taxonRank=="VARIETY","var.", taxonRank)))),
                              countryCode = countrycode::countrycode(countryCode, origin =  'iso2c', destination = 'iso3c', nomatch = NA),
                              establishmentMeans = ifelse(establishmentMeans=="INTRODUCED","Introduced",ifelse(establishmentMeans== "NATIVE", "Native", establishmentMeans)))]
           # create 'fullname', 'is_cultivated_observation', 'individualCount' and 'sourceID' columns
           cwr_gbifDATA[, `:=`(fullname = ifelse(is.na(infraspecificEpithet), paste(species), paste(species, taxonRank, infraspecificEpithet)),
                          is_cultivated_observation = "No",
                          individualCount = NA,
                          sourceID = 'CWRGBIF')]
           # delete taxonRank and infraspecificEpithet columns
           cwr_gbifDATA[, `:=`(taxonRank = NULL,
                           infraspecificEpithet = NULL)]
           # set column order
           data.table::setcolorder(cwr_gbifDATA, colNames)
           # set the key to the species column to enable fast binary search
           data.table::setkey(cwr_gbifDATA, 'species')
           cwr_db_flag = TRUE
         }else{
           cat('\n\n')
           cat('...<< Unable to read data from crop wild relative database >>...')
         }


       },

       'spLink' = {

         # reading
         spLinkDATA <- try(data.table::fread(paste("unzip -p",system.file("extdata/Occ_dir/data_sources/SpeciesLinkData.zip", package='UsefulPlants')),header=TRUE, showProgress=FALSE, encoding = 'UTF-8', na.strings=c("",NA),
                                             select= c("query",
                                                       "scientificname",
                                                       "longitude",
                                                       "latitude",
                                                       "country",
                                                       "coordinateprecision",
                                                       "yearcollected",
                                                       "individualcount",
                                                       "catalognumber",
                                                       "basisofrecord",
                                                       "institutioncode"
                                             )), silent=TRUE)

         if(!inherits(spLinkDATA,'try-error')){
          # formatting
          # change the column names
          data.table::setnames(spLinkDATA, c("species", "fullname", "decimalLatitude","decimalLongitude", "countryCode", "coordinateUncertaintyInMeters", "year", "individualCount", "gbifID", "basisOfRecord", "institutionCode"))

          spLinkDATA <- na.omit(spLinkDATA, cols= c("decimalLatitude", "decimalLongitude"))
          spLinkDATA[ ,`:=`(species = paste0(substr(species,1,1), tolower(substr(species,2,nchar(species)))),
                           countryCode = ifelse(countryCode=='Brasil','Brazil',countryCode),
                           basisOfRecord = ifelse(basisOfRecord=="S","SPECIMEN",ifelse(basisOfRecord== "O", "OBSERVATION", basisOfRecord)),
                           is_cultivated_observation = NA,
                           establishmentMeans = NA,
                           sourceID = 'spLink')]
          spLinkDATA[ ,countryCode:= countrycode::countrycode(countryCode, origin =  'country.name', destination = 'iso3c', nomatch = NA)]
          # set column order
          data.table::setcolorder(spLinkDATA, colNames)
          # set the key to the species column to enable fast binary search
          data.table::setkey(spLinkDATA, 'species')
          spLink_flag = TRUE
         }else{
           cat('\n\n')
           cat('...<< Unable to read data from spLink database >>...')
         }

       },

       'rainbio' = {

         # reading
         rainbioDATA <- try(data.table::fread(paste("unzip -p",system.file("extdata/Occ_dir/data_sources/rainbioyear.zip", package='UsefulPlants')), header=TRUE, showProgress=FALSE,
                                            select=c("tax_sp_level",
                                              "species",
                                              "decimalLatitude",
                                              "decimalLongitude",
                                              "iso3lonlat",
                                              "basisOfRecord",
                                              "institutionCode",
                                              "catalogNumber",
                                              "coly")), silent=TRUE)

         if(!inherits(rainbioDATA,'try-error')){

           #---------------
           #= formatting
           #---------------
           # change the column names
           data.table::setnames(rainbioDATA, c("species", "fullname", "decimalLatitude","decimalLongitude", "countryCode", "basisOfRecord", "institutionCode", "gbifID", "year"))
           # remove NA coordinates
           rainbioDATA <- na.omit(rainbioDATA, cols= c("decimalLatitude", "decimalLongitude"))
           # create additional columns: 'coordinateUncertaintyInMeters' , 'is_cultivated_observation', 'establishmentMeans', 'individualCount' and 'sourceID'
           rainbioDATA[ ,`:=`(coordinateUncertaintyInMeters = NA,
                              is_cultivated_observation = "No",
                              establishmentMeans = NA,
                              individualCount = NA,
                              sourceID = 'RAINBIO')]
           # set column order
           data.table::setcolorder(rainbioDATA, colNames)
           # set the key to the species column to enable fast binary search
           data.table::setkey(rainbioDATA, 'species')
           rainbio_flag = TRUE
         }else{
           cat('\n\n')
           cat('...<< Unable to read data from rainbio database >>...')
         }

       },

       'biotime' = {

          # reading
          biotimeDATA  <- try(data.table::fread(paste("unzip -p",system.file("extdata/Occ_dir/data_sources/BioTIMEQuery02_04_2018.zip", package='UsefulPlants')), header=TRUE, showProgress=FALSE,
                                                select=c("GENUS_SPECIES",
                                                         "LATITUDE",
                                                         "LONGITUDE",
                                                         "YEAR",
                                                         "sum.allrawdata.ABUNDANCE",
                                                         "SAMPLE_DESC")), silent=TRUE)

          if(!inherits(biotimeDATA,'try-error')){

            #---------------
            #= formatting
            #---------------
            # change the column names
            data.table::setnames(biotimeDATA, c("species","decimalLatitude","decimalLongitude", "year", "individualCount", "gbifID"))
            # remove NA coordinates
            biotimeDATA  <- na.omit(biotimeDATA, cols= c("decimalLatitude", "decimalLongitude"))
            # create additional columns: 'fullname', 'coordinateUncertaintyInMeters','institutionCode', 'countryCode', 'basisOfRecord', 'establishmentMeans', 'is_cultivated_observation' and 'sourceID'
            biotimeDATA[ ,`:=`(fullname = species,
                               coordinateUncertaintyInMeters = NA,
                               institutionCode = NA,
                               countryCode = NA,
                               basisOfRecord = NA,
                               establishmentMeans = NA,
                               is_cultivated_observation = "No",
                               sourceID = 'BIOTIME')]
            # set column order
            data.table::setcolorder(biotimeDATA, colNames)
            # set the key to the species column to enable fast binary search
            data.table::setkey(biotimeDATA, 'species')
            biotime_flag = TRUE
          }else{
            cat('\n\n')
            cat('...<< Unable to read data from biotime database >>...')
          }

       }

  )

}
# end of loop

if(!any(c(genesys_flag, spLink_flag, rainbio_flag, biotime_flag, cwr_db_flag))){
  cat('\n\n')
  stop('<< Impossible to read data from databases required >>')
}

cat('\n\n')

#=====================================================================
#=2.Saving offline databases (downloaded from online)
#=====================================================================
cat('> saving data...')

db_dir <- file.path(path.package('UsefulPlants'), "extdata/Occ_dir/data_formatted")
if(!dir.exists(db_dir))
  dir.create(db_dir, recursive=TRUE)

if(genesys_flag){
  cat('\n')
  cat('==> genesys: ')

  if(file.exists(file.path(db_dir,"genesysDATA.rds")) & !overwrite){
    genesysDATA0 <- readRDS(file.path(db_dir,"genesysDATA.rds"))
    genesysDATA <- rbind(genesysDATA0, genesysDATA)
  }

  if(inherits(try(saveRDS(genesysDATA, file.path(db_dir,"genesysDATA.rds")), silent=TRUE),'try-error')){
    cat('failed\n')
    warning('Unable to save genesys database')
    genesys_flag = FALSE
  }else
    cat('done\n')
}

if(cwr_db_flag){
  cat('\n')
  cat('==> cwr: ')
  if(file.exists(file.path(db_dir,"cwr_gbifDATA.rds")) & !overwrite){
    cwr_gbifDATA0 <- readRDS(file.path(db_dir,"cwr_gbifDATA.rds"))
    cwr_gbifDATA <- rbind(cwr_gbifDATA0, cwr_gbifDATA)
  }
  if(inherits(try(saveRDS(cwr_gbifDATA, file.path(db_dir,"cwr_gbifDATA.rds")), silent=TRUE),'try-error')){
    cat('failed\n')
    warning('Unable to save cwr database')
    cwr_db_flag = FALSE
  }else
    cat('done\n')
}

if(spLink_flag){
  cat('\n')
  cat('==> spLink: ')
  if(file.exists(file.path(db_dir,"spLinkDATA.rds")) & !overwrite){
    spLinkDATA0 <- readRDS(file.path(db_dir,"spLinkDATA.rds"))
    spLinkDATA <- rbind(spLinkDATA0, spLinkDATA)
  }
  if(inherits(try(saveRDS(spLinkDATA, file.path(db_dir,"spLinkDATA.rds")), silent=TRUE),'try-error')){
    cat('failed\n')
    warning('Unable to save spLink database')
    spLink_flag = FALSE
  }else
    cat('done\n')
}

if(rainbio_flag){
  cat('\n')
  cat('==> rainbio: ')
  if(file.exists(file.path(db_dir,"rainbioDATA.rds")) & !overwrite){
    rainbioDATA0 <- readRDS(file.path(db_dir,"rainbioDATA.rds"))
    rainbioDATA <- rbind(rainbioDATA0, rainbioDATA)
  }
  if(inherits(try(saveRDS(rainbioDATA, file.path(db_dir,"rainbioDATA.rds")), silent=TRUE),'error')){
    cat('failed\n')
    warning('Unable to save rainbio database')
    rainbio_flag = FALSE
  }else
    cat('done\n')
}

if(biotime_flag){
  cat('\n')
  cat('==> biotime: ')
  if(file.exists(file.path(db_dir,"biotimeDATA.rds")) & !overwrite){
    biotimeDATA0 <- readRDS(file.path(db_dir,"biotimeDATA.rds"))
    biotimeDATA <- rbind(biotimeDATA0, biotimeDATA)
  }
  if(inherits(try(saveRDS(biotimeDATA, file.path(db_dir,"biotimeDATA.rds")), silent=TRUE),'error')){
    cat('failed\n')
    warning('Unable to save biotime database')
    biotime_flag = FALSE
  }else
    cat('done\n')
}

if(!any(c(genesys_flag, spLink_flag, rainbio_flag, biotime_flag, cwr_db_flag))){
  cat('\n\n')
  stop('<< Impossible to save formatted databases >>')
}

cat('> ...data saved')

cat('\n\n')

#
# genesysDATA <- try(data.table::fread('C:/Users/io03kg/Desktop/UsefulPlantsProject/UsefulPlants/inst/extdata/Occ_dir/data_sources/GENESYSdata.csv', header=TRUE, showProgress=FALSE, na.strings=c("",NA),
#                    select= c("GENUS",
#                              "SPECIES",
#                              "SUBTAXA",
#                              "DECLONGITUDE",
#                              "DECLATITUDE",
#                              "ORIGCTY",
#                              "COLLSRC",
#                              "INSTCODE",
#                              "UUID",
#                              "COLLDATE",
#                              "COORDUNCERT",
#                              "SAMPSTAT")),silent=TRUE)
#
# spLinkDATA <- try(data.table::fread('C:/Users/io03kg/Desktop/UsefulPlantsProject/UsefulPlants/inst/extdata/Occ_dir/data_sources/SpeciesLinkData.txt', header=TRUE, showProgress=FALSE, na.strings=c("",NA)), silent=TRUE)
#
# rainbioDATA <- try(data.table::fread('C:/Users/io03kg/Desktop/UsefulPlantsProject/UsefulPlants/inst/extdata/Occ_dir/data_sources/RAINBIO.csv', header=TRUE, showProgress=FALSE,
#                                      select=c("tax_sp_level",
#                                               "species",
#                                               "decimalLatitude",
#                                               "decimalLongitude",
#                                               "iso3lonlat",
#                                               "basisOfRecord",
#                                               "institutionCode",
#                                               "catalogNumber",
#                                               "coly")), silent=TRUE)
#
# biotimeDATA <- try(data.table::fread('C:/Users/io03kg/Desktop/UsefulPlantsProject/UsefulPlants/inst/extdata/Occ_dir/data_sources/BioTIMEQuery02_04_2018.csv', header=TRUE, showProgress=FALSE,
#                                      select=c("GENUS_SPECIES",
#                                               "LATITUDE",
#                                               "LONGITUDE",
#                                               "YEAR",
#                                               "sum.allrawdata.ABUNDANCE",
#                                               "SAMPLE_DESC")), silent=TRUE)



