data_sources = commandArgs(trailingOnly=TRUE)

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

cat('> reading and formatting offline databases...')

# start of loop
for(src in data_sources) {

  switch(src,

       'genesys' = {

         #---------------
         #= reading
         #---------------
         genesysDATA <- try(data.table::fread(system.file("extdata/Occ_dir/data_sources/GENESYSdata.csv", package='UsefulPlants'), header=TRUE, showProgress=FALSE, na.strings=c("",NA),
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

         if(!inherits(genesysDATA,'error')){

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
         cwr_gbifDATA <- try(data.table::fread(system.file("extdata/Occ_dir/data_sources/CWRdatabaseGBIF.csv", package='UsefulPlants'), header=TRUE, showProgress=FALSE, quote="", na.strings=c("",NA)),
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
                                      "individualCount"), silent=TRUE)

         if(!inherits(cwr_gbifDATA,'error')){

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
           cwr_gbifDATA[,`:=`(fullname = ifelse(is.na(infraspecificEpithet), paste(species), paste(species, taxonRank, infraspecificEpithet)),
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
         spLinkDATA <- try(data.table::fread(system.file("extdata/Occ_dir/data_sources/SpeciesLink.txt", package='UsefulPlants'),header=TRUE, showProgress=FALSE, na.strings=c("",NA)), silent=TRUE)

         if(!inherits(spLinkDATA,'error')){
          # formatting
          spLinkDATA <- na.omit(spLinkDATA, cols= c("latitude", "longitude"))

          spLink_flag = TRUE
         }

       },

       'rainbio' = {

         # reading
         rainbioDATA <- try(data.table::fread(system.file("extdata/Occ_dir/data_sources/RAINBIO.csv", package='UsefulPlants'), header=TRUE, showProgress=FALSE, drop=1,
                                            select=c("tax_sp_level",
                                              "species",
                                              "decimalLatitude",
                                              "decimalLongitude",
                                              "iso3lonlat",
                                              "basisOfRecord",
                                              "institutionCode",
                                              "catalogNumber",
                                              "coly")), silent=TRUE)

         if(!inherits(rainbioDATA,'error')){

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
          biotimeDATA  <- try(data.table::fread(system.file("extdata/Occ_dir/data_sources/BioTIMEQuery02_04_2018.csv", package='UsefulPlants'), header=TRUE, showProgress=FALSE,
                                                select=c("GENUS_SPECIES",
                                                         "LATITUDE",
                                                         "LONGITUDE",
                                                         "YEAR",
                                                         "sum.allrawdata.ABUNDANCE",
                                                         "SAMPLE_DESC")), silent=TRUE)

          if(!inherits(biotimeDATA,'error')){

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

if(!all(c(genesys_flag, spLink_flag, rainbio_flag, biotime_flag))){
  cat('\n\n')
  stop('<< Impossible to read data from databases required >>')
}

cat('\n\n')

#=====================================================================
#=2.Saving offline databases (downloaded from online)
#=====================================================================
cat('> saving data...')

db_dir <- system.file("extdata/Occ_dir/data_formatted", package='UsefulPlants')
if(!dir.exists(db_dir))
  dir.create(db_dir, recursive=TRUE)

if(genesys_flag){
  cat('\n')
  cat('==> genesys: ')
  if(inherits(try(saveRDS(genesysDATA, system.file("extdata/Occ_dir/data_formatted/genesysDATA.rds", package='UsefulPlants')), silent=TRUE),'error')){
    warning('Unable to save genesys database')
    genesys_flag = FALSE
  }else
    cat('done')
}

if(cwr_db_flag){
  cat('\n')
  cat('==> cwr: ')
  if(inherits(try(saveRDS(cwr_gbifDATA, system.file("extdata/Occ_dir/data_formatted/cwr_gbifDATA.rds", package='UsefulPlants')), silent=TRUE),'error')){
    warning('Unable to save cwr database')
    cwr_db_flag = FALSE
  }else
    cat('done')
}

if(spLink_flag){
  cat('\n')
  cat('==> spLink: ')
  if(inherits(try(saveRDS(spLinkDATA, system.file("extdata/Occ_dir/data_formatted/spLinkDATA.rds", package='UsefulPlants')), silent=TRUE),'error')){
    warning('Unable to save spLink database')
    spLink_flag = FALSE
  }else
    cat('done')
}

if(rainbio_flag){
  cat('\n')
  cat('==> rainbio: ')
  if(inherits(try(saveRDS(rainbioDATA, system.file("extdata/Occ_dir/data_formatted/rainbioDATA.rds", package='UsefulPlants')), silent=TRUE),'error')){
    warning('Unable to save rainbio database')
    rainbio_flag = FALSE
  }else
    cat('done')
}

if(biotime_flag){
  cat('\n')
  cat('==> biotime: ')
  if(inherits(try(saveRDS(biotimeDATA, system.file("extdata/Occ_dir/data_formatted/biotimeDATA.rds", package='UsefulPlants')), silent=TRUE),'error')){
    warning('Unable to save biotime database')
    biotime_flag = FALSE
  }else
    cat('done')
}

if(!all(c(genesys_flag, spLink_flag, rainbio_flag, biotime_flag, cwr_db_flag))){
  cat('\n\n')
  stop('<< Impossible to save formatted databases >>')
}

cat('DONE')

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



