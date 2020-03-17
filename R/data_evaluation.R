#'                                 Cross-validation with geographically masked blocks of occurrence records for evaluation
#'
#' Partition both occurrence records and background into evaluation bins based on some spatial rules.
#'
#'
#' @param loc_dat A two-column matrix,a data.frame, a data.table or a csv file with longitude and latitude coordinates of occurrence records.
#' @param env_dat A RasterStack object or a character string specifying the path to a directory with environmental raster layers.
#' @param k A numeric integer specifying the number of cluster required.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param bg_masks A RasterStack object of the background to be partitionned.
#' @param species_name [Optional] A character string specifying the name of species. If \code{NULL}, will be set to `Unknown`.
#' @param maxent_settings A list of maxent flags.
#' @param do.mask A logical
#' @param ... Additional parameters to be passed to \code{maxent} or \code{make_geographic_block}
#' @export
block_cv_maxent <- function(loc_dat, env_dat,
                            k=3,
                            coordHeaders=NULL,
                            bg_masks=NULL,
                            outputdir,
                            species_name=NULL,
                            maxent_settings=list(),
                            eval_metrics=c("auc","omission_rate","tss"),
                            varying_parameter_name=NULL,
                            do.mask=TRUE, do.parallel=TRUE, ncores=NULL, ...){

  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  if(missing(env_dat))
    stop("Please provide a path to your environmental layers.")

  if(missing(outputdir))
    stop("Output directory is missing. Please select a directory for maxent outputs.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!loc_file_flag & !loc_data_flag)
    stop('Unable to read location data. Please provide valid location data')

  if(loc_file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(species_name)){
    if(loc_file_flag){
      species_name=gsub(".csv","",basename(loc_dat))
      if(nchar(species_name)==0L)
        species_name="Unknown"
    }
    else
      species_name="Unknown"
  }

  if(loc_file_flag)
    loc_dat = read.csv(loc_dat, header=TRUE)

  if(!dir.exists(outputdir))
    stop("The directory ",outputdir," does not exist or is not accessible.")

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }

  env_file_flag = tryCatch(file.exists(env_dat), error=function(err) FALSE) && !tryCatch(dir.exists(env_dat), error=function(err) FALSE)
  env_data_flag = !env_file_flag & any(inherits(loc_dat, c("RasterStack","RasterLayer")))

  if(!env_file_flag & !env_data_flag)
    stop('Unable to read environmental data. Please provide valid environmental data.')

  sf_loc_data <- loc_dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))

  maxent.args = bg_args = list()
  dot.args <- list(...)
  if(length(dot.args)>0L){
    arg.names <- names(dot.args)
    if(any(c("method","grid_res","coastline") %in% arg.names))
      bg_args <- dot.args[arg.names%in%c("method","grid_res","coastline")]
    if(any(c("path_to_maxent","path_to_java","coastline") %in% arg.names))
      maxent_args <- dot.args[arg.names%in%c("path_to_maxent","path_to_java","coastline")]
  }
  if(length(maxent_args)>0L)
    maxent_settings <- c(maxent_settings, list(maxent_args))

  varying_parameter = FALSE
  if(!is.null(varying_parameter_name)){
    if(!varying_parameter_name %in% names(maxent_settings)){
      warning(varying_parameter_name, "is not a valid maxent parameter name.")
    }
    if(length(maxent_settings[[varying_parameter_name]])<2){
      warning("The parameter ",varying_parameter_name, " contains one single value. ")
    }
    varying_parameter = TRUE
  }

  if(!is.null(bg_masks)){
    if(!inherits(bg_mask,'RasterStack'))
      stop('Argument bg_masks must be a RasterStack object.')
    if(nlayers(bg_masks)!=k)
      stop('Number of layers in bg_masks must match the number of cluster k.')
  }
  else {
    #==================================
    #= 0. Create a background area with
    #     geographically masked blocks
    #==================================
    if(length(bg_args)>0L)
      bg_masks <- try(do.call(make_geographic_block, list(k=k, bg_args)), silent=TRUE)
    else
      bg_masks <- try(make_geographic_block(kdat=loc_dat, k=k), silent=TRUE)

    if(inherits(bg_masks,"try-error"))
      stop("Unable to create a geographically masked background dataset.")
  }
  #================================================
  #= 1. Create directories for each block partition
  #================================================
  comb <- combn(1:k,k-1)
  k_fold_names <- apply(comb,2,function(x) paste(letters[x],collapse="_"))

  dirs_to_create 		<- file.path(outputdir,k_fold_names)
  new.dirs.created 	<- sapply(dirs_to_create, function(x) if(!dir.exists(x)) dir.create(x) else TRUE)

  if(!any(new.dirs.created))
    stop("Unable to write output directories.")

  #===========================================
  #= 2. Cross-validate the model across blocks
  #===========================================
  cv_out <- foreach::foreach(j=1:k, .packages=c("UsefulPlants","rmaxent","dplyr"), .combine=rbind) %dopar% {

    require(dplyr)
    #--------------------------
    #= a. get the training data
    #--------------------------

    # training occurrence records
    training_block <- bg_masks %>%
      subset(layer %in% comb[,j])
    is_training_loc <- suppressMessages(sf_loc_data %>%
                                          sf::st_intersects(training_block) %>%
                                          purrr::map_lgl(., function(x) length(x) > 0L))
    training_loc <- data.frame(species=rep(species_name, sum(is_training_loc))) %>%
      dplyr::bind_cols(as.data.frame(sf::st_coordinates(sf_loc_data[is_training_loc,])))

    # save samples
    trainingsampfile = file.path(dirs_to_create[j],'sptrain.csv')
    write.csv(training_loc, trainingsampfile, row.names=FALSE)

    # environmental layers
    if(env_data_flag){
      env_bg <- stars::st_as_stars(env_dat)
      names(env_bg) <- strip_extension(names(env_bg))
      training_bg <- env_bg %>%
        `[`(training_block)
    }
    else if(env_file_flag){
      env_bg <- try(env_dat %>%
                      read_layers(),silent=TRUE)
      names(env_bg) <- strip_extension(names(env_bg))
      training_bg <- env_bg %>%
                           `[`(training_block)
      if(inherits(training_bg,"try-error"))
        stop("Unable to read environmental layers from ",env_dat,". Please make sure all the layers have the same dimensions.")
    }
    dn <- file.path(dirs_to_create[j],'env')
    if(!dir.exists(dn)) dir.create(dn)

    # save layers
    for(l in 1:length(training_bg)){
      training_bg_cvrt <- as(training_bg[l],"Raster")
      raster::writeRaster(training_bg_cvrt, file.path(dirs_to_create[j],'env',paste0(names(training_bg)[l],'.asc')), format='ascii', overwrite=TRUE)
    }

    #-------------------------
    #= b. get the testing data
    #-------------------------
    testing_loc <- data.frame(species=rep(species_name, sum(!is_training_loc))) %>%
      dplyr::bind_cols(as.data.frame(sf::st_coordinates(sf_loc_data[!is_training_loc,])))

    # add environmental variables
    testing_loc <- testing_loc %>%
      dplyr::bind_cols(UsefulPlants::map_lfd(env_bg, function(x) raster::extract(as(x,"Raster"),.[,-1])))

    # save samples with data
    testingsampfile = file.path(dirs_to_create[j],'sptest.csv')
    write.csv(testing_loc, testingsampfile, row.names=FALSE)

    #-----------------------------
    #= c. run the cross-validation
    #-----------------------------
    outputdir_path = dirs_to_create[j]

    if(varying_parameter){
      # create new directories
      param_dn = paste(varying_parameter_name, maxent_settings[[varying_parameter_name]],sep="_")
      outputdir_path = file.path(outputdir_path, param_dn)
      outdirs_created = sapply(outputdir_path, function(x) if(!dir.exists(x)) dir.create(x))
      if(!all(outdirs_created)) stop("Unable to create directories", outputdir_path[!outdirs_created])
    }

    cv.model.args <-Reduce(append,list(
                                  list(sp_points=trainingsampfile,

                                  env_layers=file.path(dirs_to_create[j],'env'),

                                  outputdir=outputdir_path[1]),

                                  list(testsamplesfile=testingsampfile),

                                  replace(maxent_settings,names(maxent_settings)==varying_parameter_name,maxent_settings[[varying_parameter_name]][1]),

                                  list(threads=ncores, pictures=FALSE, warnings = FALSE)))

    # train the model
    model_train <- do.call(maxent, cv.model.args)

    # get the outputs
    model_output <- UsefulPlants::get_maxent_output(outputdir_path[1], j)

    if(varying_parameter){

      # get the previous command line
      old_command_line <- get_maxent_command_line(file.path(outputdir_path[1],"maxent.log"))

      for(i in 3:length(maxent_settings[[varying_parameter_name]])){

        # change output directory
        new_output_directory <- outputdir_path[i]
        new_command_line <- UsefulPlants::set_command_param(old_command_line, "outputdirectory", new_output_directory)

        # change parameter value
        new_command_line <- UsefulPlants::set_command_param(new_command_line, varying_parameter_name, maxent_settings[[varying_parameter_name]][i])

        # run the model
        mem <- if("memory_allocated" %in% maxent_settings) maxent_settings$memory_allocated else 512
        new_command_line <- UsefulPlants::set_command_arg(new_command_line, "density.MaxEnt", paste0('-mx',mem,'m',' -jar ',cv.model.args$path_to_maxent))
        system(new_command_line)

        # get the output
        model_output %<>%
          dplyr::bind_rows(get_maxent_output(new_output_directory, j))
      }

      # set row names
      rownames(model_output) <- paste(varying_parameter_name, maxent_settings[[varying_parameter_name]],sep="_")
    }
    model_output

  }

  return(model_output)

}

#' @export
train_geoModel <- function(loc_dat,
                           bg_dat=NULL,
                           algorithm="idw",
                           coordHeaders=NULL,
                           outputdir,
                           species_name=NULL,
                           grid_res = 0.01666667, save_bg=TRUE, ...){

  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  if(missing(outputdir))
    stop("Output directory is missing. Please select a directory for inverse-distance weighted model output.")

  loc_file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  loc_data_flag = !loc_file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!loc_file_flag & !loc_data_flag)
    stop('Unable to read location data. Please provide valid location data')

  if(loc_file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  if(is.null(species_name)){
    if(loc_file_flag){
      species_name=gsub(".csv","",basename(loc_dat))
      if(nchar(species_name)==0L)
        species_name="Unknown"
    }
    else
      species_name="Unknown"
  }

  if(loc_file_flag)
    loc_dat = read.csv(loc_dat, header=TRUE)

  if(!dir.exists(outputdir))
    stop("The directory ",outputdir," does not exist or is not accessible.")

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }

  sf_loc_data <- loc_dat %>%
    dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
    sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84')) %>%
    sf::st_coordinates()

  #==========================
  #= 1. Make the background #
  #==========================
  if(is.null(bg_dat)){

    bg_args <- list(loc_dat=sf_loc_data, do.alpha_hull=TRUE, dissolve=FALSE, land_file=worldbackground, save.outputs=FALSE, verbose=FALSE)
    dot.args <- list(...)
    if(length(dot.args)>0L){
      arg.names <- names(dot.args)
      if(any(c("land_file") %in% arg.names))
        bg_args[["land_file"]] <- dot.args[["land_file"]]
    }
    domain <- suppressMessages(try(do.call(make_geographic_domain, bg_args), silent=TRUE))

    if(inherits(domain,'try-error'))
      stop("Unable to build a background from the set of points")

    if(save_bg){
      out_bg_dn <- file.path(outdir,"modelsBackground")
      if(!dir.exists(out_bg_dn)) dir.create(out_bg_dn)
      saveRDS(domain, file.path(out_bg_dn, paste0(species_name,".rds")))
    }

    # get background coordinates
    raster_template  <- raster::raster(domain, res = grid_res)
    raster_template[] <- 1
    bg_dat <- raster::xyFromCell(raster_template, cell=raster::extract(x=raster_template, y=as(domain,"Spatial"),cellnumber=TRUE,df=TRUE)[['cell']])
  }

  #======================
  #= 2. Train the model #
  #======================
  geo_model <- switch(algorithm,
            "idw"  <- try(dismo::geoDist(p=sf_loc_data,lonlat=TRUE), silent=TRUE),
            "dist" <- try(dismo::geoIDW(p=sf_loc_data, a=bg_dat), silent=TRUE)
  )

  if(inherits(geo_model,"try-error")){
    cat(model)
    stop("Unable to train the geographic model: ", algorithm)
  }

  #=========================
  #= 3. Evaluate the model #
  #=========================
  if(inherits(model,"InvDistWeightModel"))
    list_args <- list(pr=sf_loc_data, bg=bg_dat)
  else
    list_args <- list(pr=sf_loc_data)

  dot.args <- list(...)
  if(length(dot.args)>0L){
    if("tr" %in% names(dot.args))
      list_args <- append(list_args, dot.args[names(dot.args)=="tr"])
  }
  # leave-one-out cross-validation
  eval_dist_model <- do.call(loocv_geo, list_args)

  # save outputs
  saveRDS(list(model.object=geo_model, model.eval=eval_dist_model), file.path(outputdir, paste0(species_name,".rds")))

  return(eval_dist_model)
}



#' utility function that performs a leave one out cross validation for geographic models
#' @export
loocv_geo <- function(pr, bg, tr=seq(from=1.E-3, to=0.99, by=1.E-2)){

  if(missing(pr))
    stop("Presence data locations are missing.")

  cv_eval_metrics <- vector(mode = "list", length = nrow(pr))

  for(i in 1:nrow(pr)){
    #======================
    #= 1. Train the model #
    #======================
    if(!missing(bg)){
      model <- try(dismo::geoIDW(p=pr[-i,], a=bg), silent=TRUE)
    }
    else{
      model <- try(dismo::geoDist(p=pr[-i,]), silent=TRUE)
    }

    if(inherits(model,"try-error")){
      cat(model)
      stop("Unable to train a model of class:",class(model))
    }


    #=========================
    #= 2. Evaluate the model #
    #=========================
    if(inherits(model,"InvDistWeightModel"))
      eval_geo_model <- dismo::evaluate(model=model, p=cbind(pr[i,1],pr[i,2]), a=bg, tr=tr)
    else
      eval_geo_model <- dismo::evaluate(model=model, p=cbind(pr[i,1],pr[i,2]), tr=tr)

    cv_eval_metrics[[i]] <- data.frame(AUC=eval_geo_model@auc,
                                       TSS=max(eval_geo_model@TPR + eval_geo_model@TNR - 1),
                                       OR10=calc_omission_rate(model, occ_train=pr[-i,], occ_test=data.frame(pr[i,1],pr[i,2]), threshold = 10))
  }

  colMeans(do.call(rbind,cv_eval_metrics))
}


#' utility function that calculates the omission rate from a model for different threshold values
#' @export
calc_omission_rate <- function(model,  occ_train, occ_test, threshold = 10) {

  suit_pred_train <- dismo::predict(model, occ_train)
  suit_pred_test <- dismo::predict(model, occ_test)

  om_rate <- vector("numeric", length = length(threshold))

  for (i in 1:length(threshold)) {
    val <- ceiling(length(occ_train[, 1]) * threshold[i] / 100) + 1

    quantval_or <- sort(suit_pred_train)[val]

    om_rate[i] <- as.numeric(length(suit_pred_test[suit_pred_test < quantval_or]) / length(suit_pred_test))
  }
  return(om_rate)
}


# bg_with_data <- sf::st_coordinates(env_bg) %>%
#   dplyr::bind_cols(map_lfd(env_bg, function(x) raster::extract(as(x,"Raster"),.))) %>%
#   na.omit()

# utility function: map function that transforms input by applying a function to each element and returning
# a data frame with a number of columns equal to the number of element in the input.
#' @export
map_lfd <- function(x, fun, ...){
  num_index <- length(x)
  ret <- lapply(1:num_index, FUN=function(i) fun(`[`(x,i)), ...)
  res <- as.data.frame(do.call(cbind, ret))
  names(res) <- names(x)
  res
}


# utility function: read multiple raster layers from a directory.
#' @export
read_layers <- function(path){

  if(!dir.exists(path))
    stop("Unable to find directory ",path,".")

  # Regular expression to check for rasters layer with different formats
  layer_formats <- "(*.grd$)|(*.asc$)|(*.bil$)|(*.sdat$)|(*.rst$)|(*.tif$)|(*.envi$)|(*.img$)|(*hdr.adf$)"
  list_layers <- list.files(path, pattern = layer_formats, full.names = TRUE, recursive=TRUE)

  if(length(list_layers)==0L)
    stop("No raster layers found in directory ",path)

  layers <- try(stars::read_stars(list_layers),silent=TRUE)

  if(inherits(layers,"try-error"))
    stop("Unable to read raster layers. Please check the dimensions of your layers.")

  return(layers)
}

# utility function: strip raster layers file extension.
#' @export
strip_extension <- function(fn){
  stringr::str_extract(fn, ".+?(?=\\.[a-z]{3}$)")
}




