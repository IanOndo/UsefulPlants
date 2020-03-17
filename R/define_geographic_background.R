#'                                 Estimating the potential geographic distribution area of a species
#'
#' Delineate a geographic distribution area of a species by combining a set occurrence records locations with bioregions
#'
#' @param loc_dat A two-column matrix or data.frame or data.table or a path to a csv file with longitude and latitude coordinates.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param output_dir A directory where to save the domain defined.
#' @param do.alpha_hull A logical. Should an alpha hull be built from the set of points ? Default is \code{TRUE}.
#' @param dissolve A logical. Should the borders between adjacent polygons be dissolved ? Default is \code{FALSE}.
#' @param path_to_alpha_hull A character string specifying the path to the alpha-hull shapefile. Ignored if \code{do.alpha_hull=TRUE}.
#' @param ... Parameters to be passed to \code{make_alpha_hulls} function
#' @return An sf object
#' @export
make_geographic_domain <- function(loc_dat, coordHeaders=NULL, output_dir=NULL, do.alpha_hull=TRUE, dissolve=FALSE, path_to_alpha_hull=NULL, verbose=TRUE, ...) {
  #===================
  #= 0.Check packages
  #===================
  # check if required packages are available
  pkgs.required = c('dplyr','sf','magrittr','purrr','AlphaHullRangeModeller')
  pkgs.available = sapply(pkgs.required, require, quietly = TRUE, character.only=TRUE)
  if(!any(pkgs.available)){
    warning(paste('Packages',paste(pkgs.required[!pkgs.available],collapse=', '),'failed to load'))
    stop("Try to re-install packages that failed to load")
  }
  if(verbose){
    cat('#==================\n')
    cat('#= 1. Check inputs\n')
    cat('#==================\n\n')
  }
  if(missing(loc_dat))
    stop("Please provide a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  file_flag = tryCatch(file.exists(loc_dat), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  data_flag = !file_flag & any(inherits(loc_dat, c("data.frame","matrix")))

  if(!file_flag & !data_flag)
    stop('Unable to read input data. Please provide valid input data')

  if(file_flag && length(readLines(loc_dat))==0L)
    stop(paste("The file provided:",loc_dat," is empty!"))

  # ensure that the output directory provided exists
  if(!is.null(output_dir) && !dir.exists(output_dir)){
    warning(paste("Output directory:", output_dir,"does not exist or could not be found."));flush.console()
  }

  if(!is.null(path_to_alpha_hull)){
    if(!file.exists(path_to_alpha_hull)){
      warning(paste("Path:", path_to_alpha_hull,"does not exist or could not be found."));flush.console()
      do.alpha_hull =TRUE
    }
  }

  cat('\n')
  if(verbose){
    cat('#==================\n')
    cat('#= 2. Make domain \n')
    cat('#==================\n\n')
  }
  require(magrittr)

  if(do.alpha_hull){
    if(verbose) cat('> [...building alpha-hull...]\n')
    dot.args <- list(...)
    if(length(dot.args)>0L){
      dot.args <-dot.args[!names(dot.args)=="save.outputs"]
    }
    alpha_hull <- do.call(AlphaHullRangeModeller::make_alpha_hulls,
                          c(list(loc_data=loc_dat, save.outputs=FALSE, verbose=FALSE), dot.args))
    if(is.null(alpha_hull))
      stop("Unable to build an alpha hull from this set of points.")
  }
  else{
    if(verbose) cat('> [...read alpha-hull from file...]\n')
    alpha_hull <- try(suppressMessages(sf::st_read(path_to_alpha_hull)), silent=TRUE)
    if(inherits(alpha_hull,'try-error'))
      stop("Unable to read alpha hull from this file specified.")
  }
  if(verbose) cat('> [...Done...]')
  cat('\n\n')
  if(verbose){
    cat('#----------------------------------------------\n')
    cat('#= 2.a. Step 1: Coarse scale domain estimation \n')
    cat('#----------------------------------------------\n')
  }

  if(file_flag){
    loc_dat <- try(read.csv(loc_dat,h=TRUE),silent=TRUE)
    if(inherits(loc_dat,"try-error"))
      stop("Unable to read the csv file of the input dataset.")
  }

  sf_loc_data <- loc_dat %>%
    as.data.frame(row.names=NULL)

  if(is.null(coordHeaders)){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(sf_loc_data), value=TRUE)[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(sf_loc_data), value=TRUE)[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
  }else if(length(coordHeaders)!=2){
    stop("Argument 'coordHeaders' must be of length 2.")
  }
  sf_loc_data <- sf_loc_data %>%
      dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
      sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs(alpha_hull))

  if(verbose) cat('> [...keep biomes with more than one occurrence records inside...]\n')
  num.point.by.biome <- suppressMessages(lengths(sf::st_intersects(biomes, sf_loc_data)))
  # keep biomes with more than one occurrence records inside
  biomes.with.more.than.one.point <- biomes[num.point.by.biome > 1,]
  if(verbose) cat('> [...extract polygons whose biomes contain more than one record and that intersect with the alpha hull...]\n')
  # extract polygons whose biomes contains more than one records and that intersect with the alpha hull (https://rpubs.com/sogletr/sf-ops)
  biomes_more_than_one_point_geom <- biomes.with.more.than.one.point %>%
    sf::st_cast(., "POLYGON")
  biome_alpha_hull <- sf::st_geometry(alpha_hull)
  cond <- suppressMessages(sf::st_intersects(biomes_more_than_one_point_geom, biome_alpha_hull) %>%
                      purrr::map_lgl(., function(x) length(x) == 1))
  occupied_polygons <- biomes_more_than_one_point_geom %>%
    dplyr::filter(cond) %>%
    dplyr::select(BIOME_NAME) %>%
    dplyr::rename(REGION_NAME=BIOME_NAME)

  geographic_domain <- occupied_polygons

  if(verbose) cat('> [...Step 1 [Done]...]\n')

  is.outside.geographic_domain <- suppressMessages(sf::st_intersects(sf_loc_data, geographic_domain) %>%
                                                 purrr::map_lgl(., function(x) length(x) == 0L ))
  if(sum(is.outside.geographic_domain)==0L){
    cat('\n')
    if(dissolve){
      if(verbose){
        cat('#--------------------------------------------\n')
        cat('#= 2.b Step 2: Merge estimated domain \n')
        cat('#--------------------------------------------\n')
      }
      suppressMessages(
        geographic_domain <- geographic_domain %>%
          sf::st_buffer(0) %>%
          sf::st_union()
      )
      gc()
    }
    return(geographic_domain)
  }

  cat('\n\n')
  if(verbose){
    cat('#--------------------------------------------\n')
    cat('#= 2.b. Step 2: Fine scale domain estimation \n')
    cat('#--------------------------------------------\n')
  }
  # create working location data
  working_loc_data <- sf_loc_data %>%
    filter(is.outside.geographic_domain)

  if(verbose) cat("> [...identify unique records in biomes...]\n")
  # identify unique records in biomes
  biomes.with.one.point <- biomes[num.point.by.biome == 1,]
  is.unique.point.in.biome = suppressMessages(sf::st_intersects(working_loc_data, biomes.with.one.point) %>%
                                               purrr::map_lgl(., function(x) length(x) > 0L))
  if(any(is.unique.point.in.biome)){
    # retrieve ecoregions of records flagged as unique within the biomes
    if(verbose) cat("> [...retrieve ecoregions of records flagged as unique within the biomes...]\n")
    ecoreg_from_single_record <- suppressMessages(sf::st_intersection(ecoregions, working_loc_data[is.unique.point.in.biome,]) %>%
                                                    dplyr::pull(ECO_NAME) %>%
                                                    match(., ecoregions$ECO_NAME) %>%
                                                    ecoregions[.,]) %>%
                                                    dplyr::select(ECO_NAME) %>%
                                                    dplyr::rename(REGION_NAME=ECO_NAME)

    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_single_record)
    )

    # remode those records
    working_loc_data %<>%
      dplyr::filter(!is.unique.point.in.biome)

    if(lengths(working_loc_data)==0L){
      cat('\n\n')
      if(verbose) cat('> [...Step 2 [Done]...]\n')
      if(dissolve){
        if(verbose){
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]\n')
      }
      return(geographic_domain)
    }

  }else{
    if(verbose) cat("> [...No records found...]\n")
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if(verbose) cat("> [...identify records excluded from the alpha hull...]\n")
  # identify records excluded from the alpha hull
  is.outside.alpha.hull=suppressMessages(sf::st_intersects(working_loc_data, biome_alpha_hull) %>%
                                              purrr::map_lgl(., function(x) length(x) == 0))

  if(any(is.outside.alpha.hull)){
    # do the same with records excluded from the alpha hull
    if(verbose) cat("> [...retrieve ecoregions of records flagged outside of the alpha hull...]\n")
    # retrieve ecoregions of records flagged as outside of the alpha hull
    ecoreg_from_alpha_hull_exclusion <- suppressMessages(sf::st_intersection(ecoregions, working_loc_data[is.outside.alpha.hull,]) %>%
                                                           dplyr::pull(ECO_NAME) %>%
                                                           match(., ecoregions$ECO_NAME) %>%
                                                           ecoregions[.,]) %>%
                                                           dplyr::select(ECO_NAME) %>%
                                                           dplyr::rename(REGION_NAME=ECO_NAME)
    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_alpha_hull_exclusion)
    )

    # remode those records
    working_loc_data %<>%
      dplyr::filter(!is.outside.alpha.hull)

    if(lengths(working_loc_data)==0L){
      cat('\n\n')
      if(verbose) cat('> [...Step 2 [Done]...]\n')
      if(dissolve){
        if(verbose){
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        if(verbose) cat('> [...Step 3 [Done]...]')
      }
        return(geographic_domain)
    }

  }else{
    if(verbose) cat("> [...No records found...]\n")
  }

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if(verbose) cat("> [...identify geographic outliers...]\n")
  # identify geographic outliers
  is.geographic.outlier = sf::st_coordinates(working_loc_data) %>%
      data.frame(row.names=NULL) %>%
      dplyr::mutate(kdist = geosphere::distGeo(cbind(X,Y), colMeans(cbind(X,Y)))) %>%
      dplyr::mutate_at(.vars= vars(kdist), .funs= list(is.outlier = ~ . > quantile(., probs=.75) + 1.5 * stats::IQR(.))) %$% is.outlier

  if(any(is.geographic.outlier)){
    if(verbose) cat("> [...retrieve ecoregions of records flagged as geographic outliers...]\n")
    # retrieve ecoregions of records flagged as unique within the biomes
    ecoreg_from_outliers <- suppressMessages(sf::st_intersection(ecoregions, working_loc_data[is.geographic.outlier,]) %>%
                                               dplyr::pull(ECO_NAME) %>%
                                               match(., ecoregions$ECO_NAME) %>%
                                               ecoregions[.,]) %>%
                                               dplyr::select(ECO_NAME) %>%
                                               dplyr::rename(REGION_NAME=ECO_NAME)
    # add ecoregions to the domain
    suppressMessages(
      geographic_domain <- geographic_domain %>%
        rbind(ecoreg_from_outliers)
    )

    # remode those records
    working_loc_data %<>%
      dplyr::filter(!is.geographic.outlier)

    if(lengths(working_loc_data)==0L){
      cat('\n\n')
      if(verbose) cat('> [...Step 2 [Done]...]\n')
      if(dissolve){
        if(verbose){
          cat('#--------------------------------------------\n')
          cat('#= 2.c Step 3: Merge estimated domain \n')
          cat('#--------------------------------------------\n')
        }
        suppressMessages(
          geographic_domain <- geographic_domain %>%
            sf::st_buffer(0) %>%
            sf::st_union()
        )
        gc()
        if(verbose) cat('> [...Step 3 [Done]...]\n')
      }
      return(geographic_domain)
    }
  }else{
    if(verbose) cat("> [...No records found...]\n")
    warning('[WARNING] ',lengths(working_loc_data),' records remaining outside of the geographic range !')
  }

  return(geographic_domain)
}

#'                                 Extent the bounding box of a sf point object
#' https://www.jla-data.net/eng/adjusting-bounding-box-of-a-tmap-map/
#' @param loc_dat A sf object or a two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param coordHeaders A character string vector of length two giving the names of the coordinates in \code{loc_dat}
#' @param bbox_ext A numeric value specifying how much the bounding box should be extended in percentage.
#' @export
make_bbox_domain <- function(loc_dat, coordHeaders=NULL, bbox_ext=0.5, min_xrange=0.1666667, min_yrange=0.1666667, left=TRUE, right=TRUE, bottom=TRUE, top=TRUE){

  require(magrittr)

  if(missing(loc_dat))
    stop("Occurrence date are missing.")

  if(!inherits(loc_dat, "sf")){
    if(is.null(coordHeaders)){
      id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(loc_dat), value=TRUE)[1]
      id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(loc_dat), value=TRUE)[1]
      coordHeaders <- c(id_x_lon, id_y_lat)
    }else if(length(coordHeaders)!=2 || all(!coordHeaders %in% names(loc_dat))){
      stop("coordHeaders must be a vector of length 2 of names from loc_dat.")
    }
    loc_dat %<>%
      dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
      sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs('+proj=longlat +datum=WGS84'))
  }

  bbox_new <- sf::st_bbox(loc_dat) # current bounding box

  xrange <- max(bbox_new$xmax - bbox_new$xmin, min_xrange)# range of x values
  yrange <- max(bbox_new$ymax - bbox_new$ymin, min_yrange)# range of y values

  if(left)
    bbox_new[1] <- bbox_new[1] - (bbox_ext * xrange) # xmin - left
  if(right)
    bbox_new[3] <- bbox_new[3] + (bbox_ext * xrange) # xmax - right
  if(bottom)
    bbox_new[2] <- bbox_new[2] - (bbox_ext * yrange) # ymin - bottom
  if(top)
    bbox_new[4] <- bbox_new[4] + (bbox_ext * yrange) # ymax - top

  bbox_new <- bbox_new %>%  # take the bounding box ...
    sf::st_as_sfc() # ... and make it a sf polygon

  return(bbox_new)
}
