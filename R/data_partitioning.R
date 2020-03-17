#'                                 Geographic "block" partitioning of occurrence records for evaluation
#'
#' Partition both occurrence records and background into evaluation bins based on some spatial rules.
#'
#'
#' @param kdat A two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param k A numeric integer specifying the number of cluster required.
#' @param bg A RasterLayer object of the background to be partitionned.
#' @param method A character string specifying the method to use to assign background cell to each block among "nearest.cluster" and "nearest.center" (the default).
#' @param grid_res A numeric value specifying the resolution of the output raster layer in decimal degree. Default is 0.16666667 (~ 20 km).
#' @param coastline An optional character string specifying the path to a shapefile delineating the coastline.
#' @param sf A logical. Should a sf object be returned ? Default is TRUE
#' @param ..., Additional paramters to be passed to \code{kmeansEqual}
#' @return A sf object (default) or a RasterStack of the geographically masked blocks of occurrence and background points/cells.
#' @export
make_geographic_block <- function(kdat, k = 3, bg = NULL, method=c("nearest.center"), grid_res = 0.16666667, coastline = NULL, ...){

  #==================================
  # 1. Partition occurrence records #
  #==================================
  if(verbose) cat('> [...partitioning occurrence records in ',k,' groups...]\n')
  dot.args <- list(...)
  if(length(dot.args)>0L){
    dot.args <- dot.args[!names(dot.args)=="plot"]
  }
  kdatAssigned <- do.call(kmeansEqual,
                        c(list(kdat=kdat, k=k, verbose=FALSE, plot=FALSE), dot.args))

  #==========================================
  # 2. Build convex hulls around each group #
  #==========================================
  proj <- if(is.null(bg) || is.na(sf::st_crs(bg))) '+proj=longlat +datum=WGS84' else sf::st_crs(bg)
  sf_loc_data <- kdatAssigned$Data %>%
    sf::st_as_sf(coords=c('lon', 'lat'), crs=sf::st_crs(bg))
  # create the first convex hull
  j = 1
  sf_loc_grp <- sf_loc_data %>%
    filter(assigned==j)

  sf_block <- sf::st_as_sf(data.frame(kgroup=j, geometry=sf::st_convex_hull(sf::st_union(sf_loc_grp))))

  # add other hulls
  j = j + 1
  while(j <= k){
    sf_loc_grp <- sf_loc_data %>%
      filter(assigned==j)

    sf_block <- sf_block %>%
      rbind(sf::st_sf(data.frame(kgroup=j, geometry=sf::st_convex_hull(sf::st_union(sf_loc_grp)))))
    j = j + 1
  }
  #=============================
  # 3. Build disjoint polygons #
  #=============================
  sf_block_disjoint <- sf_block
  # compute intersections between all polygons and retrieve polygons of the intersections
  sf_block_intersect <- sf::st_intersection(sf_block) %>%
    filter(n.overlaps>1) %>%
    select(kgroup)
  # substract polygons intersection to the blocks to obtain block disjoints
  if(nrow(sf_block_intersect)>0L){
    sf_block_disjoint <- sf_block_intersect %>%
      st_difference(sf_block, .) %>%
      select(kgroup)
  }

  #=========================
  # 4. Make the background #
  #=========================
  if(is.null(bg)){

    if(is.null(coastline))
      coastline = worldbackground

    domain <- try(make_geographic_domain(loc_dat=kdat, do.alpha_hull=TRUE, dissolve=FALSE, verbose=FALSE, land_file=coastline), silent=TRUE)

    if(inherits(domain,'try-error'))
      stop("Unable to build a background from the set of points")

    # fasterize the domain
    raster_template  <- raster::raster(domain, res = grid_res)
    bg <- fasterize::fasterize(sf::st_cast(domain,"MULTIPOLYGON"), raster_template) - 1 # -1 to ensure background is a distinct group
  }

  #----------------------------------------------------
  #= 4.a Assign background cells to the block cluster #
  #----------------------------------------------------
  bg_fasterized <- raster::mask(fasterize::fasterize(sf_block_disjoint, bg, field="kgroup", background=0), bg)

  # assign occurrence records cells outside of the block to their corresponding group
  # if(nrow(sf_block_intersect)>0L){
  #   cell_occ <- kdatAssigned$Data %>%
  #     sf::st_as_sf(coords=c('lon', 'lat'), crs=sf::st_crs(bg)) %>%
  #     sf::st_intersects(., sf_block_disjoint) %>%
  #     purrr::map_lgl(., function(x) length(x) == 0L ) %>%
  #     dplyr::filter(kdatAssigned$Data, .)
  # }
  tbl <- data.table::data.table(cell=which(values(bg_fasterized)>=0), key="cell")
  tbl[,`:=`(Group = as.character(bg_fasterized[cell]),
            Lon = raster::xFromCell(bg_fasterized, cell),
            Lat = raster::yFromCell(bg_fasterized, cell))]
  data.table::setkey(tbl, 'Group')

  switch(method,
       "nearest.center" = {
         centers <- kdatAssigned$Centers
         tbl[, kgroup := ifelse(Group=='0', as.character(which.min(geosphere::distGeo(cbind(Lon, Lat)), centers)), Group), by=seq_len(NROW(tbl))]
         bg_fasterized[tbl['0'][['cell']]] <- tbl['0'][['kgroup']]
       },
       "nearest.cluster" = {
         kgroup <- tbl['0'][, .(Lon,Lat)] %>%
           sf::st_as_sf(coords=c('Lon', 'Lat'), crs=sf::st_crs(bg)) %>%
           sf::st_distance(., sf_block_disjoint) %>%
           min.col() %>% as.character()
         tbl[, kgroup := replace(Group, Group=='0', kgroup)]
         bg_fasterized[tbl['0'][['cell']]] <- as.numeric(tbl['0'][['kgroup']])
       }
  )


  #-------------------------------
  #= 4.b Create background masks #
  #-------------------------------
  data.table::setkey(tbl, 'kgroup')
  comb <- combn(1:k,k-1)

  if(!sf){
    bg_masks <- stack()
    for(j in 1:k){
      blcks <- comb[,j]
      blck <- setdiff(1:k,blcks)
      cells_group <- c(tbl[as.character(blcks[1])][['cell']], tbl[as.character(blcks[2])][['cell']])
      non_cells_group <- tbl[as.character(blck)][['cell']]
      bg_group <- bg_fasterized
      bg_group[cells_group] <- "1"
      bg_group[non_cells_group] <- "0"
      bg_masks <- raster::addLayer(bg_masks, bg_group)
      names(bg_masks) <- apply(comb,2,function(x) paste(letters[x],collapse="_"))
    }
  }else{
    bg_masks <- bg_fasterized %>%
      stars::st_as_stars() %>%
      sf::st_as_sf(as_points=FALSE, merge=TRUE) %>%
      sf::st_cast("MULTIPOLYGON")
  }
  # return multi-layer background mask
  return(bg_masks)
}

#'                                 Kmeans clustering with near equal group size
#'
#' Modified from https://rviews.rstudio.com/2019/06/13/equal-size-kmeans/
#'
#'
#' @param kdat A two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param k A numeric integer specifying the number of cluster required.
#' @param iter_max A numeric integer specifying the maximum number of iterations to run.
#' @return A list with two elements:
#'         \code{Data} A two-column matrix or data.frame or data.table with longitude and latitude coordinates and a column 'assigned' specifying
#'         which of the k cluster each point belongs to.
#'         \code{Centers} A two-column data.frame with longitude and latitude coordinates of the cluster centers.
#' @export
kmeansEqual <- function(kdat, k = 3, iter_max = 30, random.seed=1234, method=c("iqr.outliers","top.outliers"), nstart=20, n.outliers=10, plot=TRUE){

  #---------------------------------------------------------------------------
  #=step 0: Check input data
  #---------------------------------------------------------------------------
  if(!any(inherits(kdat, c("matrix","data.frame","data.table"))))
    stop("Argument kdat must be a matrix, a data.frame or a data.table")
  if(nrow(kdat) < 2)
    stop("Input data must have at least 2 coordinates.")

  if(!is.numeric(k))
    stop("Argument 'k' must be numeric.")
  if(k<0)
    stop("The number of cluster must be > 0.")
  if(k<2)
    warning("The number of group requested is < 2.")

  if(!is.numeric(iter_max))
    stop("Argument 'iter_max' must be numeric")
  if(iter_max<0)
    stop("The number of iterations must be > 0.")

  #---------------------------------------------------------------------------
  #=step 1: Initialize "naive" clusters with a basic kmeans
  #---------------------------------------------------------------------------

  set.seed(random.seed)

  kclust <- kdat %>%
    stats::kmeans(k, nstart = nstart)

  # initial clusters
  centers = kclust$centers

  # retrieve coordinates if needed
  if(!all(c("lon","lat") %in% colnames(kdat))){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(kdat))[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(kdat))[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
    kdat %>%
      dplyr::rename(lon = coordHeaders[1], lat = coordHeaders[2])
  }
  if(is.matrix(kdat))
    kdat <- data.frame(kdat)

  # start iterations
  iter = 0
  while(iter <= iter_max){

    kclust <- kdat %>%
      select(lon,lat) %>%
      stats::kmeans(centers)

    kdat$assigned = kclust$cluster
    #---------------------------------------------------------------------------
    #=step 2: Compute the distance matrix between each point and the centroids
    #---------------------------------------------------------------------------
    # calculate the distance to the center
    for(j in 1:nrow(centers)){
      kdat <- kdat %<>%
        dplyr::mutate(kdist = geosphere::distGeo(cbind(lon,lat), centers[j,]))
      # rename the column
      kdat <- kdat %>%
        dplyr::rename(!!paste0('kdist',j) := kdist)
    }
    if(iter > 0){
      # assign the outliers from each group to the closest cluster center
      switch(method[1],

             "top.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::arrange_at(paste0("kdist", j), desc) %>%
                 dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(row_number() <= n.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::arrange_at(paste0("kdist", j), desc) %>%
                           dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(row_number() <= n.outliers, is_kdist_min, assigned))
                   )
               }
             },

             "iqr.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate_at(.vars= vars(!!paste0("kdist", j)),
                                  .funs= list(is.outliers = ~ . > quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                 dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate_at(.vars= vars(!!paste0("kdist", j)),
                                            .funs= list(is.outliers = ~ . > quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                           dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             }
      )
      kdat.copy %<>% data.frame(row.names=NULL)
      kdat <- kdat.copy

    }
    #---------------------------------------------------------------------------
    #=step 3: Re-assign points to clusters
    #---------------------------------------------------------------------------
    kdat$index = 1:nrow(kdat)
    working = kdat
    ReAssignment = nrow(kdat) - (nrow(kdat) %% k)

    for(i in 1:ReAssignment){
      #cluster counts can be off by 1 due to uneven multiples of k.
      j = if(i %% k == 0) k else (i %% k)
      itemloc =  working$index[which.min(working[,(paste0("kdist", j))])][1]
      kdat$assigned[kdat$index == itemloc] = j
      working %<>%
        dplyr::filter(!index == itemloc)
    }
    # if there are some leftover points, assign to whoever's closest, without regard to k
    if(length(working$index)>0L){
      for(i in working$index){
        kdat$assigned[kdat$index == i] = which.min(working[working$index == i,grepl("kdist",colnames(working))])
      }
    }
    # assign the outliers from each group to the closest cluster center
    switch(method[1],

           "top.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::arrange_at(paste0("kdist", j), desc) %>%
               dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(row_number() <= n.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::arrange_at(paste0("kdist", j), desc) %>%
                         dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(row_number() <= n.outliers, is_kdist_min, assigned))
                 )
             }
           },

           "iqr.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate_at(.vars= vars(!!paste0("kdist", j)),
                                .funs= list(is.outliers = ~ . > quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
               dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate_at(.vars= vars(!!paste0("kdist", j)),
                                          .funs= list(is.outliers = ~ . > quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                         dplyr::mutate(is_kdist_min =  select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           }
    )
    kdat.copy %<>% data.frame(row.names=NULL)
    kdat <- kdat.copy
    #---------------------------------------------------------------------------
    #=step 4: Recalculate the centroids
    #---------------------------------------------------------------------------
    j = 1
    NewCenters <- kdat %>%
      dplyr::filter(assigned == j) %>%
      dplyr::select(lon, lat) %>%
      kmeans(1) %$% centers
    while(j < k){
      j = j+1
      NewCenters %<>%
        rbind(kdat %>%
                dplyr::filter(assigned == j) %>%
                dplyr::select(lon, lat) %>%
                kmeans(1) %$% centers)
    }
    # New centroids
    NewCenters %<>% data.frame(row.names=NULL)
    centers <- NewCenters

    # keep assigments lon, lat and assigments only
    kdat <- kdat %>%
      dplyr::select(lon, lat, assigned)

    if(plot){
      if(!is.factor(kdat$assigned))
        kdat$assigned %<>% as.factor()
      print(
        kdat %>% ggplot2::ggplot(aes(x = lon, y = lat, color = assigned)) +
          theme_minimal() + geom_point()  +
          geom_point(data = centers, aes(x = lon, y = lat), color = "black", size = 4)
      )
    }

    iter = iter + 1
  }
  # end iterations

  # assignment to factor
  if(!is.factor(kdat$assigned))
    kdat$assigned %<>% as.factor()

  return(list(Data=kdat, Centers=centers))
}

# A helper function that erases all of y from x:
# st_erase = function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y)))

# utility function: find the minimum position for each row of a matrix, breaking ties at random.
#' @export
min.col <- function(m, ...) max.col(-m,...)
