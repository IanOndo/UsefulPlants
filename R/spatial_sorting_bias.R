
#' @export
ssb <- function(p, bg, reference, lonlat=TRUE, avg=TRUE) {

  if(!require(magrittr))
    stop("Package 'magrittr' must be installed.")

  #=============================================
  #= Calculate the pairwise distance matrices
  #=============================================
  pwd_p <- distFast(reference, p)
  pwd_bg <- distFast(reference, bg)

  #---------------------------------------------
  #= Get the distance from testing/background
  #  points to the nearest training presence
  #---------------------------------------------
  Dp <- pwd_p[ , .SD[which.min(distance)], by=orig][['distance']]
  Da <- pwd_bg[ , .SD[which.min(distance)], by=orig][['distance']]

  rm(pwd_p)
  rm(pwd_bg)
  gc()

  if(avg){
    # Get the mean distance from testing/background points to the nearest
    # training presence
    Dp %<>%
      mean()

    Da %<>%
      mean()
    return( stats::setNames(c(Dp, Da),c('Dp','Da')) )

  }else{
    return( list(Dp, Da) )
  }

}


#' @export
distFast <- function(X, Y=NULL, nthreads=NULL, lonlat=TRUE, max_size = 2^31 - 1, chunk_size=1000){

  if(!require(dplyr))
    stop("[distFast] requires package 'dplyr' to work.")

  if(!inherits(X,c("data.frame","data.table")))
    stop("Argument X must be a data.frame or a data.table")

  if(inherits(X,'data.frame'))
    data.table::setDT(X)

  if(!is.null(Y)){
    if(!inherits(Y,c("data.frame","data.table")))
      stop("Argument Y must be a data.frame or a data.table")
    if(inherits(Y,'data.frame'))
      data.table::setDT(Y)
  }

  if(is.null(nthreads))
    nthreads <- parallel::detectCores()-1
  else
    nthreads <- min(nthreads, parallel::detectCores())

  #====================================================
  #= 1. Create a join table with possible combinations
  #     of origin and destination
  #====================================================
  id <- as.character(1:nrow(X))

  if(!is.null(Y)){
    id.y <- as.character((nrow(X)+1):(nrow(X)+nrow(Y)))
    W <- optiRum::CJ.dt(data.table::data.table(id),data.table::data.table(id.y)) %>%
      data.table::setnames(c("orig","dest")) %>%
      data.table::setkeyv(c("orig","dest"))
    id <- as.character(1:(nrow(X)+nrow(Y)))
  }
  else{
    nComb <- RcppAlgos::comboCount(id, m = 2)

    if(nComb > max_size){
      ids <- split2(id, chunk_size)

      require(doParallel)
      doParallel::registerDoParallel(nthreads)
      W <- foreach::foreach(k = 1:length(ids), .packages=c('data.table','RcppAlgos','dplyr'), .combine=rbind) %dopar% {

        gc()

        W <- data.table::data.table(RcppAlgos::comboGeneral(ids[[k]], m = 2, freqs=,nThreads = nthreads-1)) %>%
          data.table::setnames(c("orig","dest")) %>%
          data.table::setkeyv(c("orig","dest"))
        W
      }
      doParallel::stopImplicitCluster()

    }else{
      W <- data.table::data.table(RcppAlgos::comboGeneral(id, m = 2, nThreads = nthreads)) %>%
        data.table::setnames(c("orig","dest")) %>%
        data.table::setkeyv(c("orig","dest"))
    }
  }

  #====================================================
  #= 2. Add the same id column to the original dataset
  #     and use it as a key
  #====================================================
  if(!is.null(Y)){
    X <- X %>%
      data.table::setnames(colnames(Y)) %>%
      dplyr::bind_rows(Y)

    X[, orig:= id]
    data.table::setkey(X, "orig")
  }else{
    X[, orig:= id]
    data.table::setkey(X, "orig")
  }

  #====================================================
  #= 3. Combine the two datasets alternating inner join
  #     between the origin and the destination field
  #====================================================
  Z <- W[X, nomatch=0]
  data.table::setkey(Z, "dest")
  Z <- Z[X, nomatch=0]
  data.table::setkey(Z, "orig")

  #====================================================
  #= 4. Compute the distance between pairs of locations
  #
  #====================================================
  distPlane <- function (x, y) {
    dfun <- function(x, y) {
      sqrt( (x[,1] - y[,1])^2 + (x[,2] - y[,2])^2 )
    }
    n = nrow(x)
    m = nrow(y)
    dm = matrix(ncol = m, nrow = n)
    for (i in 1:n) {
      dm[i, ] = dfun(x[i, ,drop=FALSE], y)
    }
    return(dm)
  }

  data.table::setnames(Z, 3:6, c('x', 'y', 'i.x', 'i.y'))
  out <- try(
    if(lonlat)
      Z[, distance := spatialrisk::haversine(y, x, i.y, i.x)]
    else
      Z[, distance := distPlane(c(x, y), c(i.x, i.y))]
    , silent=TRUE)

  if(inherits(out,"try-error")){
    cat(out)
    stop("Unable to compute the distance matrix.")
  }

  return(Z)
}




#' helper function that split a vector x in chunks of maximum size d
#' @export
split2 <- function(x,d) split(x, ceiling(seq_along(x)/d))

