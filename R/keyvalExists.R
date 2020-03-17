#' Check for a key value whether or not a key value exist
#'
#' Check whether or not a key value exist in a data.table
#'
#' @param DT A data.table object. DT must have a key set.
#' @param keyval A character/numeric vector of key values.
#' @param how A character string specifying the type of logical output to return comparison to perform. Either \code{'all'} or \code{'any'}
#' @return A logical value if \code{how} is set to 'all' or 'any'. If \code{how} is set to NULL (default), a logical vector is returned. TRUE if the key value exists in DT, otherwise FALSE.
#' @export
# https://stackoverflow.com/questions/17331684/fast-exists-in-data-table
is.keyval.exists <- function(DT, keyval, how=NULL){

  if(!inherits(DT,"data.table"))
    stop("Argument DT must be a data.table object")
  if(!data.table::haskey(DT))
    stop("DT must have a key !")
  # get the key
  keyDT <- data.table::key(DT)

  if(is.null(how))
    return(keyval %in% .subset2(DT[J(keyval), mult="first", nomatch=0], keyDT))

  switch(how,
         "any" = any(keyval %in% .subset2(DT[J(keyval), mult="first", nomatch=0], keyDT)),
         "all" = all(keyval %in% .subset2(DT[J(keyval), mult="first", nomatch=0], keyDT))
  )
}
