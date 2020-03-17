#'                                                     Project a fitted Maxent model
#'
#' This function is a wrapper for the function `project` from the `rmaxent` package
#' It projects a fitted Maxent model to new environmental data enabling for post-processing of maxent predictions.
#'
#' @param lambdas Either (1) a `MaxEnt` fitted model object (fitted with
#'   the `maxent` function in the `dismo` package), (2) a file path to
#'   a Maxent .lambdas file, or (3) a `lambdas` object returned by
#'   [parse_lambdas()].
#' @param newdata A `RasterStack`, `RasterBrick`, `list`,
#'   `data.frame`, `data.table`, or `matrix` that has
#'   layers/elements/columns whose names correspond to the names of predictors
#'   used to fit the model. These layers/elements/columns must all have the same
#'   length.
#' @param training_data [Optional] Same as newdata. Only required for post-processing options 2 and 3
#' @param outType A integer number between 0-3 specifying which output is to be calculated and returned. Default is 0 ('raw', 'logistic' and 'cloglog' not adjusted).
#'                A number between 1-3 involved post-processing: 1 = raw adjusted, 2 = cloglog adjusted, 3 = logistic adjusted.
#' @param estim_prevalence A numeric value estimating the prevalence of the species.
#' @param return_lfx Logical. Should `Raster` layers be returned giving
#'   lambda*feature values for each feature with a non-zero lambda? Currently
#'   ignored if `newdata` is not a `Raster*` object.
#' @param mask (Optional; requires that `newdata` is a `Raster*`
#'   object.) A `Raster` object with `NA` values in cells for which
#'   the model should _not_ be projected. These cells will be assigned
#'   `NA` in the returned output.
#' @param quiet Logical. Should projection progress be reported?
#' @return If `newdata` is a `RasterStack` or `RasterBrick`, a
#'   list with three elements:
#' * `prediction_raw`: a `Raster` layer giving the raw Maxent
#'   prediction;
#' * `prediction_logistic`: a `Raster` layer giving the
#'   logistic Maxent prediction; and
#' * `prediction_cloglog`: a `Raster` layer giving the
#'   cloglog Maxent prediction.
#' If `newdata` is _not_ a `RasterStack` or `RasterBrick`,
#' the raster layers will be replaced with `data.table`s in the returned
#' list.
#' @export
project_maxent <- function(lambdas, newdata, training_data, outType=0, estim_prevalence, maxIter=10, return_lfx=FALSE, mask, quiet=FALSE){

  if(!inherits(lambdas,"lambdas"))
    lambdas <- rmaxent::parse_lambdas(lambdas)

  pred_rmaxent <- rmaxent::project(lambdas, newdata, return_lfx, mask, quiet)

  if(!is.numeric(outType))
    stop("Argument 'outType' must be numeric.")
  if(!outType %>% dplyr::between(0,3))
    stop("Argument 'outType' must be between 0 and 3.")

  is_raster <- inherits(newdata,c("RasterStack","RasterBrick"))

  if(outType>0){

    if(missing(estim_prevalence))
      stop("Argument 'estim_prevalence' should be provided when post-processing is required.")

    if(!is.numeric(estim_prevalence) || is.na(estim_prevalence))
      stop("estim_prevalence must be numeric.")

    if(!estim_prevalence %>% dplyr::between(0,1))
      stop("estim_prevalence must be between 0 and 1.")

    pred_rmaxent <- switch(outType,

                          {
                            # raw adjusted
                              raw_pred  <- pred_rmaxent$prediction_raw
                              num_bg    <- lbds$numBackgroundPoints
                              new_pred  <- raw_pred * estim_prevalence * num_bg
                              new_pred[new_pred>1] <- 1
                              new_pred
                          },

                          {
                            # cloglog adjusted (find 'c' so that prevalence over predictions equals estim_prevalence )
                            if(missing(training_data))
                              stop("training_data should be provided for the cloglog post-processing type")

                            eta <- log(rmaxent::project(lambdas, training_data, return_lfx=FALSE, mask, quiet=TRUE)$prediction_raw)# linear predictor from Maxent (training data)

                            # Finding c by direct search to match observed prevalence (method proposed by Steven Phillips)
                            if(is_raster){
                              eta <- na.omit(values(eta))
                            }
                            cL <- ifelse( psi.estim(c = -5, eta = eta) > estim_prevalence, -20, -5) # select appropriate starting points for search
                            cH <- ifelse( psi.estim(c = 5, eta = eta) < estim_prevalence, 20, 5)
                            c.estim <- search.c(cL=cL, cH=cH, prev=estim_prevalence, niter=maxIter, eta=eta)  # find best c

                            # ALTERNATIVELY, find c as the intercept of a glm model with cloglog link and 'eta' as offset (method proposed by Bob O'Hara)
                            # c.estim <- coef(glm(M[[k]][[r]]$pa ~ 1,  offset = eta,  family=binomial(link="cloglog")))

                            etaT <- log(pred_rmaxent$prediction_raw) # linear predictor from Maxent (test data)

                            1-exp(-exp(etaT+c.estim))
                          },

                          {
                            if(missing(training_data))
                              stop("training_data should be provided for the logistic post-processing type")

                            psi.logi <- rmaxent::project(lambdas, training_data, return_lfx=FALSE, mask, quiet=TRUE)$prediction_logistic # get default logistic output for training data (tau=0.5)
                            if(is_raster){
                              psi.logi <- na.omit(values(psi.logi))
                            }
                            tau.estim <- search.tau(tauL=0, tauH=1, prev=estim_prevalence, niter=10, psi.logi=psi.logi) # find best tau
                            psi.logi.T <- pred_rmaxent$prediction_logistic # get default logistic output for test data (tau=0.5)

                            1/(1+((1-tau.estim)/tau.estim)*((1-psi.logi.T)/psi.logi.T))  # map to chosen tau

                          })


  }

  return(pred_rmaxent)
}


#'                                          Post-processing of maxent raw outputs
#'
#' Adjust maxent raw outputs to take into account species prevalence
#'
#' @param training_raw_maxent A vector of maxent raw predictions from training data
#' @param testing_raw_maxent A vector of maxent raw predictions from testing data
#' @param estim_prev A numeric value estimating the prevalence of the species
#' @param outType A integer number between 0-3 specifying which output is to be returned. Default is 0 ('raw', 'logistic' and 'cloglog' not adjusted).
#'                A number between 1-3 involves post-processing: 1 = raw adjusted, 2 = cloglog adjusted, 3 = logistic adjusted.
#' @param num_bg A numeric integer specifying the number of background points used for training
#' @param H A numeric value specifying the maximum entropy value reached by the model.Only required when \code{outType=1}.
#' @param maxIter A numeric integer specifying the number of iterations to run for searching parameter `c` or `tau`. Only require when \code{outType=2} or \code{outType=3}.
#' @export
maxent_post_process <- function(training_raw_maxent, testing_raw_maxent, estim_prev, outType=1, num_bg, H, maxIter=10){

  if(missing(training_raw_maxent))
    stop("Raw maxent output from training data is missing.")

  if(!is.numeric(outType))
    stop("Argument 'outType' must be numeric.")
  if(!outType %>% dplyr::between(1,3))
    stop("Argument 'outType' must be between 1 and 3.")

  if(missing(estim_prev))
    stop("Argument 'estim_prevalence' should be provided when post-processing is required.")

  if(!is.numeric(estim_prev) || is.na(estim_prevalence))
    stop("estim_prevalence must be numeric.")

  if(!estim_prev %>% dplyr::between(0,1))
    stop("estim_prevalence must be between 0 and 1.")


  pred_maxent <- switch(outType,

             {
               # raw adjusted
               if(missing(num_bg))
                 stop("num_bg should be provided for the raw post-processing type")

               new_pred  <- training_raw_maxent * estim_prev * num_bg
               new_pred[new_pred>1] <- 1
               new_pred
             },

             {
               # cloglog adjusted (find 'c' so that prevalence over predictions equals estim_prevalence )
               if(missing(testing_raw_maxent))
                 stop("testing_raw_maxent should be provided for the cloglog post-processing type")

               eta <- log(training_raw_maxent) # linear predictor from Maxent (training data)

               # Finding c by direct search to match prevalence (method proposed by Steven Phillips)
               cL <- ifelse( psi.estim(c = -5, eta = eta) > estim_prev, -20, -5) # select appropriate starting points for search
               cH <- ifelse( psi.estim(c = 5, eta = eta) < estim_prev, 20, 5)
               c.estim <- search.c(cL=cL, cH=cH, prev=estim_prev, niter=maxIter, eta=eta)  # find best c

               etaT <- log(testing_raw_maxent) # linear predictor from Maxent (test data)

               1-exp(-exp(etaT+c.estim))
             },

             {
               if(missing(testing_raw_maxent))
                 stop("testing_raw_maxent should be provided for the logistic post-processing type")

               if(missing(H))
                 stop('When logistic post-processing type is required, H must not be missing.')

               psi.logi <- rmaxent::to_logistic(training_raw_maxent, from='raw', H = H) # get default logistic output for training data (tau=0.5)

               tau.estim <- search.tau(tauL=0, tauH=1, prev=estim_prev, niter=maxIter, psi.logi=psi.logi) # find best tau

               psi.logi.T <- rmaxent::to_logistic(testing_raw_maxent, from='raw', H = H) # get default logistic output for test data (tau=0.5)

               1/(1+((1-tau.estim)/tau.estim)*((1-psi.logi.T)/psi.logi.T))  # map to chosen tau

             })

  return(pred_maxent)
}

#'                                                     Project a fitted geographic distance-based model
#'
#' This function is a wrapper for the function `predict` from the `raster` package
#'
#' @param domain An object of class `RasterLayer` representing the spatial grid domain where to interpolate the model
#' @param model An object of class `GeographicDistance` or `InvDistWeightModel`
#' @param output_directory An optional character string specifying the directory where to save the raster of interpolated values
#' @param save_output A logical. Should the output raster be saved ? Ignored if `output_directory` is not provided. Default is `TRUE`.
#' @param return An object of class `RasterLayer` with values interpolated by the model
#' @export
project_geoModel <- function(domain, model, output_directory, output_name=NULL, save_output=TRUE){

  if(missing(domain))
    stop("Domain is missing.")

  if(missing(model))
    stop("Model is missing.")

  if(save_output & missing(output_directory))
    stop("output directory is missing.")

  if(!inherits(model,c("GeographicDistance","InvDistWeightModel")))
    stop("Model must be of class `GeographicDistance` or `InvDistWeightModel`")

  out <- raster::predict(domain, model, mask=TRUE, xyOnly=TRUE)

  if(save_output){
    if(!dir.exists(output_directory))
      stop("Unable to find directory :", output_directory)
    if(is.null(output_name)){
      warning("output_name is null. Will assign absolute data-time as output name.")
      output_name=gsub(" ","_",Sys.time())
    }
    raster::writeRaster(out, file.path(output_directory, paste0(output_name,".tif")))
  }

  return(out)
}



