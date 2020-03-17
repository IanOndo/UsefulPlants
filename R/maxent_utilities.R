#' From a list to a command-line for MaxEnt
#'
#' Reads a named list and returns an executing command line for MaxEnt
#'
#' @param argList A list object whose elements refer to Maxent arguments
#'
#' @export
List2MaxEntCommand <- function(argList){

  if(missing(argList))
    stop("Argument argList is missing")

  if(sum(lengths(argList))==0L)
    stop("Argument argList is of length 0L")

  if(all(sapply(names(argList),nchar,USE.NAMES=F)==0L))
    stop("All elements in argList must be named")

  cmd = ''
  arg_names  = names(argList)

  for(arg in arg_names){

    flag = NULL

    flag = switch(arg,

                  'responsecurves' 	=	if(argList[[arg]])'-P' else '',
                  'pictures' 			=	paste0(arg,'=',argList[[arg]]),
                  'jackknife' 		=	if(argList[[arg]]) '-J' else '',
                  'outputformat' 		=	paste0(arg,'=',argList[[arg]]),
                  'outputfiletype' 	=	paste0(arg,'=',argList[[arg]]),
                  #outputdirectory' 	=	paste('-o',argList[[arg]]),
                  'projectionlayers' 	=	paste('-j',argList[[arg]]),


                  #samplesfile' 		=	paste('-s',argList[[arg]])
                  #environmentallayers'=	paste('-e',argList[[arg]])

                  'randomseed'	 	=	paste0(arg,'=',argList[[arg]]),

                  'logscale' 			=	paste0(arg,'=',argList[[arg]]),
                  'warnings' 			=	paste0(arg,'=',argList[[arg]]),

                  'tooltips' 			=	paste0(arg,'=',argList[[arg]]),
                  #askoverwrite'		= 	if(argList[[arg]]) '-r' else '',

                  'skipifexists'		= 	if(argList[[arg]]) '-S' else '',

                  'removeduplicates'	=	paste0(arg,'=',argList[[arg]]),


                  'writeclampgrid'	= 	paste0(arg,'=',argList[[arg]]),

                  'writemess'			= 	paste0(arg,'=',argList[[arg]]),

                  'randomtestpoints'	=	paste('-X',argList[[arg]]),
                  'betamultiplier' 	=	paste('-b',argList[[arg]]),
                  'maximumbackground'	=	paste('-B',argList[[arg]]),
                  'biasfile'			=	paste0(arg,'=',argList[[arg]]),


                  'testsamplesfile'	= 	paste('-T',argList[[arg]]),


                  'replicates'		=	paste0(arg,'=',argList[[arg]]),
                  'replicatetype'		=	paste0(arg,'=',argList[[arg]]),



                  'perspeciesresults'	=	paste0(arg,'=',argList[[arg]]),
                  'writebackgroundpredictions'=	paste0(arg,'=',argList[[arg]]),
                  'responsecurvesexponent'=	paste0(arg,'=',argList[[arg]]),
                  'linear'			=	if(!argList[[arg]])'-l' else '',
                  'quadratic'			= 	if(!argList[[arg]])'-q' else '',
                  'product' 			=	if(!argList[[arg]])'-p' else '',
                  'threshold'			=	if(!argList[[arg]]) 'nothreshold' else '',
                  'hinge'				= 	if(!argList[[arg]])'-h' else '',
                  'addsamplestobackground'=	if(!argList[[arg]])'-d' else '',
                  'addallsamplestobackground'= paste0(arg,'=',argList[[arg]]),
                  'autorun'			= 	if(argList[[arg]])'-a' else '',
                  'writeplotdata'		=	paste0(arg,'=',argList[[arg]]),
                  'fadebyclamping'	= 	paste0(arg,'=',argList[[arg]]),

                  'extrapolate'		=	paste0(arg,'=',argList[[arg]]),
                  'visible'			=	if(!argList[[arg]]) '-z' else '',
                  'autofeature'		= 	if(!argList[[arg]])'-A' else '',
                  'doclamp'			= 	paste0(arg,'=',argList[[arg]]),
                  'outputgrids'		= 	if(!argList[[arg]]) '-x' else '',
                  'plots'				= 	if(!argList[[arg]]) 'noplots' else '',
                  'appendtoresultsfile'= 	paste0(arg,'=',argList[[arg]]),
                  'maximumiterations'=	paste('-m',argList[[arg]]),
                  'convergencethreshold'= paste('-c',argList[[arg]]),
                  'adjustsampleradius'= 	paste0(arg,'=',argList[[arg]]),

                  'threads'			= paste0(arg,'=',argList[[arg]]),
                  'lq2lqptthreshold'	= paste0(arg,'=',argList[[arg]]),
                  'l2lqthreshold'		= paste0(arg,'=',argList[[arg]]),
                  'hingethreshold'	= paste0(arg,'=',argList[[arg]]),
                  'beta_threshold'	= paste0(arg,'=',argList[[arg]]),
                  'beta_categorical'	= paste0(arg,'=',argList[[arg]]),
                  'beta_lqp'			= paste0(arg,'=',argList[[arg]]),
                  'beta_hinge'		= paste0(arg,'=',argList[[arg]]),
                  'logfile'			= paste0(arg,'=',argList[[arg]]),
                  'cache'				= paste0(arg,'=',argList[[arg]]),
                  'defaultprevalence'	= paste0(arg,'=',argList[[arg]]),

                  'applythresholdrule'= shQuote(paste0(arg,'=',argList[[arg]]),type='sh'),
                  'togglelayertype'	=	paste('-t',argList[[arg]]),
                  'togglespeciesselected'= paste('-E',argList[[arg]]),
                  'togglelayerselected'= 	paste('-N',argList[[arg]]),
                  'verbose'			= 	if(argList[[arg]])	'-v' else '',
                  'allowpartialdata' = paste0(arg,'=',argList[[arg]]),
                  'prefixes'			= paste0(arg,'=',argList[[arg]]),
                  'nodata'			=	paste('-n',argList[[arg]])
    )
    if(is.null(flag))
      stop(paste('Argument',arg,'does not exist'))

    if(arg=='visible'){
      if(!argList[[arg]])
        flag <- paste(flag,'-a')
    }
    cmd = paste(cmd,flag)
  }

  return(cmd)
}

#' From an output directory
#'
#' Reads a named list and returns an executing command line for MaxEnt
#'
#' @param  maxent_outputdir list object whose elements refer to Maxent arguments
#' @export
get_maxent_output <- function(maxent_outputdir, index=1){

  # get evaluation data
  samplePredictions     <- read.csv(list.files(maxent_outputdir,pattern=".+samplePredictions.csv$",full.names=TRUE),h=TRUE)
  backgroundPredictions <- read.csv(list.files(maxent_outputdir,pattern=".+backgroundPredictions.csv$",full.names=TRUE),h=TRUE)
  maxentResults         <- read.csv(file.path(maxent_outputdir,"maxentResults.csv"),h=TRUE)
  lbds <- rmaxent::parse_lambdas(list.files(maxent_outputdir,pattern="\\.lambdas",full.names=TRUE))

  training_raw_maxent <- samplePredictions %>%
    dplyr::filter(Test.or.train=='train') %>%
    dplyr::pull(Raw.prediction)

  testing_raw_maxent <- samplePredictions %>%
    dplyr::filter(Test.or.train=='test') %>%
    dplyr::pull(Raw.prediction)

  bg_raw_maxent <- backgroundPredictions %>%
    dplyr::pull(raw)

  n_trainSamples <- length(training_raw_maxent)
  n_testSamples <- length(testing_raw_maxent)
  n_samplesBG <- lbds$numBackgroundPoints

  trained_pred_at_occ <- rmaxent::to_cloglog(training_raw_maxent, from='raw', H=lbds$entropy) # maxent_post_process(training_raw_maxent, testing_raw_maxent, estim_prev = 0.3, outType = 1, num_bg = n_samplesBG, H = lbds$entropy)
  test_pred_at_occ    <- rmaxent::to_cloglog(testing_raw_maxent, from='raw', H=lbds$entropy)  # maxent_post_process(testing_raw_maxent, testing_raw_maxent, estim_prev = 0.3, outType = 1, num_bg = n_samplesBG, H = lbds$entropy)
  pred_at_bg          <- rmaxent::to_cloglog(bg_raw_maxent, from='raw', H=lbds$entropy)       # maxent_post_process(bg_raw_maxent, testing_raw_maxent, estim_prev = 0.3, outType = 1, num_bg = n_samplesBG, H = lbds$entropy)

  # evaluating data
  train_predicted <- c(trained_pred_at_occ, pred_at_bg)
  train_labels    <- c(rep(1, n_trainSamples),rep(0, length(pred_at_bg)))

  test_predicted <- c(test_pred_at_occ, pred_at_bg)
  test_labels  <- c(rep(1, n_testSamples),rep(0, length(pred_at_bg)))

  pred_test <- ROCR::prediction(test_predicted, labels = test_labels)

  # evaluation metrics
  auc_train <- ROCR::performance(ROCR::prediction(train_predicted, labels = train_labels), measure = "auc")@y.values[[1]]
  auc_test  <- ROCR::performance(pred_test, measure = "auc")@y.values[[1]]
  tss_test  <- max(ROCR::performance(pred_test, measure = "sens")@y.values[[1]] + ROCR::performance(pred_test, measure = "spec")@y.values[[1]] - 1)
  or10_test <- maxentResults[['X10.percentile.training.presence.test.omission']]

  # evaluation output
  data.frame(Species=species_name,
             Bin=index,
             n_train_samples=n_trainSamples, n_test_samples=n_testSamples,
             TSS_test=tss_test, AUC_train=auc_train, AUC_test=auc_test, AUC_diff=auc_test-auc_train, OR10_test=or10_test)
}


# utility function: get command line to repeat from maxent output
#' @export
get_maxent_command_line <- function(fn){
  pattern <-"^Command line to repeat\\: "
  line <- grep(pattern, base::readLines(fn), value=TRUE)
  command_line <- gsub(pattern, "", line)
  return(command_line)
}

# utility function: modify maxent command line with new argument
#' @export
set_command_param <- function(cmd, param_name, new_value){

  pattern = paste0(param_name,"=.[^ ]+")
  replacement = paste0(param_name,"=",new_value)
  if(grepl(pattern, cmd))
    new_command_line <- stringr::str_replace(cmd, pattern, replacement)
  else
    new_command_line <- paste(cmd, replacement)

  return(new_command_line)
}

# utility function: modify maxent command line with new argument
#' @export
set_command_arg <- function(cmd, arg_name, new_arg){

  new_command_line <- stringr::str_replace(cmd, arg_name, new_arg)

  return(new_command_line)
}

#' function that calculates the mean of the cloglog output, based on linear predictor 'eta' and constant 'c'
#' @export
psi.estim <- function(c, eta) {mean(1-exp(-exp(eta+c)))}

#' function to search for the value of 'c' in cloglog that approximates the prevalence
#' @export
search.c <- function (cL, cH, prev, niter, eta){
  for (ii in 1:niter){
    cM <- cL+(cH-cL)/2   # find mid-value
    if( psi.estim(cM,eta) > prev){cH <- cM}
    else {cL <- cM}
  }
  if(abs((psi.estim((cL+cH)/2,eta))-prev)>0.1) {
    print("WARNING: 'c' obtained with not enough resolution, increase n.iter")
  }
  return((cL+cH)/2)
}

#' function that implements the logistic output, based on linear predictor 'eta[i]', the relative entropy 'r' and parameter 'tau'
#' @export
psi.estim.logistic <- function(tau, psi.logi) {mean(1/(1+((1-tau)/tau)*((1-psi.logi)/psi.logi)))}

#' function to search for the value of 'tau' that approximates the observed prevalence given in prev
#' @export
search.tau <- function (tauL, tauH, prev, niter, psi.logi){
  for (ii in 1:niter){
    tauM <- tauL+(tauH-tauL)/2   # find mid-value
    if( psi.estim.logistic(tauM,psi.logi) > prev){tauH <- tauM}
    else {tauL <- tauM}
  }
  if(abs((psi.estim.logistic((tauL+tauH)/2,psi.logi))-prev)>0.1) {
    print("WARNING: 'tau' obtained with not enough resolution, increase n.iter")
  }
  return((tauL+tauH)/2)
}
