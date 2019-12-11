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

  if(all(sapply(names(L),nchar,USE.NAMES=F)==0L))
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
                  'maximumbackground'	=	paste('-MB',argList[[arg]]),
                  'biasfile'			=	paste0(arg,'=',argList[[arg]]),


                  'testsamplesfile'	= 	paste('-T',argList[[arg]]),


                  'replicates'		=	paste0(arg,'=',argList[[arg]]),
                  'replicatetype'		=	paste0(arg,'=',argList[[arg]]),



                  'perspeciesresults'	=	paste0(arg,'=',argList[[arg]]),
                  'writebackgroundpredictions'=	paste0(arg,'=',argList[[arg]]),
                  'responsecurvesexponent'=	paste0(arg,'=',argList[[arg]]),
                  'linear'			=	if(argList[[arg]])'-l' else '',
                  'quadratic'			= 	if(argList[[arg]])'-q' else '',
                  'product' 			=	if(argList[[arg]])'-p' else '',
                  'threshold'			=	paste0(arg,'=',argList[[arg]]),
                  'hinge'				= 	if(argList[[arg]])'-h' else '',
                  'addsamplestobackground'=	if(argList[[arg]])'-d' else '',
                  'addallsamplestobackground'= paste0(arg,'=',argList[[arg]]),
                  #autorun'			= 	if(argList[[arg]])'-a' else '',
                  'writeplotdata'		=	paste0(arg,'=',argList[[arg]]),
                  'fadebyclamping'	= 	paste0(arg,'=',argList[[arg]]),

                  'extrapolate'		=	paste0(arg,'=',argList[[arg]]),
                  #visible'			=	if(argList[[arg]])'-z'
                  'autofeature'		= 	if(argList[[arg]])'-A' else '',
                  'doclamp'			= 	paste0(arg,'=',argList[[arg]]),
                  'outputgrids'		= 	paste('-x',argList[[arg]]),
                  'plots'				= 	paste0(arg,'=',argList[[arg]]),
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
                  'verbose'			= 	if(argList[[arg]])	'-v',
                  'allowpartialdata' = paste0(arg,'=',argList[[arg]]),
                  'prefixes'			= paste0(arg,'=',argList[[arg]]),
                  'nodata'			=	paste('-n',argList[[arg]])
    )
    if(is.null(flag))
      stop(paste('Argument',arg,'does not exist'))

    cmd = paste(cmd,flag)
  }

  return(cmd)
}

#' Performing Species Distribution Modelling with MaxEnt
#'
#' Runs MaxEnt on a species occurrence dataset and with given environmental variables
#'
#' @param sp_points A character string specifying the path to the species occurrence csv files
#' @param env_layers A character string specifying the path to the directory containing environmental predictors layers
#' @param outputdir A character string specifying the path to the output directory for maxent results.
#' @param path_to_maxent A character string specifying the path to the maxent .jar file
#' @param path_to_java A character string specifying the path to a java interpreter
#' @param wait A logical. Should the function wait until the job is finished ? Default is TRUE.
#' @param ... Additional parameters to be passed to maxent.
#'
#' @export
maXent <- function(sp_points, env_layers, outputdir, path_to_maxent=NULL, path_to_java=NULL, wait=TRUE, ...){

  # check if a java interpreter is available
  java_interpreter = 'java'
  if(!file.exists(Sys.which(java_interpreter))){
    if(is.null(path_to_java)){
      cat('JAVA interpreter does not exist in your PATH environment variables')
      cat('You must provide a path to your java interpreter')
      return()
    }
    else{
      if(!file.exists(path_to_java))
        stop(paste('Unable to find: ',path_to_java,'. This file does not exist or is not accessible'))
      java_interpreter = path_to_java
    }
  }
  # check if the maxent file .jar exists
  maxent_jar <- paste(system.file(package="dismo"),"/java/maxent.jar", sep='')
  if(!file.exists(maxent_jar)){
    if(is.null(path_to_maxent)){
      stop('- MAXENT jar file has not been found on your computer -\n - You must provide a path to your maxent program')
    }
    else{
      if(!file.exists(path_to_maxent))
        stop(paste('File:',path_to_maxent,'does not exist or is not accessible'))
      maxent_jar = path_to_maxent
    }
  }
  main_cmd = paste(java_interpreter,'-jar', paste0('"',maxent_jar,'"'),'-a -r','-e',env_layers,'-s', sp_points,'-o',outputdir)

  # add extra arguments
  extra_cmd = NULL
  argList = list(...)
  if(length(argList)>0L)
    extra_cmd = gsub("'",'"',List2MaxEntCommand(argList))

  maxent_cmd = paste(main_cmd,extra_cmd)

  exec <- system(maxent_cmd,wait=wait)
}

