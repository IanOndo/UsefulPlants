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
maxent <- function(sp_points, env_layers, outputdir, path_to_maxent=NULL, path_to_java=NULL, memory_allocated=512, wait=TRUE, ...){

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
  main_cmd = paste(java_interpreter,paste0('-mx',memory_allocated,'m'),'-jar', paste0('"',maxent_jar,'"'),'-a -r','-e',env_layers,'-s', sp_points,'-o',outputdir)

  # add extra arguments
  extra_cmd = NULL
  argList = list(...)
  if(length(argList)>0L)
    extra_cmd = gsub("'",'"',List2MaxEntCommand(argList))

  maxent_cmd = paste(main_cmd,extra_cmd)

  exec <- system(maxent_cmd,wait=wait)
}
