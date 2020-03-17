#' Running SDM Workflow
#'
#' Runs Useful Plants Modelling workflow
#'
#' @param base_directory A character string specifying the directory where to store the outputs of the workflow
#' @param species_directory A character string specifying the directory where are stored the species occurrence records.
#' @param env_directory A character string specifying the directory where are stored the environmental layers.
#' @param which_uses A character string specifying the name of the plant use category to model among: `AnimalFood`, `HumanFood`,`InvertebratesFood`,`Medicinals`
#' `EnvironmentalUses`,`Materials`,`SocialUses`,`GeneSources`,`Fuels`,`Poisons`. Default is `AllUses`, which means that outputs of all categories will be aggregated.
#' @return None
#' @export
#' @examples
#' run_SDMWorkflow()
run_SDMWorkflow <-function(base_directory,
                           species_directory,
                           env_directory,
                           which_uses="AllUses",
                           use_default_settings=TRUE){

  UsefulPlantsWorkflow <- system.file("exdata/", "1UsefulPlants_Workflow.R",package="UsefulPlants")
  if(UsefulPlantsWorkflow == ""){
    stop("Couldn't find the modelling workflow script. Try re-installing UsefulPlants package", call.=FALSE)
  }

  #=================
  #= 0.Load packages
  #=================
  pkgs.to.load  <- c('dplyr')
  pkgs.loaded   <- sapply(pkgs.to.load,require,character.only=TRUE)
  if(!any(pkgs.loaded)){
    warning(paste('Packages',paste(pkgs.to.load[!pkgs.loaded],collapse=', '),'failed to load'))
    stop("Try to re-install packages that failed to load")
  }

  #=================
  #= 1. Check inputs
  #=================
  if(missing(base_directory)){
    warning("Base directory is missing, and will be set to the current working directory");flush.console()
    base_directory = getwd()
  }
  if(missing(species_directory))
    stop("Species directory is missing.")

  if(missing(env_directory))
    env_directory=NULL

  if(!is.null(mc_cores) & !is.numeric(mc_cores))
    stop("Argument 'mc_cores' must be a numeric integer")

  if(is.null(mc_cores)) mc_cores = parallel::detectCores()-1 else min(mc_cores, parallel::detectCores())

  #===================
  #= 2. Set parameters
  #===================
  params = c(
    'base_dir' 			= base_directory,

    'species_dir' 	= species_directory,		# path to species csv files directory

    'env_dir' 			= env_directory,			# path to environmental predictors directory

    'which_uses'    = which_uses,

    'points_proj' 	= points_proj,		# character string of occurrence points projection

    'run_name'			= which_uses,

    'mc_cores' 			= round(mc_cores)  # number of cores to use
  )

  #==================
  #== 3. Run workflow
  #==================
  commandArgs <- function(...) paste0(names(params),'=',params)
  source(system.file("UsefulPlants", "1UsefulPlants_Workflow.R",package="UsefulPlants"))
}
