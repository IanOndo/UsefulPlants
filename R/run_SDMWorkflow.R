#' Running SDM Workflow
#'
#' Runs Useful Plants Modelling workflow
#'
#' @param use_category A character string specifying which category of plant uses the modelling should be performed on
#' @param use_default_settings A logical. Should the models be run with default settings. Default is TRUE.
#' @return None
#' @export
#' @examples
#' run_SDMWorkflow()
run_SDMWorkflow <-function(use_category="All", use_default_settings=TRUE){
  UsefulPlantsWorkflow <- system.file("exdata/", "1UsefulPlants_Workflow.r",package="UsefulPlants")
  if(UsefulPlantsWorkflow == ""){
    stop("Couldn't find the modelling workflow script. Try re-installing UsefulPlants package", call.=FALSE)
  }
  commandArgs <- function(...) paste0(c("use_category","run.with.default.settings"),"=",c(use_category,run.with.default.settings))
  source(system.file("UsefulPlants", "1UsefulPlants_Workflow.r",package="UsefulPlants"))
}
