.onLoad <- function(...) {
  usefulplants <- data.table::fread(system.file('extdata','utilised_plants_species_list.csv',package="UsefulPlants"), key="binomial_acc_name")
  # load the datasets into the package environment
  assign("usefulplants", usefulplants, envir=parent.env(environment()))
}
