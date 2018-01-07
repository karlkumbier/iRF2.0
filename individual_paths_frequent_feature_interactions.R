individual_paths_frequent_feature_interactions <- function(irf_output, sample_index, verbose = T){
  # Generate the frequent used feature interactions for an individual sample
  # Args:
  #   irf_output: the output of irf package
  #   sample_index: integer, which sample to look at
  # Return:
  #   interactions: a list of frequent interactions for the individual
  require('iRF')
  require('arules')
  require('arulesViz')
  # read last forest
  rf_last <- readForest(tail(irf_output$rf.list, n = 1)[[1]], x=x, return.node.obs=TRUE)
  # get the id of the target sample
  id <- rf_last$node.obs[,sample_index]
  # get the paths for the specific sample
  paths <- subsetReadForest(rf_last, id)$node.feature
  # find frequent interactions
  out = apriori(t(paths), parameter = list(supp = .1, conf = .4, target = "rules"))
  if (verbose) {
    inspect(out)
  }
  return(as.character(inspect(items(out))$items))
}
