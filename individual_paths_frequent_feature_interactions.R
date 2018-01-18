individual.paths.frequent.feature.interactions <- function(irf.output, 
                                                           data, 
                                                           sample.indexes, 
                                                           method = 'apriori',
                                                           verbose = F){
  # Generate the frequent used feature interactions for an individual sample
  # Args:
  #   irf.output     : list, the output of iRF2.0 package
  #   data           : matrix, whole sample matrix
  #   sample.indexes : array of integers, a set of samples to look at
  #   method         : string, either 'RIT' or 'apriori'
  #   verbose        : boolean, either to make plots of the final outputs
  # Return:
  #   interactions   : a list of frequent interactions for 
  require('iRF')
  require('arules')
  require('arulesViz')
  # read last forest
  ind <- length(irf.output$rf.list)
  rf.last <- readForest(irf.output$rf.list[[ind]], x=data, return.node.obs=TRUE)
  # get the ids of the target sample
  if (sum(sample.indexes) == 1) {
    # when we only need one sample, ids is easy to obtian.
    ids <- rf.last$node.obs[,sample.indexes]
  } else {
    ids <- apply(rf.last$node.obs[,sample.indexes], FUN = any, MARGIN = 1)
  }
  
  # get the paths for the specific sample
  paths <- subsetReadForest(rf.last, ids)$node.feature
  # find frequent interactions
  if (method == 'apriori') {
    apriori.out <- apriori(t(paths), parameter = list(supp = .1, conf = .4, target = "rules"))
    if (verbose) {
      inspect(apriori.out)
    }
    interaction <- sapply(LIST(items(apriori.out)), function(x) paste(x, collapse = ' '))
    prevalence <- quality(apriori.out)$support
    out <- unique(data.frame(Interaction = interaction, Prevalence = prevalence))
    
  } else if (method == 'RIT') {
    # to be implemented
    out <- RIT(paths)
    #out <- NULL
  } else {
    sprintf('method parameter is not available.')
    out <- NULL
  }
  
  return(out)
}
