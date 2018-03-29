load('results_enhancer.Rdata')
load('enhancer.Rdata')
ls()
source('individual_paths_frequent_feature_interactions.R')
out0 <- individual.paths.frequent.feature.interactions(fit, X[train.id,], Y[train.id]==0, method = 'apriori')
out1 <- individual.paths.frequent.feature.interactions(fit, X[train.id,], Y[train.id]==1, method = 'apriori')
index2names <- function(strs, name.list) {
  strs %>% strsplit(' ') %>%
    lapply(FUN = as.numeric) %>% 
    lapply(FUN = function(x) name.list[x]) %>% 
    sapply(FUN = function(x) paste(x, collapse = ' '))
}
out0 %>% mutate(Interaction = index2names(Interaction, varnames.all[[2]]))
require(dplyr)
tmp <- out0 %>% left_join(out1, by = 'Interaction') 
plot(tmp$Prevalence.x, tmp$Prevalence.y)

rf <- fit$rf.list[[3]]
records <- readForest(rf, x = X[train.id,], return.node.obs = TRUE)

#####################################
# compute individual RIT for each   #
#####################################
out.list = list()
for (i in 1:length(train.id)) {
  out.list[[i]] <- individual.paths.frequent.feature.interactions(fit, X[train.id,], train.id==train.id[i], method = 'apriori')
#  if (i > 3) {
#    break
#  }
}
save.image('individual_enhancer.RData')
