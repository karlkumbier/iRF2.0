# analyze the individual RIT
#   the range of each interaction, which has the largest variance?
load('./individual_enhancer.RData')
interaction.freq = list()
for (out in out.list) {
  out$Interaction <- as.character(out$Interaction)
  for (row.index in 1:nrow(out)) {
    comb <- out$Interaction[row.index]
    if (comb %in% names(interaction.freq)) {
      interaction.freq[[comb]] = c(interaction.freq[[comb]], out$Prevalence[row.index])
    } else {
      interaction.freq[[comb]] = c(out$Prevalence[row.index])
    }
  }
}
tmp <- lapply(interaction.freq, function(x) max(x) - mean(x))
max(do.call('rbind', tmp), na.rm = T)
head(which(tmp > .23))
hist(interaction.freq[['46 54 64']])
varnames.all$Predictor_collapsed[c(46, 64)]
