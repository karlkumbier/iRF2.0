stabilityScore <- function(fit, x, y, iter, bs.sample, ints.eval, weights,
                           varnames.grp, rit.param, signed, n.core, ntree,
                           ...) {
  # Wrapper function for bsgRIT. Calcuates stability of importance
  # metrics across bootstrap samples.
  out <- list()
  for (i in 1:length(bs.sample)) {

    sample.id <- bs.sample[[i]]
    out[[i]] <- bsgRIT(fit, x, y, iter, sample.id, ints.eval=ints.eval, 
                       ntree=ntree, weights=weights, rit.param=rit.param,
                       varnames.grp=varnames.grp, signed=signed, n.core=n.core,
                       ...)

  }

  # Summarize stability and importance metrics across bootstrap replicates
  out <- summarizeInteract(out)
  return(out)
}


bsgRIT <- function(fit, x, y, iter, sample.id, ints.eval, weights, ntree,
                   varnames.grp, rit.param, signed, n.core, ...) {
  # Fit RF on single bootstrap replicate and evalutes interaction importance 
  # metrics for fitted RF. 

  if (iter == 1)
    mtry.select.prob <- rep(1, ncol(x))
  else
    mtry.select.prob <- fit[[iter - 1]]$importance

  # Fit random forest on bootstrap sample
  rf <- parRF(x=x[sample.id,], y=y[sample.id], xtest=xtest, ytest=ytest,
              mtry.select.prob=mtry.select.prob, ntree=ntree,
              n.core=n.core,
              ...)

  # Run generalized RIT on rf.b to learn interactions
  ints <- gRIT(rand.forest=rf, x=x, y=y,
               weights=weights,
               varnames.grp=varnames.grp,
               rit.param=rit.param,
               signed=signed,
               ints.full=ints.eval,
               n.core=n.core)

  return(ints)
}

summarizeInteract <- function(x) {
  # Summarize interaction importance metrics across bootstrap samples
  n.bootstrap <- length(x)
  x <- rbindlist(x)

  if (nrow(x) > 0) {
    imp <- mutate(x, diff=(prev1-prev0)) %>%
      group_by(int) %>%
      summarize(prevalence.diff=mean(diff),
                sta.diff=mean(diff > 0),
                independence=mean(prev.test),
                sta.independence=mean(prev.test > 0),
                precision=mean(prec),
                sta.precision=mean(prec.test > 0),
                stability=mean(recovered)) %>%
      arrange(desc(prevalence.diff))
  } else {
    imp <- data.table(prevalence.diff=numeric(0),
                      sta.diff=numeric(0),
                      independence=numeric(0),
                      sta.independence=numeric(0),
                      precision=numeric(0),
                      sat.precision=numeric(0),
                      stability=numeric(0))
  }

  return(data.table(imp))
}

