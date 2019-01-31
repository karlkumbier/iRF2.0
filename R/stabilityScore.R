stabilityScore <- function(fit, iter, bs.sample, ints.eval, x, y, weights, 
                           varnames.grp, rit.param, signed, n.core, ntree, 
                           ...) {
  # Wrapper function for stabilityScore_. Calcuates stability of importance
  # metrics across bootstrap samples.

  interact <- list()
  importance <- list()
  for (i in 1:length(bs.sample)) {
    out <- stabilityScore_(sample.id=bs.sample[[i]], fit=fit, iter=iter,
                           ints.eval=ints.eval, x=x, y=y, ntree=ntree,
                           weights=weights, rit.param=rit.param, 
                           varnames.grp=varnames.grp, signed=signed, 
                           n.core=n.core, ...)
    interact[[i]] <- out$stab
    importance[[i]] <- out$imp
  }
  stab <- summarizeInteract(interact, ints.eval)
  imp <- summarizeImp(importance)
  return(list(stab=stab, imp=imp))
}

stabilityScore_ <- function(fit, iter, sample.id, ints.eval, x, y, weights,
                            varnames.grp, rit.param, signed, n.core, ntree,
                            ...) {    
  # Evalutes interaction importance metrics across a single bootstrap replicate.
  
  if (iter == 1)
    mtry.select.prob <- rep(1, ncol(x))
  else
    mtry.select.prob <- fit[[iter - 1]]$importance

  # Fit random forest on bootstrap sample
  rf <- parRF(x=x[sample.id,], y=y[sample.id], xtest=xtest, ytest=ytest,
              mtry.select.prob=mtry.select.prob, ntree=ntree, n.core=n.core)  

  # Run generalized RIT on rf.b to learn interactions
  ints <- gRIT(rand.forest=rf, x=x, y=y,
               weights=weights,
               varnames.grp=varnames.grp,
               rit.param=rit.param,
               signed=signed,
               ints.full=ints.eval,
               n.core=n.core)

  return(list(stab=ints$int, imp=ints$imp))
}

summarizeInteract <- function(x, ints.full){
  # Aggregate interactions across bootstrap samples

  n.bootstrap <- length(x)
  x <- unlist(x)
  out <- rep(0, length(ints.full))
  names(out) <- ints.full

  if (length(x) >= 1){
    int.tbl <- sort(c(table(x)), decreasing=TRUE)
    int.tbl <- int.tbl / n.bootstrap
    out[names(int.tbl)] <- int.tbl
    return(out)
  } else {
    return(c(interaction=numeric(0)))
  }
}

summarizeImp <- function(imp) {
  # Summarize interaction importance metrics across bootstrap samples
  imp <- rbindlist(imp)

  if (nrow(imp) > 0) {
    imp <- mutate(imp, diff=(prev1-prev0)) %>%
    group_by(int) %>%
    summarize(prevalence.diff=mean(diff),
              sta.diff=mean(diff > 0),
              independence=mean(prev.test),
              sta.independence=mean(prev.test > 0),
              precision=mean(prec),
              sta.precision=mean(prec.test > 0)) %>%
    arrange(desc(prevalence.diff))
  } else {
    imp <- data.table(prevalence.diff=numeric(0),
                      sta.diff=numeric(0),
                      independence=numeric(0),
                      sta.independence=numeric(0),
                      precision=numeric(0),
                      sta.precision=numeric(0))
  }

  return(data.table(imp))
}

