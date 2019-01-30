stabilityScore <- function(fit, iter, bs.sample, ints.eval, x, y, weights, 
                           varnames.grp, rit.param, signed, n.core, ...) {
  # Wrapper function for stabilityScore_. Calcuates stability of importance
  # metrics across bootstrap samples.

  interact <- list()
  importance <- list()
  for (i in 1:length(bs.sample)) {
    out <- stabilityScore_(sample.id=bs.sample[i], fit=fit, iter=iter,
                           ints.eval=ints.eval, x=x, y=y, 
                           weights=weights, rit.param=rit.param, 
                           varnames.grp=varnames.grp, signed=signed, 
                           n.core=n.core, ...)
    interact[[i]] <- out$interact
    importance[[i]] <- out$importance
  }
  stab <- summarizeInteract(interact.list, ints.full$int)
  imp <- summarizeImp(imp.list)
  return(list(stab=stab, imp=imp))
}

stabilityScore_ <- function(fit, iter, sample.id, ints.eval, x, y, weights,
                            varnames.grp, rit.param, signed, n.core, ...) {    
  # Evalutes interaction importance metrics across a single bootstrap replicate.
  
  if (iter == 1)
    mtry.select.prob <- rep(1, ncol(x))
  else
    mtry.select.prob <- rf.list[[iter - 1]]$importance

  # Fit random forest on bootstrap sample
  rf <- parRF(x=x[sample.id,], y=y[sample.id], xtest=xtest, ytest=ytest,
              mtry.select.prob=mtry.select.prob, ntree=ntree)  

  # Run generalized RIT on rf.b to learn interactions
  ints <- gRIT(rand.forest=rf, x=x, y=y,
               weights=weights,
               varnames.grp=varnames.grp,
               rit.param=rit.param,
               signed=signed,
               ints.full=ints.eval,
               n.core=n.core)

  return(interact=ints$int, imp=ints$imp)
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
    summarize(sta.diff=mean(diff > 0),
              diff=mean(diff),
              sta.prev=mean(prev.test
                            > 0),
              prev1=mean(prev1),
              prev0=mean(prev0),
              sta.prec=mean(prec.test
                            > 0),
              prec=mean(prec)) %>%
    arrange(desc(diff))
  } else {
    imp <- data.table(sta.diff=numeric(0),
                      diff=numeric(0),
                      sta.prev=numeric(0),
                      prev1=numeric(0),
                      prev0=numeric(0),
                      sta.prec=numeric(0),
                      prec=numeric(0))
  }

  return(data.table(imp))
}

