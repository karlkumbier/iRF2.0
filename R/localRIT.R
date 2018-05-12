localORIT <- function(idcs, rf, cl=1) {
  # Run observation interacion RIT over
  # class-cl leaf nodes containing 
  # observations indicated by idcs
 
  id.cl <- rf$tree.info$prediction == cl + 1
  id.nd <- apply(as.matrix(rf$node.obs[,idcs]) != 0, MAR=1, all)
  os1 <- rf$node.obs[id.cl & id.nd,]
  ti1 <- rf$tree.info[id.cl & id.nd,]
  rit <- RIT(os1, weights=ti1$size.node)$Interaction
  rit <- unname(gsub(' ', '_', rit))
  prev1 <- sapply(rit, prevalence, nf=os1, wt=ti1$size.node)

  if (sum(!id.cl & id.nd) > 0) {
    os0 <- rf$node.obs[!id.cl & id.nd,]
    ti0 <- rf$tree.info[!id.cl & id.nd,]
    prev0 <- sapply(rit, prevalence, nf=os0, wt=ti0$size.node)
  } else {
    prev0 <-  rep(0, length(prev1))
  }

  pd <- prev1 - prev0  
  d <- data.frame(int=rit, p1=prev1, p0=prev0, d=pd)
  rownames(d) <- NULL 
  return(d)
}

localFRIT <-  function(idcs, rf, cl=1, varnames=NULL) {
  # Run feature interaction RIT over 
  # class-cl leaf nodes containing 
  # observations indicated by idcs
  id.cl <- rf$tree.info$prediction == cl + 1
  id.nd <- apply(as.matrix(rf$node.obs[,idcs]) != 0, MAR=1, all)
  fs1 <- rf$node.feat[id.cl & id.nd,]
  ti1 <- rf$tree.info[id.cl & id.nd,]
  rit <- RIT(fs1, weights=ti1$size.node)$Interaction
  rit <- unname(gsub(' ', '_', rit))
  prev1 <- sapply(rit, prevalence, nf=fs1, wt=ti1$size.node)

  if (sum(!id.cl & id.nd) > 0) {
    fs0 <- rf$node.feat[!id.cl & id.nd,]
    ti0 <- rf$tree.info[!id.cl & id.nd,]
    prev0 <- sapply(rit, prevalence, nf=fs0, wt=ti0$size.node)
  } else {
    prev0 <-  rep(0, length(prev1))
  }
  
  if (!is.null(varnames)) rit <- nameInts(rit, varnames)
  pd <- prev1 - prev0
  d <- data.frame(int=rit, p1=prev1, p0=prev0, d=pd)
  rownames(d) <- NULL
  return(d)
}

