nfFilterInt <- function(nf, interact) {
  # Return leaf node indices where interact falls on decision path
  # args:
  #   nf: binary node feature matrix indicating variables on decision paths 
  #     for each leaf node
  #   interact: numeric vector indicating interaction variables
  is.int <- apply(nf[,interact], MAR=1, function(z) all(as.logical(z)))
  return(is.int)
}

filterInt <- function(paths, interact) {
  # Filter paths to include only splits with given interaction variables
  # args:
  #   paths: a data frame containing metadata on each decision path
  #   interact: a numeric vector indicating interaciton variables
  require(dplyr)
  filtered <- filter(paths, vars %in% interact)
  return(filtered)
}

getPathsTree <- function(rf.obj, nodes, tree.idx, interact, 
                         varnames.grp, min.node.size=5) {
  # Extract decision paths across all leaf nodes in a given tree.
  # args:
  #   tree.info: data frame as returned by readForest
  #   tree.idx: index of the tree to return paths for
  #   rf.obj: a random forest object as returned by randomForest
  #   interact: a numeric vector indicating interaciton variables
  #   varnames.grp: feature grouping vector
  require(data.table)
  require(dplyr)
  
  tree.info <- as.data.frame(getTree(rf.obj, tree.idx))
  tree.info$node.idx <- 1:nrow(tree.info)
  
  # get node size and purity attributes for current tree
  which.leaf <- data.table(table(nodes))
  names(which.leaf) <- c('node.idx', 'size.node')
  which.leaf$node.idx <- as.numeric(which.leaf$node.idx)
  
  out <- rbindlist(getPaths(tree.info, varnames.grp))
  out <- filterInt(out, interact)
  out$tree <- tree.idx

  out <- select(tree.info, prediction, node.idx) %>% 
    right_join(out, by='node.idx') %>%
    right_join(which.leaf, by='node.idx')
  
  if (min.node.size < 1 & rf.obj$type == 'classification') { 
    # if min.node.size < 1, it is assumed to be given as quantile
    m1 <- quantile(out$size.node[out$prediction == 2], min.node.size)
    m0 <- quantile(out$size.node[out$prediction == 1], min.node.size)
    out <- filter(out, size.node > m1 & prediction == 2 | 
                    size.node > m0 & prediction == 1)
  } else if (min.node.size < 1) {
    m1 <- quantile(out$size.node, min.node.size)
    out <- filter(out, size.node > m)
  } else {
    out <- filter(out, size.node > min.node.size)
  }
  
  return(out)
}

getPaths <- function(tree.info, varnames.grp) {
  # Extract decision path information for select leaf nodes.
  # args:
  #   tree.info: tree info matrix as returned by readForest
  #   varnames.grp: feature grouping vector
  
  paths <- getAncestorPath(tree.info, varnames.grp, split.pt=TRUE)
  p <- nrow(paths) / 2
  out <- lapply(as.character(which(tree.info$status == -1)), function(z) {
    zz <- paths[,z]
    path.splits <- zz[zz != 0]
    path.vars <- which(zz != 0) 
    path.dirs <- ifelse(path.vars > p, 1, 0)  
    path.vars <- path.vars %% p
    path.vars[path.vars == 0] <- p
    out <- data.table(dirs=path.dirs, vars=path.vars, 
                      splits=path.splits, node.idx=as.numeric(z))
    return(out)
  })
  
  return(out)
}

interactHyperrctangles <- function(paths, interact, x) {
  # Extracts hyperrectangle information specific to a given interaction for 
  # each leaf node.
  # args:
  #   paths: a data frame containing metadata on each decision path
  #   interact: numeric vector indicating interaction variables
  #   x: data matrix used to determine range of interacting features
  require(dplyr)
  
  # Filter out interaction variables and determine range of features for each
  paths <- group_by(paths, vars) %>% 
    mutate(vmin=min(x[,vars]), vmax=max(x[,vars])) %>%
    ungroup()
    
  # Generate rectangle objects for each leaf node
  suppressWarnings(
    paths.rect <- arrange(paths, vars) %>%
      group_by(tree, node.idx) %>%
      do(rect=makeRectangle(.), int=interact) %>%
      ungroup()
  )
  
  # Set leaf node properties: prediction, number observations
  suppressMessages(
    paths.rect <- select(paths, tree, node.idx, prediction, size.node) %>%
      right_join(paths.rect, by=c('tree', 'node.idx'))
  )
  paths.rect <- paths.rect[seq(1, nrow(paths.rect), by=2),]
  return(paths.rect)
}

forestHyperrectangle <- function(rf.obj, x, y, interact, varnames.grp=NULL, 
                                 groupFUN=mean, min.node.size=1, n.cores=1) {
  # Read hyperrectangles from RF for a specified interactin
  # args:
  #   rf.obj: randomForest object
  #   read.forest: output of readForest
  #   x: data matrix
  #   y: response vector
  #   interact: vector of features for which to generate hyperrectangles
  #   min.node.size: filter all leaf nodes that are not larger than specified 
  #     value for faster processing
  require(parallel)
  
  interact <- strsplit(interact, '_')[[1]]
  if (!is.null(varnames.grp))
    interact <- which(unique(varnames.grp) %in% interact)
  else if (is.character(interact)) 
    interact <- which(rownames(rf.obj$importance) %in% interact)
  if (length(interact) == 0)
    stop ('interaction features not in rf.obj')
  
  print('Extracting paths...')
  prf <- predict(rf.obj, newdata=x, nodes=TRUE)
  nodes <- attr(prf, 'nodes')
  
  tree.paths <- mclapply(1:rf.obj$ntree, function(k) {
    getPathsTree(rf.obj=rf.obj, nodes=nodes[,k], tree.idx=k,
                 min.node.size=min.node.size,
                 interact=interact, varnames.grp=varnames.grp)
  }, mc.cores=n.cores)
  tree.paths <- rbindlist(tree.paths)
  
  if (is.factor(y)) 
    tree.paths <- mutate(tree.paths, prediction = prediction - 1)
  
  tree.paths <- group_by(tree.paths, tree, node.idx) %>%
    mutate(int.size=length(unique(vars))) %>%
    filter(int.size == length(interact)) %>%
    ungroup()
  
  return(tree.paths)
}

plotInt2 <- function(rectangles, interact.plot, x, y, 
                     grid.size=100, col.pal=c('blue', 'yellow'), 
                     xlab=NULL, ylab=NULL, zlab=NULL,
                     range.color=NULL,
                     z.range=c(0,1),
                     n.cols=100, grids=NULL, axes=TRUE,
                     varnames.grp=colnames(x),
                     groupFUN=mean,
                     pred.prob=FALSE, min.node.size=10) {
  # Generates surface map of order-2 interaction
  # args:
  #   rectangles: hyprerectangle list as generated by forestHyperrectangle
  #   interact.plot: order-2 interaction to plot
  #   grid.size: number of bins to plot surface map over
  #   col.pal: color palette of surface map
  #   xlab, ylab, zlab: axis labels
  #   range.color: range of color values
  #   z.range: range for response axis
  #   n.cols: number of colors in color pal
  #   grids: user generated grid to plot over. If supplied, grid.size is 
  #     ignored
  #   axes: T/F indicating whether axes should be plotted
  #   varnames.group: character vector indicating feature grouping
  #   groupFUN: function to use for grouping columns of x
  #   pred.prob: T/F indicating whether the z axis should indicate predicted
  #     probability from the forest or raw data distribution
  #   min.node.size: minimum leaf node size to use for grid
  require(rgl)

  n <- nrow(x)
  interact.plot <- strsplit(interact.plot, '_')[[1]]
  
  if (!is.null(varnames.grp)) 
    interact.plot <- which(unique(varnames.grp) %in% interact.plot)
  else
    interact.plot <- as.numeric(interact.plot)
  
  stopifnot(length(interact.plot) == 2) 

  # Group data matrix if grouping features indicated
  if (!is.null(varnames.grp)) {
    stopifnot(length(varnames.grp) == ncol(x))
    x <- sapply(unique(varnames.grp), function(g) {
      apply(x[,varnames.grp == g], MAR=1, groupFUN)
    })
  }
  
  if (is.factor(y)) {
    y <- as.numeric(y) - 1
    rectangles$prediction <- rectangles$prediction - 1
  }
  
  # generate grid to plot surface over either as raw values or quantiles
  intp1 <- interact.plot[1]
  intp2 <- interact.plot[2]
  
  if (is.null(grids)) {
    g1 <- seq(min(x[,intp1]), max(x[,intp1]), length.out=grid.size)
    g2 <- seq(min(x[,intp1]), max(x[,intp2]), length.out=grid.size)
    g1.n <- g1
    g2.n <- g2
  } else {
    g1 <- grids$g1
    g2 <- grids$g2
    grid.size <- length(g1)
    g1.n <- as.numeric(gsub('%', '', names(g1)))
    g2.n <- as.numeric(gsub('%', '', names(g2)))
    
    g1.n <- seq(0, 1, length.out=length(g1))
    g2.n <- seq(0, 1, length.out=length(g2))
  }
  
  # Evaluate responses over grid based on RF hyperrectangles
  grid <- matrix(0, nrow=grid.size, ncol=grid.size)
  for (i in seq(1, nrow(rectangles), by=2)) {

    rtemp <- rectangles[c(i, i+1),]
    tt1 <- rtemp$splits[rtemp$vars == intp1]
    tt2 <- rtemp$splits[rtemp$vars == intp2]
    
    if (rtemp$dirs[rtemp$vars == intp1] == 0) {
      idcs.1 <- g1 <= tt1
      x1 <- x[,intp1] <= tt1
    } else {
      idcs.1 <- g1 >= tt1
      x1 <- x[,intp1] >= tt1
    }
    
    if (rtemp$dirs[rtemp$vars == intp2] == 0) {
      idcs.2 <- g2 <= tt2
      x2 <- x[,intp2] <= tt2
    } else {
      idcs.2 <- g2 >= tt2
      x2 <- x[,intp2] >= tt2
    }
  
    if (pred.prob) {
      yy <- rectangles$prediction[i]
      grid[idcs.1, idcs.2] <- grid[idcs.1, idcs.2] + yy * rectangles$size.node[i]
    } else {
      # Evaluate average response in partitioned space
      if (sum(x1 & x2) < min.node.size) next
    
      grid[idcs.1, idcs.2] <- grid[idcs.1, idcs.2] + 
        ifelse(any(x1 & x2), mean(y[x1 & x2]), 0) * rectangles$size.node[i]
      grid[!idcs.1, idcs.2] <- grid[!idcs.1, idcs.2] + 
        ifelse(any(!x1 & x2), mean(y[!x1 & x2]), 0) * rectangles$size.node[i]
      grid[idcs.1, !idcs.2] <- grid[idcs.1, !idcs.2] +
        ifelse(any(x1 & !x2), mean(y[x1 & !x2]), 0) * rectangles$size.node[i]
      grid[!idcs.1, !idcs.2] <- grid[!idcs.1, !idcs.2] + 
        ifelse(any(!x1 & !x2), mean(y[!x1 & !x2]), 0) * rectangles$size.node[i]
    }
  }
  
  # rescale surface for node size 
  grid <- 2 * (grid / sum(rectangles$size.node[rectangles$size.node >= min.node.size]))
  
  if (is.null(range.color)) range.color <- range(grid)
  palette <- colorRampPalette(col.pal)
  colors <- palette(n.cols)
  facet.col <- cut(c(range.color, as.vector(grid)), n.cols)[-seq(2)]

  persp3d(x=g1.n, y=g2.n, z=grid, 
          xlab=xlab,
          ylab=ylab,
          zlab=zlab, zlim=z.range, 
          col=colors[facet.col],
          axes=axes)
}

quantileGrid <- function(x, grid.size, interact) {
  stopifnot(length(interact) == 2)
  grid.size <- grid.size + 1
  grids <- list()
  grids$g1 <- quantile(x[,interact[1]], probs=seq(0, 1, length.out=grid.size))
  grids$g2 <- quantile(x[,interact[2]], probs=seq(0, 1, length.out=grid.size))
  return(grids)
}
