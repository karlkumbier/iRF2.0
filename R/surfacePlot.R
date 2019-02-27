plotInt <- function(x, y, int, 
                    varnames=NULL,
                    read.forest=NULL,
                    qcut=0.5,
                    col.pal=c('#1c3f66', '#306aab', 
                      '#6e96c4', '#ffb003',  '#ff8300'), 
                    xlab=NULL, ylab=NULL, zlab=NULL,
                    range.color=NULL,
                    grid.size=50,
                    min.surface=100,
                    min.node=5,
                    pred.prob=FALSE,
                    main=NULL,
                    plot.dir=NULL) {
  # Generates surface plots for first two interacting features by levels 
  # of the remaining features
  # args:
  #   x: data matrix, with replicate features grouped
  #   y: response vector, binary or continuous
  #   int: signed interaction to plot, formatted as <xi>_<xj>  
  #   varnames: vector of feature names corresponding to columns of x
  #   read.forest: list output of readForest, required if rectangles is null
  #   qcut: quantile to define low/high levels of additional features 
  #     (when order int > 2).
  #   xlab, ylab, zlab: axis labels
  #   min.surface: minimum number of observations required to generate surface 
  #     map plot.
  #   min.node: minimum leaf node size to use for grid
  #   pred.prob: T/F indicating whether the z axis should indicate predicted
  #     probability from the forest or raw data distribution
  #   plot.dir: directory to write plots to
  if (! 'rgl' %in% rownames(installed.packages()))
    stop('Surface map plots require rgl installation')

  if (is.null(varnames)) {
    if (is.null(colnames(x))) varnames <- 1:ncol(x)
    else varnames <- colnames(x)
  }
  
  int <- str_split(int, '_')[[1]]
  int.nf <- int2Id(int, varnames, split=TRUE, signed=TRUE)
  int.x <- int.nf %% ncol(x) + ncol(x) * (int.nf %% ncol(x) == 0)
  
  if (length(int) > 2) {
    # Evaluate thresholds for high/low levels of additional features
    nf <- read.forest$node.feature
    int.id <- Matrix::rowMeans(nf[,int.nf] != 0) == 1
    qCut <- function(x) quantile(x, probs = qcut)
    thr <- apply(nf[int.id, int.nf], MAR=2, qCut)
    sgn <- intSign(int)
    
    # Determine which plot each observation is assigned to
    id.plot <- 3:length(int)
    int.clean <- str_remove_all(int, '[-\\+]')
    xx <- t(x[,int.x[id.plot]]) * sgn[id.plot]
    thr[id.plot] <- thr[id.plot] * sgn[id.plot]
    id <- xx > thr[id.plot]
    id <- apply(id, MAR=2, paste, collapse='_')
  } else {
    int.clean <- str_remove_all(int, '[-\\+]')
    id <- rep(1, nrow(x))
    id.plot <- 1:2
  }
  
  # Generate grid of x/y values for surface maps
  grids <- quantileGrid(x, grid.size, int.x[1:2])
  
  # Extract hyperrectangles from read forest output
  rectangles <- forestHR(read.forest, int.nf, min.node)
  
  # Generate surface maps for each group of observations
  ids <- lapply(unique(id), '==', id)
  surfaces <- lapply(ids, function(ii) {
    if (sum(ii) < min.surface) return(NULL)
    genSurface(x[ii,], y[ii], int.nf[1:2], varnames=varnames, 
               rectangles=rectangles, min.node=min.node, grids=grids)
  })
  
  # Set color range for surface map plots
  if (is.null(range.color)) range.color <- range(unlist(surfaces))

  # Iterate over observation groups to generate response surfaces
  for (i in 1:length(unique(id))) {
    
    open3d()
    if (is.null(surfaces[[i]])) next
    
    # Generate title for surface map
    ii <- unique(id)[i]
    ii <- str_replace_all(ii, 'TRUE', 'High')
    ii <- str_replace_all(ii, 'FALSE', 'Low')
    if (length(int) > 2 & !is.null(main)) {
      group <- paste(int.clean[3:length(int)], ii, sep='-')
      main <- paste0(main, '; ', group)
    }
    
    # Generate response surface for curent group
    plotInt2(surfaces[[i]], xlab=xlab, ylab=ylab, zlab=zlab, main=main,
             col.pal=col.pal, range.color=range.color)
    
    # Write plot to output directory
    if (!is.null(plot.dir)) {
      pp <- paste(int.clean[id.plot], collapse='_')
      sub.dir <- paste0(plot.dir, pp, '-', ii)
      dir.create(sub.dir, recursive=TRUE, showWarnings=FALSE)
      par3d(windowRect = c(20, 30, 1000, 1000))
      rgl.viewpoint(zoom=0.85, theta=0, phi=-75)
      movie3d(spin3d(axis = c(0,0,1), rpm = 10), duration=6,  type="png", 
              dir=sub.dir, convert=FALSE, clean=FALSE)
    }
  }
}

plotInt2 <- function(surface,
                     col.pal=c('#1c3f66', '#306aab', 
                       '#6e96c4', '#ffb003',  '#ff8300'), 
                     xlab=NULL, ylab=NULL, zlab=NULL,
                     main=NULL,
                     range.color=NULL,
                     z.range=c(0,1),
                     n.cols=100, 
                     axes=TRUE) {
  # Generates surface map plot of order-2 interaction
  # args:
  #   surface: response surface matrix, output of genSurface
  #   col.pal: color palette of surface map
  #   xlab, ylab, zlab: axis labels
  #   range.color: range of color values
  #   z.range: range for response axis
  #   n.cols: number of colors in color pal
  #   grids: user generated grid to plot over. If supplied, grid.size is 
  #     ignored
  #   axes: T/F indicating whether axes should be plotted
  range.color <- c(range.color, as.vector(surface))
  if (length(unique(range.color)) == 1) n.cols <- 1
  
  palette <- colorRampPalette(col.pal)
  colors <- palette(n.cols)
  facet.col <- cut(range.color, n.cols)[-seq(2)]

  # Plot interaction response surface
  par3d(cex=1.5)
  if (is.null(zlab)) zlab <- ''
  g1n <- as.numeric(rownames(surface))
  g2n <- as.numeric(colnames(surface))
  persp3d(x=g1n, y=g2n, z=surface, xlab=xlab, ylab=ylab, zlab=zlab, 
          zlim=z.range, col=colors[facet.col], axes=axes, main=main)
  
}

genSurface <- function(x, y, int, 
                       varnames=NULL,
                       rectangles=NULL,
                       read.forest=NULL,
                       grid.size=100, 
                       grids=NULL, 
                       pred.prob=FALSE, 
                       min.node=5) {
  # Generates surface map of order-2 interaction
  # args:
  #   x: design matrix
  #   y: response vector
  #   int: order-2 signed interaction to plot, formatted as <xi>_<xj>
  #   varnames: character vector indicating feature grouping
  #   rectangles: hyprerectangle list as generated by forestHR
  #   read.forest: list output of readForest, required if rectangles is null
  #   grid.size: number of bins to plot surface map over
  #   col.pal: color palette of surface map
  #   xlab, ylab, zlab: axis labels
  #   range.color: range of color values
  #   z.range: range for response axis
  #   n.cols: number of colors in color pal
  #   grids: user generated grid to plot over. If supplied, grid.size is 
  #     ignored
  #   axes: T/F indicating whether axes should be plotted
  #   pred.prob: T/F indicating whether the z axis should indicate predicted
  #     probability from the forest or raw data distribution
  #   min.node: minimum leaf node size to use for grid
  stopifnot(!is.null(rectangles) | !is.null(read.forest))
  
  if (is.null(varnames)) {
    if (is.null(colnames(x))) varnames <- 1:ncol(x)
    else varnames <- colnames(x)
  }
  
  if (is.factor(y)) y <- as.numeric(y) - 1
  n <- nrow(x)
  p <- ncol(x)
  
  stopifnot(length(int) == 2) 
  i1 <- int[1] %% p + p * (int[1] %% p == 0)
  i2 <- int[2] %% p + p * (int[2] %% p == 0)
  
  # generate grid to plot surface over either as raw values or quantiles
  if (is.null(grids)) {
    g1 <- seq(min(x[,i1]), max(x[,i1]), length.out=grid.size)
    g2 <- seq(min(x[,i2]), max(x[,i2]), length.out=grid.size)
    g1n <- round(g1, 2)
    g2n <- round(g2, 2)
  } else {
    g1 <- grids$g1
    g2 <- grids$g2
    g1n <- seq(0, 1, length.out=length(g1))
    g2n <- seq(0, 1, length.out=length(g2))
    grid.size <- length(g1)
  }
  
  # Evaluate responses over grid based on RF hyperrectangles
  if (is.null(rectangles)) rectangles <- forestHR(read.forest, int, min.node)
  grid <- matrix(0, nrow=grid.size, ncol=grid.size)
  removed <- rep(FALSE, nrow(rectangles))
  for (i in 1:nrow(rectangles)) {
    
    rcur <- rectangles[i,]
    wt <- rectangles$size.node[i]
    tt1 <- unlist(rcur$splits)[unlist(rcur$vars) == int[1]]
    tt2 <- unlist(rcur$splits)[unlist(rcur$vars) == int[2]]
    
    # Evalaute which observations/grid elements correspond to current HR
    sgn <- ifelse(int[1] <= p, -1, 1)
    idcs1 <- sgn * g1 >= sgn * tt1
    x1 <- sgn * x[,i1] >= sgn * tt1
    
    sgn <- ifelse(int[2] <= p, -1, 1)
    idcs2 <- sgn * g2 >= sgn * tt2
    x2 <- sgn * x[,i2] >= sgn * tt2
    
    if (pred.prob) {
      # Evaluate RF predictions for region corresponding to current HR
      yy <- rectangles$prediction[i]
      grid[idcs1, idcs2] <- grid[idcs1, idcs2] + yy * wt
    } else {
      # If a region contains no observations, move to next hyperrectangle
      if (!any(x1 & x2) | !any(!x1 & !x2) | !any(x1 & !x2) | !any(!x1 & x2)) {
        removed[i] <- TRUE
        next
      }
      
      # Evaluate average response value in regions corresponding to HR
      wt <- rectangles$size.node[i]
      
      grid[idcs1, idcs2] <- grid[idcs1, idcs2] +  mean(y[x1 & x2]) * wt
      grid[!idcs1, idcs2] <- grid[!idcs1, idcs2] +  mean(y[!x1 & x2]) * wt
      grid[idcs1, !idcs2] <- grid[idcs1, !idcs2] +  mean(y[x1 & !x2]) * wt
      grid[!idcs1, !idcs2] <- grid[!idcs1, !idcs2] +  mean(y[!x1 & !x2]) * wt
    }
  }
  
  # rescale surface for node size and generate corresponding color palette
  grid <- grid / sum(rectangles$size.node[!removed])
  if (all(grid == 0)) grid <- grid + 1e-3
  rownames(grid) <- g1n
  colnames(grid) <- g2n
  return(grid)
}


getPathsTree <- function(read.forest, int) {
  # Extract decision paths across all leaf nodes in a given tree.
  # args:
  #   read.forest: list as returned by readForest, including node.feature and
  #     tree.info entries
  #   int: a numeric vector indicating interaciton variables
  nf <- read.forest$node.feature
  out <- select(read.forest$tree.info, prediction, node.idx, tree, size.node)
  
  # Filter nodes to active interactions
  id <- Matrix::rowMeans(nf[,int] != 0) == 1
  out <- out[id,]
  nf <- nf[id,]

  # Extract splitting features and thresholds from decision paths
  d <- data.table(row=nf@i + 1, thr=nf@x,
                  col=findInterval(seq(nf@x) - 1, nf@p[-1]) + 1) %>%
    group_by(row) %>%
    summarize(vars=list(col), splits=list(thr))

  out$vars <- d$vars
  out$splits <- d$splits
  return(out)
}

forestHR<- function(read.forest, int, min.node=1) {
  # Read hyperrectangles from RF for a specified interactin
  # args:
  #   read.forest: list as returned by readForest, including node.feature and
  #     tree.info entries
  #   varnames.group: character vector indicating feature grouping
  #   int: vector of features for which to generate hyperrectangles
  #   min.node: filter all leaf nodes that are not larger than specified 
  #     value for faster processing
  idrm <- read.forest$tree.info$size.node < min.node
  read.forest <- subsetReadForest(read.forest, !idrm)
  out <- getPathsTree(read.forest, int) 
  
  return(out)
}

quantileGrid <- function(x, grid.size, int) {
  stopifnot(length(int) == 2)
  grid.size <- grid.size + 1
  grids <- list()
  grids$g1 <- quantile(x[,int[1]], probs=seq(0, 1, length.out=grid.size))
  grids$g2 <- quantile(x[,int[2]], probs=seq(0, 1, length.out=grid.size))
  return(grids)
}
