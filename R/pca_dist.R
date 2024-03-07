#' Function to calculate PCA distance between individuals
#' @description Calculate PCA distance
#'
#' @param gi Object of class `genind` from the `adegenet` package
#' @param n_axes Number of principal component axes to retain to calculate pairwise distances
#' @param scale (Default = TRUE) Calculated distance will be rescaled 0-1
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' ## Not run:
#'df <- data.frame(locusA = c("11","11","12","32"),
#'                 locusB = c("37","34","55","15"),
#'                 locusC = c("22","22","21","22"))
#'df
#'
#'obj <- adegenet::df2genind(df,
#'                           ploidy = 2,
#'                           ncode = 1)
#'
#'pca_dist(gi = obj,
#'         n_axes = 4)
#'
#'pca_dist(gi = obj,
#'         n_axes = 4,
#'         scale = FALSE)
#' ## End (Not run)

pca_dist <- function(gi,
                     n_axes = 16,
                     scale = TRUE,
                     dudi.pca = TRUE){
  a_tab <- adegenet::tab(gi)
  if(n_axes > dim(a_tab)[1] | n_axes > dim(a_tab)[2]){
    stop("Number of axes must be less than the number of individuals and loci in the data set.")
  }

  if(isTRUE(dudi.pca)){
    gi_scale <- scaleGen(gi)
    pc <- dudi.pca(gi_scale,
                   cent = FALSE,
                   scale = FALSE,
                   scannf = FALSE,
                   nf = n_axes)
    pc_ <- pc$li
  } else {
    pc <- prcomp(a_tab)
    pc_ <- pc$x
  }

  if(isTRUE(scale)){
    pc_dist <- as.matrix(scale_01(dist(pc_[,1:n_axes])))

  } else {
    pc_dist <- as.matrix(dist(pc_[,1:n_axes]))
  }

  return(pc_dist)
}
