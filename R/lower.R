#' Make a vector of the lower half of a square distance matrix
#'
#' Extract lower half of distance matrix
#'
#' @param matrix Square distance matrix with no row names.
#' @return A vector of the lower half of the matrix
#'
#' @details This is a convenience function to obtain the lower half of a matrix

#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' ## Not run:
#' ## m <- matrix(rnorm(25), 5)
#' ## l_m <- lower(m)
#'
#' ## End (Not run)

lower <- function(matrix) {
  if (is.vector(matrix) == TRUE ||
      dim(matrix)[1] != dim(matrix)[2]) {
    warning("Must provide square distance matrix with no column or row names")
  }
  lm <- matrix[lower.tri(matrix)]
  return(lm)
}
