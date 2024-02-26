#' Scale vector of values to range 0-1
#' @author Bill Peterman
#' @title Scale vector
#' @description Scale vector to range from 0-1
#'
#' @param x Vector to be scaled
#' @export
#' @examples
#' x <- runif(5, -10, 10)
#' scale_01(x)
scale_01 <- function(x) { #scales a given vector so that min(x)==0 and max(x)==1
  (x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))
}
