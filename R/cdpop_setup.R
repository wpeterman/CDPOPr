#' Execute to run CDPOP from R
#' @description Function to format CDPOP to run with R
#'
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' ## Not run:
#' ## cdpop_setup()
#'
#' ## End (Not run)

cdpop_setup <- function(){
  if(dir.exists(paste0(tempdir(),"\\cdpop"))){
    to <- normalizePath(paste0(tempdir(),"\\cdpop\\"))

  } else {
    from <- system.file("cdpop", package="CDPOPr")
    dir.create(paste0(tempdir(),"\\cdpop"), recursive = T, showWarnings = F)
    to <- normalizePath(paste0(tempdir(),"\\cdpop\\"))

    file.copy(list.files(from,
                         # recursive = TRUE,
                         full.names = TRUE),
              to,
              recursive = TRUE,
              overwrite = T)
  }

  CDPOP.py <- normalizePath(list.files(to,
                                       full.names = T,
                                       recursive = T,
                                       pattern = 'CDPOP.py'))
  return(CDPOP.py)
}
