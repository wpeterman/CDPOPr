## Internal functions for processing CDPOP files

# Import CDPOP ------------------------------------------------------------

read.grid <- function(grid,
                      pops = NULL){
  suppressWarnings(
    cdpop_out <- read_csv(grid,
                          col_types = cols(Subpopulation = col_skip(),
                                           XCOORD = col_skip(), YCOORD = col_skip(),
                                           sex = col_skip(), age = col_skip(),
                                           infection = col_skip(), DisperseCDist = col_skip(),
                                           hindex = col_skip()))
  )

  occ_pop <- which(cdpop_out$ID != "OPEN")

  if(!is.null(pops)) {
    return(occ_pop)
  } else {

    cd_df <- as.data.frame(cdpop_out[occ_pop,-1])

    ncode <- 1
    gi <- adegenet::df2genind(cd_df,
                              ncode = ncode)
    return(gi)
  }
}


# CDPOP to adegenet -------------------------------------------------------

cdpop2genind <- function(cdpop_out,
                         ncode = 1){

}



# Sample gi ---------------------------------------------------------------


## Randomly select populations and individuals from within populations

gi_samp <- function(gi,
                    n_ind = 100) {
  ind_samp <- sort(sample(1:n_ind(gi), n_ind))
  gi_s <- gi[ind_samp]
  # gi_s <- lapply(gi_s, function(x) x[sample(1:nrow(x$tab), n_ind)])
  # gi_s <- repool(gi_s)

  out <- list(genind = gi_s,
              pop_samp = ind_samp)
}
