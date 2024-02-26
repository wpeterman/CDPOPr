#' Function to facilitate the running of CDPOP software from R
#' @description Function to run CDPOP from R
#'
#' @param sim_name Name for simulation results. Defaults to 'output'
#' @param pts Spatial points object
#' @param sim_dir Directory where simulation results will be written. Be sure there are no spaces or dashes in this directory path.
#' @param resist_rast Resistance surface
#' @param agefilename Path to age file. Default will create and use a non-overlapping generations file.
#' @param mcruns Default = 1
#' @param looptime Default = 401; Number generations to conduct simulation
#' @param output_years Default = 50; Interval to write out simulation results
#' @param gridformat Default = 'genepop'; c('genepop', 'genalex', 'structure', 'cdpop')
#' @param output_unicor Default = 'N'. ‘Y’ – create X,Y files corresponding to the output_years in CSV format compatible with the landscape connectivity program UNICOR.
#' @param cdclimgentime Default = 0. To initiate the CDClimate module, this is the generation/year that the next effective distance matrix will be read in at. You can specify multiple generations by separating each generation to read in the next cost distance matrix by ‘|’.  Then in the following surface columns, a separate file can be given for each generation.
#' @param matemoveno Default = 2; Uses Inverse Square (1 / (Cost Distance^2)). This function gets rescaled to min and threshold of the inverse square cost distance.
#' @param matemoveparA Not used with inverse square movement
#' @param matemoveparB Not used with inverse square movement
#' @param matemoveparC Not used with inverse square movement
#' @param matemovethresh Default = 'max'; The maximum movement is the maximum resistance distance
#' @param output_matedistance Should mate distance matrix be saved. Default = 'N'
#' @param sexans Default = 'N'; No selfing
#' @param Freplace Default = 'Y'; Females mate without replacement
#' @param Mreplace Default = 'Y'; Males mate without replacement
#' @param philopatry Default = 'N';
#' @param multiple_paternity Default = 'N'; No philopatry of Males or Females
#' @param selfans Default = 'N'; No selfing
#' @param mateFrequency Default = 1. Probability with which the female will mate. E.g., 0.5 would mean that a mature female will mate every other year.
#' @param Fdispmoveno Default = NULL. Will be set equal to matemoveno
#' @param FdispmoveparA Not used with inverse square movement
#' @param FdispmoveparB Not used with inverse square movement
#' @param FdispmoveparC Not used with inverse square movement
#' @param Fdispmovethresh Default = NULL. Will be set to matemovethresh
#' @param Mdispmoveno Default = NULL. Will be set equal to matemoveno
#' @param MdispmoveparA Not used with inverse square movement
#' @param MdispmoveparB Not used with inverse square movement
#' @param MdispmoveparC Not used with inverse square movement
#' @param Mdispmovethresh Default = NULL. Will be set to matemovethresh
#' @param offno Default = 2; Poisson draw around ‘mean fecundity’
#' @param MeanFecundity Default = 5; Specifies mean fecundity in age variable file.
#' @param Femalepercent Default = 50
#' @param EqualsexratioBirth Default = 'N'
#' @param TwinningPercent Default = 0
#' @param popModel Default = 'exp'
#' @param r Population growth rate. No applicable when using exponential growth rate
#' @param K_env Equal to the number of individuals simulated
#' @param subpopmortperc Default = 0|0|0|0; Not using subpopulation features
#' @param muterate Default = 0.0005
#' @param mutationtype Default = 'forward'
#' @param loci Default = 1000; For simulating SNP-like markers
#' @param intgenesans Default = 'random'; Random initiation of alleles
#' @param allefreqfilename Default = 'N'
#' @param alleles Default = 2; For simulating SNP-like markers
#' @param mtdna Default = 'N'
#' @param startGenes Default = 0;
#' @param cdevolveans Default = 'N'; No loci are under selection
#' @param startSelection Default = 0; No selection
#' @param betaFile_selection Default = 'N'; No selection
#' @param epistasis Default = 'N'; No epigenetics
#' @param epigeneans Default = 'N'; No epigenetics
#' @param startEpigene Default = 0; No epigenetics
#' @param betaFile_epigene Default = 'N'; No epigenetics
#' @param cdinfect Default = 'N'; No epigenetics
#' @param transmissionprob Default = 0; No epigenetics
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' ## Not run:
#' ## *** TO BE COMPLETED *** ##
#'
#' ## End (Not run)

cdpop <- function(sim_name = 'output_',
                       pts,
                       sim_dir,
                       resist_rast,
                       agefilename = NULL,
                       mcruns = 1,
                       looptime = 400,
                       output_years = 50,
                       gridformat = 'genepop',
                       output_unicor = 'N',
                       cdclimgentime = 0,
                       matemoveno = 2,
                       matemoveparA = 0,
                       matemoveparB = 0,
                       matemoveparC = 0,
                       matemovethresh = 'max',
                       output_matedistance = 'N',
                       sexans = 'Y',
                       Freplace = 'Y',
                       Mreplace = 'Y',
                       philopatry = 'N',
                       multiple_paternity = 'N',
                       selfans = 'N',
                       mateFrequency = 1,
                       Fdispmoveno = NULL,
                       FdispmoveparA = 0,
                       FdispmoveparB = 0,
                       FdispmoveparC = 0,
                       Fdispmovethresh = NULL,
                       Mdispmoveno = NULL,
                       MdispmoveparA = 0,
                       MdispmoveparB = 0,
                       MdispmoveparC = 0,
                       Mdispmovethresh = NULL,
                       offno = 2,
                       MeanFecundity = 5,
                       Femalepercent = 50,
                       EqualsexratioBirth = 'N',
                       TwinningPercent = 0,
                       popModel = 'exp',
                       r = 1,
                       K_env = length(pts),
                       subpopmortperc = 0,
                       muterate = 0.0005,
                       mutationtype = 'random',
                       loci = 1000,
                       intgenesans = 'random',
                       allefreqfilename = 'N',
                       alleles = 2,
                       mtdna = 'N',
                       startGenes = 0,
                       cdevolveans = 'N',
                       startSelection = 0,
                       betaFile_selection = 'N',
                       epistasis = 'N',
                       epigeneans = 'N',
                       startEpigene = 0,
                       betaFile_epigene = 'N',
                       cdinfect = 'N',
                       transmissionprob = 0){

  # CDPOP Setup -------------------------------------------------------------

  CDPOP.py <- cdpop_setup()

  # Create directories ------------------------------------------------------

  if(!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
  suppressWarnings(
    dir.create(paste0(sim_dir,"/data/"), recursive = TRUE)
  )
  data_dir <- paste0(sim_dir,"/data/")


  # Fill NULL ---------------------------------------------------------------

  if(matemoveno == 9){
    if(class(resist_rast) == "RasterLayer"){
      stop('Specify a probability matrix instead of a raster layer!')
    }
    write.table(resist_rast,
                paste0(data_dir, "move_prob.csv"),
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)

    cdmat <- 'move_prob'
    # write.table(matemoveno, paste0(data_dir, 'DispProb.csv'),
    #             sep = ",",
    #             row.names = F)
    # matemoveno <- 'DispProb'
  }

  if(is.null(Fdispmoveno)){
    Fdispmoveno <- matemoveno
  }

  if(is.null(Mdispmoveno)){
    Mdispmoveno <- matemoveno
  }

  if(is.null(Fdispmovethresh)){
    Fdispmovethresh <- matemovethresh
  }

  if(is.null(Mdispmovethresh)){
    Mdispmovethresh <- matemovethresh
  }

  # Age file ----------------------------------------------------------------

  if(is.null(agefilename)){
    age_df <- data.frame(`Age class` = c(0,1),
                         Distribution = c(0,1),
                         `Male Mortality` = c(0,100),
                         `Female Mortality` = c(0,100),
                         `Mean Fecundity` = c(0,MeanFecundity),
                         `Std Fecundity` = c(0,0),
                         `Male Maturation` = c(0,1),
                         `Female Maturation` = c(0,1),
                         check.names = F)
    write.table(age_df, paste0(data_dir, 'AgeVars.csv'),
                sep = ",",
                row.names = F)
    # age_file <- paste0(data_dir, 'AgeVars.csv')
    age_file <- 'AgeVars.csv'

  }


  # XY File -----------------------------------------------------------------

  xyFile_df <- data.frame(Subpopulation = rep(1, length(pts)),
                          XCOORD = pts@coords[,1],
                          YCOORD = pts@coords[,2],
                          Subpop_mortperc = subpopmortperc,
                          ID = paste0('initial',1:length(pts) - 1),
                          sex = sample(c(0,1),
                                       replace = T,
                                       size = length(pts)),
                          Fitness_AA = rep(0, length(pts)),
                          Fitness_Aa = rep(0, length(pts)),
                          Fitness_aa = rep(0, length(pts)),
                          Fitness_AABB = rep(0, length(pts)),
                          Fitness_AaBB = rep(0, length(pts)),
                          Fitness_aaBB = rep(0, length(pts)),
                          Fitness_AABb = rep(0, length(pts)),
                          Fitness_AaBb = rep(0, length(pts)),
                          Fitness_aaBb = rep(0, length(pts)),
                          Fitness_AAbb = rep(0, length(pts)),
                          Fitness_Aabb = rep(0, length(pts)),
                          Fitness_aabb = rep(0, length(pts))
  )

  write.table(xyFile_df, paste0(data_dir, 'xyFile.csv'),
              sep = ',',
              row.names = F)
  # xyFile <- paste0(data_dir, 'xyFile')
  xyFile <- 'xyFile'


  # # Convert points to file --------------------------------------------------
  #
  # pts.tab <- data.frame(site = 1:length(pts),
  #                       x = pts@coords[,1],
  #                       y = pts@coords[,2])
  #
  # write.table(pts.tab,
  #             paste0(sim_dir, "xy.csv"),
  #             sep = ",",
  #             row.names = FALSE,
  #             col.names = TRUE)


  # Resistance Distance -----------------------------------------------------
  if(matemoveno != 9){
    print("Calculating resistance distance with `gdistance`...")

    trans <- transition(x = resist_rast,
                        transitionFunction = function(x)  1 / mean(x),
                        directions = 8)

    trR <- geoCorrection(trans, "r", scl = T)
    resist_mat <- as.matrix(commuteDistance(trR, pts) / 1000)

    ## Check file format, row/col names?
    write.table(resist_mat,
                paste0(data_dir, "resist_mat.csv"),
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)

    # cdmat <- paste0(data_dir, "resist_mat.csv")
    cdmat <- 'resist_mat'
  }

  # CDPOP input ----------------------------------------------------------

  cdpop_df <- data.frame(xyfilename = xyFile,
                         agefilename = age_file,
                         mcruns = mcruns,
                         looptime = looptime,
                         output_years = output_years,
                         gridformat = gridformat,
                         output_unicor = output_unicor,
                         cdclimgentime = cdclimgentime,
                         matecdmat = cdmat,
                         dispcdmat = cdmat,
                         matemoveno = matemoveno,
                         matemoveparA = matemoveparA,
                         matemoveparB = matemoveparB,
                         matemoveparC = matemoveparC,
                         matemovethresh = matemovethresh,
                         output_matedistance = output_matedistance,
                         sexans = sexans,
                         Freplace = Freplace,
                         Mreplace = Mreplace,
                         philopatry = philopatry,
                         multiple_paternity = multiple_paternity,
                         selfans = selfans,
                         mateFrequency = mateFrequency,
                         Fdispmoveno = Fdispmoveno,
                         FdispmoveparA = FdispmoveparA,
                         FdispmoveparB = FdispmoveparB,
                         FdispmoveparC = FdispmoveparC,
                         Fdispmovethresh = Fdispmovethresh,
                         Mdispmoveno = Mdispmoveno,
                         MdispmoveparA = MdispmoveparA,
                         MdispmoveparB = MdispmoveparB,
                         MdispmoveparC = MdispmoveparC,
                         Mdispmovethresh = Mdispmovethresh,
                         offno = offno,
                         Femalepercent = Femalepercent,
                         EqualsexratioBirth = EqualsexratioBirth,
                         TwinningPercent = TwinningPercent,
                         popModel = popModel,
                         r = r,
                         K_env = K_env,
                         subpopmortperc = subpopmortperc,
                         muterate = muterate,
                         mutationtype = mutationtype,
                         loci = loci,
                         intgenesans = intgenesans,
                         allefreqfilename = allefreqfilename,
                         alleles = alleles,
                         mtdna = mtdna,
                         startGenes = startGenes,
                         cdevolveans = cdevolveans,
                         startSelection = startSelection,
                         betaFile_selection = betaFile_selection,
                         epistasis = epistasis,
                         epigeneans = epigeneans,
                         startEpigene = startEpigene,
                         betaFile_epigene = betaFile_epigene,
                         cdinfect = cdinfect,
                         transmissionprob = transmissionprob,
                         check.names = F)
  # colnames(cdpop_df) <- input_names[[1]]

  write.table(cdpop_df,
              paste0(data_dir, "CDPOP_inputs.csv"),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE,
              quote = F)


  # Run CDPOP ---------------------------------------------------------------
  print("Running CDPOP...")

  system(paste("python", CDPOP.py, data_dir, "CDPOP_inputs.csv", sim_name))


  # Import Results ----------------------------------------------------------

  fi <- file.info(list.files(path = sim_dir,
                             pattern = "grid",
                             recursive = T,
                             full.names = T))

  ## Get latest simulation results
  newest_sim <- dirname(rownames(fi)[which.max(fi$mtime)])

  grid_dir <- list.files(path = newest_sim,
                         pattern = "grid",
                         recursive = T,
                         full.names = T)

  grid_list <- lapply(grid_dir, read.grid)
  pop_list <- lapply(grid_dir, read.grid, pops = TRUE)

  gens <- basename(grid_dir) %>% sub('.csv', '', .) %>%
    sub('grid', '',.) %>% as.numeric()

  grid_list <- grid_list[order(gens)]
  pop_list <- pop_list[order(gens)]

  names(pop_list) <- names(grid_list) <- paste0('gen_', sort(gens))


  # Wrap-up -----------------------------------------------------------------
  out <- list(grid_list = grid_list,
              pop_list = pop_list)

  # unlink(paste0(tempdir(),"\\cdpop"), recursive = T)
  return(out)
}
