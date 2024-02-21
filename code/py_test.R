

# Make Python function ----------------------------------------------------
## Basic run
system(paste("python C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/ inputvars.csv output_test"))

## Test with all files in same directory
system(paste("python C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py C:/Users/peterman.73/Box/R/CDPOP/R_Test/ XX_inputvars.csv test_"))

system(paste("python C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/ inputvars_mod.csv test_"))

input.csv <- read.csv('C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/inputvars.csv')
agefile <- read.csv('C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/agevars/Agevars_nonOverlap.csv')
xyFile <- read.csv('C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/xyfiles/xyED16.csv')

CDPOP.py <- 'C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py'

sim_name <- 'output'
n.pops <- 25
set.seed(123)
pts <- SpatialPoints(cbind(runif(n.pops, 20,80), runif(n.pops, 20,80)))
sim_dir <- "C:/Users/peterman.73/Box/R/CDPOP/test_sim"
resist_rast <- exp(NLMR::nlm_gaussianfield(100, 100) * 1.5)
plot(resist_rast)
plot(pts, add = T, pch = 19)

res <- readLines('C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/inputvars.csv')[1]
input_names <- strsplit(res, split = ',')

test_in <- read.csv('C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/inputvars_mod.csv',
                    as.is = T)
write.table(test_in,
            'C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/inputvars_out.csv',
            sep = ",", row.names = F)
writeLines(as.character(cdpop_df), 'C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/inputvars_out.csv'
)

# Compare files -----------------------------------------------------------

fn_cdpop <- read.csv("C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/CDPOP_inputs.csv")
names(fn_cdpop) == names(input.csv)

fn_xy <- read.csv("C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/xyFile.csv")
names(fn_xy) == names(xyFile)

# Function params ---------------------------------------------------------

agefilename = NULL
mcruns = 1
looptime = 400
output_years = 50
gridformat = 'genepop'
cdclimgentime = 0
matemoveno = 2
matemoveparA = 0
matemoveparB = 0
matemoveparC = 0
matemovethresh = '20max'
output_matedistance = 'N'
sexans = 'Y'
Freplace = 'N'
Mreplace = 'N'
philopatry = 'N'
multiple_paternity = 'N'
selfans = 'N'
Fdispmoveno = 2
FdispmoveparA = 0
FdispmoveparB = 0
FdispmoveparC = 0
Fdispmovethresh = '20max'
Mdispmoveno = 2
MdispmoveparA = 0
MdispmoveparB = 0
MdispmoveparC = 0
Mdispmovethresh = '20max'
offno = 2
MeanFecundity = 10
Femalepercent = 50
EqualsexratioBirth = 'AtBirth'
TwinningPercent = 0
popModel = 'exp'
r = 1
K_env = length(pts)
subpopmortperc = 0
muterate = 0.0005
mutationtype = 'random'
loci = 1000
intgenesans = 'random'
allefreqfilename = 'N'
alleles = 2
mtdna = 'N'
startGenes = 0
cdevolveans = 'N'
startSelection = 0
betaFile_selection = 'N'
epistasis = 'N'
epigeneans = 'N'
startEpigene = 0
betaFile_epigene = 'N'
cdinfect = 'N'
transmissionprob = 0


system(paste("python", CDPOP.py, data_dir, "inputvars_mod.csv", sim_name))
system(paste("python", CDPOP.py, data_dir, "inputvars_out.csv", sim_name))


# Test Function -----------------------------------------------------------

source('C:/Users/peterman.73/Box/R/CDPOP/cdpop_function.R', encoding = 'UTF-8')

sim_name <- 'short_'
n.pops <- 10
# set.seed(123)
pts <- SpatialPoints(cbind(runif(n.pops, 20,80), runif(n.pops, 20,80)))
sim_dir <- "C:/Users/peterman.73/Box/R/CDPOP/test_sim"
resist_rast <- exp(NLMR::nlm_gaussianfield(100, 100) * 1.5)
plot(resist_rast)
plot(pts, add = T, pch = 19)


cdpop(CDPOP.py = 'C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py',
      sim_name = sim_name,
      pts = pts,
      resist_rast = resist_rast,
      sim_dir = "C:/Users/peterman.73/Box/R/CDPOP/test_sim",
      looptime = 25,
      gridformat = 'cdpop',
      loci = 5,
      alleles = 2,
      Fdispmovethresh = 'max',
      Mdispmovethresh = 'max',
      matemovethresh = 'max')

cdpop_out <- read_csv("test_sim/data/output_1582908946/batchrun0mcrun0/grid11.csv", 
                      col_types = cols(Subpopulation = col_skip(), 
                                       XCOORD = col_skip(), YCOORD = col_skip(), 
                                       sex = col_skip(), age = col_skip(), 
                                       infection = col_skip(), DisperseCDist = col_skip(), 
                                       hindex = col_skip()))
cd.genind <- df2genind(cdpop_out)
