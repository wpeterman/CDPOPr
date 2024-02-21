
# Install Reticulate ------------------------------------------------------
## Directions here
# https://ttimbers.github.io/intro-to-reticulate/setup-instructions/setup-after-installing-python.html

## Packages

list.of.packages <- c("reticulate",
                      "png",
                      "adegenet",
                      "readr",
                      "terra",
                      "raster")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(reticulate)
py_available()

# Test Function -----------------------------------------------------------

source('code/cdpop_function.R', encoding = 'UTF-8')
library(ResistanceGA)

sim_name <- 'short_'
n.pops <- 500
r.dim <- 100
suit.hab <- 0.5
# set.seed(123)
pts <- unique(floor(cbind(runif(5000, 0.15 * r.dim, 0.85 * r.dim), 
                          runif(5000, 0.15 * r.dim, 0.85 * r.dim))))
# pts <- SpatialPoints(cbind(runif(n.pops, 25,225), runif(n.pops, 25,225)))
sim_dir <- "C:/Users/peterman.73/OneDrive - The Ohio State University/R/CDPOP/test_sim/"
resist_rast <- exp(NLMR::nlm_gaussianfield(r.dim, r.dim, rescale = F) * 0.5)
hab.thresh <- quantile(resist_rast, suit.hab)

sample.extract <- extract(resist_rast, pts)
sample.suit <- pts[sample.extract <= hab.thresh,]
sample.extract <- sample.extract[sample.extract <= hab.thresh]
# prob <- ResistanceGA:::SCALE.vector(1/sample.extract,0,1)
pts <- SpatialPoints(sample.suit[sample(nrow(sample.suit), n.pops, replace = F),])
# pts <- SpatialPoints(sample.suit[sample(nrow(sample.suit), n.pops, replace = F, prob = prob),])

plot(resist_rast)
plot(pts, add = T, pch = 19)

# trans <- gdistance::transition(x = resist_rast,
#                                transitionFunction = function(x)  1 / mean(x),
#                                directions = 8)
# 
# trR <- gdistance::geoCorrection(trans, "r", scl = T)
# resist_mat <- as.matrix(commuteDistance(trR, pts) / 1000) 
# 
# mx <- 0.15
# alpha <- 1/(mx * max(lower(resist_mat)))
# 
# # disp_prob <- exp(-(1 / mean(lower(resist_mat)))*resist_mat)
# disp_prob <- exp(-alpha * resist_mat)
# plot(lower(disp_prob) ~ lower(resist_mat))
# 

dist_rast <- (resist_rast * 0) + 1

cdpop_sim <- cdpop(CDPOP.py = 'C:/Users/peterman.73/OneDrive - The Ohio State University/R/CDPOP/CDPOP-master/src/CDPOP.py',
                   sim_name = sim_name,
                   pts = pts,
                   # resist_rast = disp_prob,
                   resist_rast = resist_rast,
                   sim_dir = "C:/Users/peterman.73/OneDrive - The Ohio State University/R/CDPOP/test_sim",
                   looptime = 101,
                   output_years = 100,
                   gridformat = 'cdpop',
                   loci = 30,
                   alleles = 30,
                   matemoveno = 2, ## 1 = Linear, 5 = Neg exp; 9 = custom prob matrix
                   matemovethresh = '20max',
                   MeanFecundity = 4)

resist_mat <- read.csv("C:/Users/peterman.73/Box/R/CDPOP/test_sim/data/resist_mat.csv",
                       header = F)
hist(lower(resist_mat))

cdpop_grid <- cdpop_sim$grid_list
pops <- cdpop_sim$pop_list[[length(cdpop_sim$pop_list)]]

Dps <- 1-propShared(cdpop_grid[[length(cdpop_grid)]])
pc <- pca_dist(cdpop_grid[[length(cdpop_grid)]], n_axes = 16)
ed <- as.matrix(poppr::edwards.dist(cdpop_grid[[length(cdpop_grid)]]))

plot(lower(Dps) ~ c(dist(pts@coords[pops])))
plot(lower(Dps) ~ lower(resist_mat[pops,pops]))
ecodist::mantel(lower(Dps) ~ lower(resist_mat[pops,pops]))

plot(lower(pc) ~ c(dist(pts@coords[pops])))
plot(lower(pc) ~ lower(resist_mat[pops,pops]))
ecodist::mantel(lower(pc) ~ lower(resist_mat[pops,pops]))

plot(lower(ed) ~ lower(resist_mat[pops,pops]))
ecodist::mantel(lower(ed) ~ lower(resist_mat[pops,pops]))


# Subsample ---------------------------------------------------------------

ind_samp <- sort(sample(1:nInd(cdpop_grid[[length(cdpop_grid)]]), 100))
pops_ <- pops[ind_samp]

plot(lower(Dps[ind_samp,ind_samp]) ~ c(dist(pts@coords[pops_])))
plot(lower(Dps[ind_samp,ind_samp]) ~ lower(resist_mat[pops_,pops_]))
ecodist::mantel(lower(Dps[ind_samp,ind_samp]) ~ lower(resist_mat[pops_,pops_]))

plot(lower(pc[ind_samp,ind_samp]) ~ c(dist(pts@coords[pops_])))
plot(lower(pc[ind_samp,ind_samp]) ~ lower(resist_mat[pops_,pops_]))
ecodist::mantel(lower(pc[ind_samp,ind_samp]) ~  c(dist(pts@coords[pops_])))
ecodist::mantel(lower(pc[ind_samp,ind_samp]) ~ lower(resist_mat[pops_,pops_]))

# 7 = Gaussian function: A * exp ( - (Cost Distance - B)^2 / (2*C^2))
# plot(0.05 * exp(-(lower(resist_mat) - 0.1)^2 / (2*0.1^2)) ~ lower(resist_mat))
plot(1 / lower(resist_mat)^2 ~ lower(resist_mat))

