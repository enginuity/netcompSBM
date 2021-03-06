
## Load libraries
library(netcompLib)
library(netcompSBM)
library(diagram)
library(abind)
library(igraph)


gsize = 60
netm = NetworkModelSBM(Nnodes = gsize, model_param = set_model_param(block_assign = rep(1:2, each = gsize/2), block_probs = matrix(c(.7, .2, .2, .4), nrow = 2)))
#|----##Function parameters changed -- only model_params --Thu Jul 30 20:22:23 2015--
adjm = sampleNetwork(netm)[,,1]

test = fit_SBM_v2(adjm = adjm, Nclass = 2, Niter = 100, method = "mf", verbose = 2)
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
test = fit_SBM_v2(adjm = adjm, Nclass = 2, Niter = 50, method = "bp", verbose = 2)
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
test = fit_SBM_v2(adjm = adjm, Nclass = 2, Niter = 100, method = "bp_scaled", verbose = 2)
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--


test = fit_SBM_BP(adjm = adjm, Nclass = 2, Niter = 100, stop_thres = 0.001)

load("test.RData")

diag(test.mat) = 0
fit_SBM_v2(adjm = test.mat, Nclass = 3, Niter = 100, method = "mf", verbose = 2)
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
