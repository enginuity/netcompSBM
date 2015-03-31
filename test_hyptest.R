
## Load libraries
library(netcompLib)
library(netcompSBM)
library(diagram)
library(abind)
library(igraph)
load("../../network-comparison/netcomp-project/data/method_data/small_samp_DFcorr.Rdata")



# Test code for inverting hidden edges ------------------------------------

invert_hidden_edges = function(orig_adjm, hidden_adjm) {
  res = orig_adjm
  res[!is.na(hidden_adjm)] = NA
  diag(res) = 0
  return(res)
}

gsize = 6
netm = NetworkModelSBM(Nnodes = gsize, model_param = set_model_param(block_assign = rep(1:2, each = gsize/2), block_probs = matrix(c(.7, .2, .2, .4), nrow = 2)))
adjm = sampleNetwork(netm)[,,1]

adjm_hid = hide_edges(adjm, frac = 0.4)
invert_hidden_edges(adjm, adjm_hid)


#  ------------------------------------------------------------------------


find_best_model = function(adjm, Ntries = 20) {
  bestLL = -2394792
  bestmod = NULL
  for(j in 1:Ntries) {
    cat(".")
    
    test = try(fit_SBM_v2(adjm = adjm_hid, Nobs = 2, Nclass = 3, Niter = 100, method = "mf", verbose = 0), silent = TRUE)
    if (class(test) != "try-error") {
      loglik = compute_sbm_loglik(test$classes, adjm_hid, test$nodeps, test$edgeps, Nobs = 2)
      if (loglik > bestLL) {
        bestmod = test
        bestLL = loglik
      }
    }
  }
  return(bestmod)
}

renumber_classes = function(v) {
  # if vector is c(1,1,2,4) -> c(1,1,2,3)
  if (max(v) == length(unique(v))) { return(v) }
  while(max(v) > length(unique(v))) {
    m = min(setdiff(1:max(v), v))
    v[v == max(v)] = m
  }
  return(v)
}

gsize = 60
netm = NetworkModelSBM(Nnodes = gsize, model_param = set_model_param(block_assign = rep(1:3, each = gsize/3), block_probs = matrix(c(.5, .2, .1, .2, .3, .6, .1, .6, .7), nrow = 3)))


pval_sim = rep(0, times = 200) 

for(j in 1:200) {
  print(j)
  adjm1 = sampleNetwork(netm)[,,1]
  adjm2 = sampleNetwork(netm)[,,1]
  adjm = adjm1 + adjm2
  adjm_hid = hide_edges(adjm, frac = 0.5)
  
  res = find_best_model(adjm_hid)
  res$classes = renumber_classes(res$classes)
  adjm_test = invert_hidden_edges(orig_adjm = adjm, adjm_hid)
  
  pval_sim[j] = computePval(NetS = NetworkStructSBM(Nnodes = 60, model_param = set_model_param(block_assign = res$classes, block_nclass = length(unique(res$classes)))), adja1 = invert_hidden_edges(adjm1, adjm_hid), adja2 = invert_hidden_edges(adjm2, adjm_hid), Nobs = 1, pl = set_sim_param(n_models = 1, cc_adj = 2, thres_ignore = 5), mode = "default")
}

hist(pval_sim, freq = TRUE, breaks = 20)

hist(runif(200), breaks = 20)
#  ------------------------------------------------------------------------
# 
# a1 = t(adjm1) %*% adjm1
# eigen(a1)
# round(eigen(adjm1)$values)
# plot(1:60, round(eigen(adjm1)$vectors[,1], 2))
# eig
# 
# #test = fit_SBM_v2(adjm = adjm, Nclass = 2, Niter = 100, method = "mf", verbose = 2)
