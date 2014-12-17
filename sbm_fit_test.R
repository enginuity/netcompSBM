source("setup.R")
source("src_r/sbm_fit_functions.R")


gsize = 100
g = sample_generating_models(Nnodes = gsize, mode = "block")
g$m1$assign = rep(1:2, each = gsize/2)
g$m1$probmat = matrix(c(.7, .3, .3, .5), nrow = 2)
graph = sample_network_pair(g)[[1]][,,1]




fl = fit_SBM(A = graph, q = 2)
#|----##--Reparameterizing this function --Wed Dec 17 15:10:23 2014--



test = search_best_SBM(A = graph, q = 2, Nfits = 50)

gsize = 100
g = sample_generating_models(Nnodes = gsize, mode = "block")
g$m1$assign = rep(1:2, each = gsize/2)
g$m1$probmat = matrix(c(.7, .3, .3, .5), nrow = 2)

graph2 = apply(replicate(sample_network_pair(g)[[1]][,,1], n = 10), c(1,2), sum)
test = search_best_SBM(A = graph2, q = 2, Nfits = 50, Nobs = 10)
test = fit_SBM(A = graph2, q = 2, Nobs = 10)
#|----##--Reparameterizing this function --Wed Dec 17 15:10:23 2014--




graph3 = hide_edges(graph)
test3 = search_best_SBM(A = graph3, full_A = graph, q = 2, Nfits = 50)
SBM_likelihood_fit(A = graph, fl = test3, hidden = TRUE, partial_A = graph3)





CV_Q = rep(0, times = 5)
for(Q in 1:5) {
  print("-----------------------------")
  print(Q)
  print("-----------------------------")
  test4 = CV_SBM(A = graph, q = Q, Nobs = 1)
  CV_Q[Q] = mean(test4)
}

plot(1:5, CV_Q, xlab = "Number of classes", ylab = "Log-lik on test data", type = "b", main = "CV on 2-class data")

gsize = 100
g = sample_generating_models(Nnodes = gsize, mode = "block")
g$m1$assign = c(rep(1:2, each = 30), rep(3, each = 40))
g$m1$probmat = matrix(c(.2, .1, .13, .1, .4, .28, .13, .28, .3), nrow = 3)
graph4 = sample_network_pair(g)[[1]][,,1]




## Comparing differences between two fits from same network
g1 = sample_network_pair(g)[[1]][,,1]
g2 = sample_network_pair(g)[[1]][,,1]

best1 = search_best_SBM(A = g1, q = 3, Nobs = 1, Nfits = 50)
best2 = search_best_SBM(A = g2, q = 3, Nobs = 1, Nfits = 50)

find_KL_dist = function(fit1, fit2) {
  ## Finds nodal KL distance
  KL = function(p, q) {return(p * log(p / q) + (1 - p) * log( (1 - p) / (1 - q))) }
  avKL = function(p,q) { 0.5 * sum(KL(p,q) + KL(q,p)) }
  N = length(fit1$classes)
  sapply(1:N, function(x) {avKL(
    fit1$edgeprobs[fit1$classes[x],fit1$classes],
    fit2$edgeprobs[fit2$classes[x],fit2$classes]
  )}) -> res
  return(res)
}

cbind(1:100, find_KL_dist(best1, best2)) ->  KLdists
colnames(KLdists) = c("Node", "KL distance")
KLdists[order(KLdists[,2],decreasing = TRUE),]

library(clusteval)
cluster_similarity(best1$classes, best2$classes)


gsize = 100
g = sample_generating_models(Nnodes = gsize, mode = "block")
g$m1$assign = c(rep(1:2, each = 30), rep(3, each = 40))
g$m1$probmat = matrix(c(.2, .1, .13, .1, .4, .78, .13, .28, .3), nrow = 3) * .1
graph = sample_network_pair(g)[[1]][,,1]

test = CV_SBM(A = graph, qs = 2:5, Nfits = 20, Nobs = 1, CV_folds = 10, verbose = TRUE)


## TRy combining matrices

comb_mat_diag = function(m1, m2) {
  res = array(NA, dim = dim(m1) + dim(m2))
  res[1:nrow(m1), 1:ncol(m1)] = m1
  res[-1:-nrow(m1),-1:-ncol(m1)] = m2
  return(res)
}

gsize = 50
g = sample_generating_models(Nnodes = gsize, mode = "block")
g$m1$assign = c(rep(1:2, each = 25))
g$m1$probmat = matrix(c(.2, .4, .4, .7), nrow = 2)
graph = comb_mat_diag(sample_network_pair(g)[[1]][,,1], sample_network_pair(g)[[1]][,,1])

test1 = search_best_SBM(A = graph, q = 3, Nfits = 50)




