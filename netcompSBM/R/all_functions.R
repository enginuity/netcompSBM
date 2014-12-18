##@S Currently stores ALL sbm fitting code. This will be split up later. 


#' Run the EM algorithm to fit a SBM
#' 
#' @param adjm Input adjacency matrix. This can contain NAs, values > 1, as a total count over a number of networks.
#' @param Nobs Number of network observations
#' @param Nclass Number of classes in the SBM
#' @param Niter Max number of EM steps
#' @param verbose Lots of output?
#' @param stop_thres Stopping threshold (will stop before max Niter, if change in probabilty estimate matrix is smaller than this value)
#' 
#' @return List of SBM parameters
#' 
#' @export
#' 
fit_SBM = function(adjm, Nobs = 1, Nclass, Niter = 100, verbose = 1, stop_thres = 0.000001) {
  if (FALSE) {
    ## Testing parameters
    adjm = graph; Nclass = 2; Niter = 10; Nobs = 1; verbose = 1; stop_thres = 0.0001
  }
  
  ## Initialize all elements
  N = nrow(adjm)
  nodeps = rep(1/Nclass, length = Nclass)
  edgeps = matrix(sum(adjm, na.rm = TRUE) / (N * (N-1) / 2) * runif(Nclass * Nclass, min = 0.1, max = 0.9), nrow = Nclass)
  edgeps = symmetrize_mat(edgeps)
  
  nodeps_new = nodeps
  edgeps_new = edgeps
  
  H = matrix(0, nrow = N, ncol = Nclass)
  PHI = matrix(runif(Nclass*nrow(adjm), min = 0.1, max = 0.9), nrow = N)
  PHI = PHI / rowSums(PHI)
  
  ## Do iterations
  for(I in 1:Niter) {
    
    ## E step
    H = matrix(0, nrow = N, ncol = Nclass)
    for(r in 1:Nclass) {
      for(s in 1:Nclass) {
        h1 = adjm * log(edgeps[r,s]) + (Nobs - adjm) * log(1-edgeps[r,s])
        H[,r] = H[,r] + sapply(1:N, function(x) {sum(h1[x,-x] * PHI[-x,s], na.rm = TRUE)})
      }
    }
    H = H + abs(max(H))
    peH = t(nodeps * t(exp(H)))
    PHI = peH / rowSums(peH)
    
    ## M step
    nodeps_new = apply(PHI, 2, sum) / N
    for(r in 1:Nclass) {
      for(s in r:Nclass) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num = adjm * Psq
        den = matrix(as.numeric(!is.na(adjm)), nrow = N) * Nobs * Psq
        edgeps_new[r,s] = sum(num[lower.tri(num)], na.rm = TRUE) / sum(den[lower.tri(den)])
      }
    }
    edgeps_new = symmetrize_mat(edgeps_new)
    
    ## Compute change in edge probability matrix
    delta = sum(abs(edgeps_new - edgeps))
    if (verbose > 0) { cat("Iteration ", I, " ----- change in epmat = ", delta, "\n", sep = "") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeprobs = nodeps, edgeprobs = edgeps, classes = sapply(1:N, function(x) {
    order(PHI[x,], decreasing = TRUE)[1]
  })))
}


#' Symmetrizes matrix by using upper-triangular portion
#' 
#' @param mat Input matrix
#' 
#' @return Filled matrix (with lower triangular portion replaced by the upper triangular portion). 
#' 
#' @export
#' 
symmetrize_mat = function(mat) {
  diag_mat = diag(mat)
  mat[lower.tri(mat, diag = TRUE)] = 0
  mat = mat + t(mat)
  diag(mat) = diag_mat
  return(mat)
}


#' Compute the likelihood of a SBM fit
#' 
#' @param A Input adjacency matrix
#' @param fl Output of fit_SBM
#' @param Nobs Number of network observations
#' @param hidden Only return fit on 'hidden' nodes?
#' @param partial_A Partial adjacency matrix (with some nodes NA'd out)
#' 
#' @return Likelihood of the SBM fit
#' 
#' @export
#' 
SBM_likelihood_fit = function(A, fl, Nobs = 1, hidden = FALSE, partial_A) {
  PM = outer(fl$classes, fl$classes, FUN = function(x,y) { mapply(FUN = function(a,b) {fl$edgeprobs[a,b]}, x,y)})
  RM = A * log(PM) + (Nobs - A) * log(1 - PM)
  diag(RM) = 0
  if (hidden) { return(sum(RM[is.na(partial_A)]))}
  return(sum(RM, na.rm = T))
}



#' Compute the likelihood of a SBM fit
#' 
#' @param adjm Input adjacency matrix
#' @param fitl Output of fit_SBM
#' @param Nobs Number of network observations
#' @param mode Cases -- 'full', 'hidden', 'known'
#' -- full gets likelihood of entire input adjm (ignores part_adjm)
#' -- hidden gets only likelihood on hidden nodes (used for CV)
#' -- known gets likelihood on non-hidden nodes (used for fitting on known nodes)
#' @param part_adjm Partial adjacency matrix (with some nodes NA'd out)
#' 
#' @return Likelihood of the SBM fit
#' 
#' @export
#' 
SBM_likelihood_fit_v2 = function(adjm, part_adjm, fitl, mode, Nobs = 1) {
  ## TODO: [Test] this function
  ## Compute edge probability matrix for model 
  PM = outer(fitl$classes, fitl$classes, FUN = function(x,y) { mapply(FUN = function(a,b) {fitl$edgeprobs[a,b]}, x,y)})
  
  RM = adjm * log(PM) + (Nobs - adjm) * log(1 - PM)
  diag(RM) = 0
  
  if (mode == "full") return(sum(RM))
  if (mode == "hidden") return(sum(RM[is.na(partial_A)]))
  if (mode == "known") return(sum(RM[!is.na(partial_A)]))
  stop("Invalid mode")
}


#' Fit a few SBMs, and return the one with best likelihood
#' 
#' @param A Input adjacency matrix
#' @param full_A Full adjacency matrix (if we want to do cv, A would be a matrix with missing values)
#' @param q Number of classes
#' @param Nfits Number of fits to try
#' @param Nobs Number of network observatins
#' @param hidden Use hidden edges to compute likelihood?
#' @param verbose Lots of output?
#' 
#' @return Best model
#' 
#' @export
#' 
search_best_SBM = function(A, full_A = A, q, Nfits, Nobs = 1, hidden = FALSE, verbose = 1) {
  bestlik = -Inf
  bestmod = NULL
  for(j in 1:Nfits) {
    fl = try(fit_SBM(adjm = A, Nclass = q, verbose = verbose - 1, Nobs = Nobs), silent = TRUE)
    if (class(fl) != "try-error") {
      
      newlik = SBM_likelihood_fit(A = full_A, partial_A = A, hidden = hidden, fl = fl)
      if (verbose > 0) {cat("Fit #", j, "-- likelihood = " , newlik, "\n")}
      if (newlik > bestlik) { 
        bestmod = fl 
        bestlik = newlik
      }
    } else {
      if (verbose > 0) {cat("Fit #", j, "-- FAILED" ,"\n")}
    }
  }
  return(bestmod)
}


#' Fit a few SBMs, and return the one with best likelihood
#' 
#' @param adjm Input adjacency matrix
#' @param full_adjm If null, adjm is the `full' adjacency matrix. If non-null, this is a CV iteration, and full(er) adjacency matrix is input here.  
#' @param Nclass Number of classes
#' @param Nfits Number of fits to try
#' @param Nobs Number of network observatins
#' @param do_cv Do cross-validation?
#' @param verbose Lots of output?
#' 
#' @return Best model
#' 
#' @export
#' 
search_best_SBM_v2 = function(adjm, full_adjm = NULL, Nclass, Nfits = 50, Nobs = 1, do_cv = FALSE, verbose = 1, Niter = 100, stop_thres = 0.000001) {
  bestlik = -Inf
  bestmod = NULL
  
  for (F in 1:Nfits) {
    fl = try(fit_SBM(adjm = adjm, Nobs = Nobs, Nclass = Nclass, Niter = Niter, verbose = verbose - 1, stop_thres = stop_thres))
    if (class(fl) != "try-error") {
      if (do_cv) {
        newlik = SBM_likelihood_fit_v2(adjm = full_adjm, part_adjm = adjm, fitl = fl, mode = "known", Nobs = Nobs)
      } else {
        newlik = SBM_likelihood_fit_v2(adjm = full_adjm, fitl = fl, mode = "full", Nobs = Nobs)
      }
      
      if (verbose > 0) { cat("Fit #", j, "-- likelihood = " , newlik, "\n") }
      
      ## Update best model if new model is better 
      if (newlik > bestlik) { 
        bestmod = fl
        bestlik = newlik
      }
    } else {
      if (verbose > 0) { cat("Fit #", j, "-- FAILED" ,"\n") }
    }
  }
  return(bestmod)
}



#' Hides a random set of edges
#' 
#' @param adjm Input adjacency matrix
#' @param frac Fraction of edges to hide
#' 
#' @return Adjacency matrix with hidden edges
#' 
#' @export
#' 
hide_edges = function(adjm, frac = 0.1) {
  tm = matrix(1:(nrow(adjm)^2), nrow = nrow(adjm))
  vals = tm[upper.tri(tm)]
  vals = sample(vals, size = floor(length(vals) * frac))
  adjm[vals] = NA
  return(symmetrize_mat(adjm))
}


#' Run CV for the SBM
#' 
#' @param A Input adjacency matrix
#' @param qs Vector of class sizes
#' @param Nfits Number of fits to find 'optimal'
#' @param Nobs Number of network observations
#' @param CV_folds Number of CV runs per number of classes
#' @param verbose Lots of output?
#' 
#' @return CV error (likelihoods for link prediction)
#' 
#' @export
#' 
CV_SBM = function(A, qs, Nfits = 50, Nobs, CV_folds = 10, verbose = 1) {
  liks_mat = matrix(0, nrow = CV_folds, ncol = length(qs))
  for(j in 1:CV_folds) {
    if (verbose > 0) {cat("\nCV fold number", j, date(), "\n")}
    subA = hide_edges(A)
#|----##--Reparameterizing this function --Thu Dec 18 00:37:59 2014--
    for(Q in seq_along(qs)) {
      if (verbose > 0) { cat(qs[Q], "-") }
      bestfit = search_best_SBM(A = subA, full_A = A, q = qs[Q], Nfits = Nfits, Nobs = Nobs, verbose = verbose - 1)
      liks_mat[j,Q] = SBM_likelihood_fit(A = A, fl = bestfit, hidden = TRUE, partial_A = subA)
    }
  }
  return(liks_mat)
}


