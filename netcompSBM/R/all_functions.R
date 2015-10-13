##@S Currently stores ALL sbm fitting code. This will be split up later. 


#' Compute the likelihood of a SBM fit
#' 
#' @param A Input adjacency matrix
#' @param fl Output of fit_SBM
#|----##**fit_SBM will be replaced by fit_SBM_v2; update this code to correspond with the new code... --Tue May 26 21:00:00 2015--
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
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



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (SBM_likelihood_fit_v2)
#' Compute the likelihood of a SBM fit
#' 
#' @param adjm Input adjacency matrix
#' @param part_adjm Partial adjacency matrix (with some nodes NA'd out)
#' @param fitl Output of fit_SBM
#|----##**fit_SBM will be replaced by fit_SBM_v2; update this code to correspond with the new code... --Tue May 26 21:00:00 2015--
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
#' @param mode Cases -- 'full', 'hidden', 'known'
#' -- full gets likelihood of entire input adjm (ignores part_adjm)
#' -- hidden gets only likelihood on hidden nodes (used for CV)
#' -- known gets likelihood on non-hidden nodes (used for fitting on known nodes)
#' @param Nobs Number of network observations
#' 
#' @return Likelihood of the SBM fit
#' 
#' @export
#' 
SBM_likelihood_fit_v2 = function(adjm, part_adjm, fitl, mode, Nobs = 1) {
  ## TODO: [Test] this function
  ## Compute edge probability matrix for model 
  PM = outer(fitl$classes, fitl$classes, FUN = function(x,y) { mapply(FUN = function(a,b) {fitl$edgeps[a,b]}, x,y)})
  ## NOTE CHANGED HERE -- wont work for old code anymore
  
  RM = adjm * log(PM) + (Nobs - adjm) * log(1 - PM)
  diag(RM) = 0
  
  if (mode == "full") return(sum(RM))
  if (mode == "hidden") return(sum(RM[is.na(part_adjm)]))
  if (mode == "known") return(sum(RM[!is.na(part_adjm)]))
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
#|----##**fit_SBM will be replaced by fit_SBM_v2; update this code to correspond with the new code... --Tue May 26 21:00:00 2015--
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
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


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (search_best_SBM_v2)
#' Fit a few SBMs, and return the one with best likelihood
#' 
#' @param adjm Input adjacency matrix
#' @param full_adjm If null, adjm is the `full' adjacency matrix. If non-null, this is a CV iteration, and full(er) adjacency matrix is input here.  
#' @param Nclass Number of classes
#' @param Nfits Number of fits to try
#' @param Nobs Number of network observatins
#' @param do_cv Do cross-validation?
#' @param verbose Lots of output?
#' @param Niter temp
#' @param stop_thres temp
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
#|----##**fit_SBM will be replaced by fit_SBM_v2; update this code to correspond with the new code... --Tue May 26 21:00:00 2015--
#|----##Renamed to fit_SBM --Tue May 26 21:03:28 2015--
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
    for(Q in seq_along(qs)) {
      if (verbose > 0) { cat(qs[Q], "-") }
      bestfit = search_best_SBM(A = subA, full_A = A, q = qs[Q], Nfits = Nfits, Nobs = Nobs, verbose = verbose - 1)
      liks_mat[j,Q] = SBM_likelihood_fit(A = A, fl = bestfit, hidden = TRUE, partial_A = subA)
    }
  }
  return(liks_mat)
}




## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (CV_SBM_v2)
#' Run CV for the SBM
#' 
#' @param adjm Input adjacency matrix
#' @param Nclass_v Vector of class sizes
#' @param Nfits Number of fits to find 'optimal'
#' @param Nobs Number of network observations
#' @param NCV_folds Number of CV runs per number of classes
#' @param verbose Lots of output?
#' @param Niter temp
#' @param stop_thres temp
#' 
#' @return CV error (likelihoods for link prediction)
#' 
#' @export
#' 
CV_SBM_v2 = function(adjm, Nclass_v, Nfits = 50, Nobs = 1, NCV_folds = 10, verbose = 1, Niter = 100, stop_thres = 0.000001) {
  liks_mat = matrix(0, nrow = NCV_folds, ncol = length(Nclass_v))
  for(j in 1:NCV_folds) {
    if (verbose > 0) {cat("\nCV fold number", j, date(), "\n")}
    sub_adjm = hide_edges(adjm)
    for(Q in seq_along(Nclass_v)) {
      if (verbose > 0) { cat(Nclass_v[Q], "-") }
      bestfit = search_best_SBM_v2(adjm = sub_adjm, full_adjm = adjm, Nclass = Nclass_v[Q], Nfits = Nfits, Nobs = Nobs, do_cv = TRUE, verbose = verbose - 1, Niter = Niter, stop_thres = stop_thres)
      liks_mat[j,Q] = SBM_likelihood_fit_v2(adjm = adjm, part_adjm = sub_adjm, fitl = bestfit, mode = "hidden", Nobs = Nobs)
    }
  }
  return(liks_mat)
}


