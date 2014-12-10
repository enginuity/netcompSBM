
fit_SBM = function(A, Nobs = 1, q, niter = 100, verbose = TRUE, stop_thres = 0.000001) {
  ## A = adjacency matrix ( can contain NAs, can contain double (multiple-count) edges)
  ## Nobs = # of added adjacency matrices in (A)
  if (FALSE) {
    A = graph; q = 2; niter = 10; Nobs = 1; verbose = TRUE; stop_thres = 0.0001
  }
  ## Initialize all elements
  N = nrow(A)
  nodeps = rep(1/q, length = q)
  edgeps = matrix(sum(A, na.rm = TRUE) / (N * (N-1) / 2) * runif(q * q, min = 0.1, max = 0.9), nrow = q)
  newedgep = edgeps
  edgeps = symmetrize_mat(edgeps)
  
  H = matrix(0, nrow = N, ncol = q)
  PHI = matrix(runif(q*nrow(A), min = 0.1, max = 0.9), nrow = N)
  PHI = PHI / rowSums(PHI)
  
  
  ## Do iterations
  for(I in 1:niter) {
    
    ## E step
    
    H = matrix(0, nrow = N, ncol = q)
    for(r in 1:q) {
      for(s in 1:q) {
        h1 = A * log(edgeps[r,s]) + (Nobs - A) * log(1-edgeps[r,s])
        H[,r] = H[,r] + sapply(1:N, function(x) {sum(h1[x,-x] * PHI[-x,s], na.rm = TRUE)})
      }
    }
    H = H + abs(max(H))
    peH = t(nodeps * t(exp(H)))
    PHI = peH / rowSums(peH)
    
    ## M step
    nodeps = apply(PHI, 2, sum) / N
    for(r in 1:q) {
      for(s in r:q) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num = A * Psq
        den = matrix(as.numeric(!is.na(A)), nrow = N) * Nobs * Psq
        newedgep[r,s] = sum(num[lower.tri(num)], na.rm = TRUE) / sum(den[lower.tri(den)])
      }
    }
    newedgep = symmetrize_mat(newedgep)
    
    ## Compute change in edge probability matrix
    delta = sum(abs(newedgep - edgeps))
    edgeps = newedgep
    if (verbose) { cat("Iteration ", I, " ----- change in epmat = ", delta, "\n", sep = "") }
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeprobs = nodeps, edgeprobs = edgeps, classes = sapply(1:N, function(x) {
    order(PHI[x,], decreasing = TRUE)[1]
  })))
}

symmetrize_mat = function(m) {
  dm = diag(m)
  m[lower.tri(m, diag = TRUE)] = 0
  m = m + t(m)
  diag(m) = dm
  return(m)
  
}


SBM_likelihood_fit = function(A, fl, Nobs = 1, hidden = FALSE, partial_A) {
  ## fl = fit list (output of fit_SBM)
  PM = outer(fl$classes, fl$classes, FUN = function(x,y) { mapply(FUN = function(a,b) {fl$edgeprobs[a,b]}, x,y)})
  RM = A * log(PM) + (Nobs - A) * log(1 - PM)
  diag(RM) = 0
  if (hidden) { return(sum(RM[is.na(partial_A)]))}
  return(sum(RM, na.rm = T))
}

search_best_SBM = function(A, full_A = A, q, Nfits, Nobs = 1, hidden = FALSE, verbose = TRUE) {
  bestlik = -Inf
  bestmod = NULL
  for(j in 1:Nfits) {
    fl = try(fit_SBM(A = A, q = q, verbose = FALSE, Nobs = Nobs), silent = TRUE)
    if (class(fl) != "try-error") {
      
      newlik = SBM_likelihood_fit(A = full_A, partial_A = A, hidden = hidden, fl = fl)
      if (verbose) {cat("Fit #", j, "-- likelihood = " , newlik, "\n")}
      if (newlik > bestlik) { 
        bestmod = fl 
        bestlik = newlik
      }
    } else {
      if (verbose) {cat("Fit #", j, "-- FAILED" ,"\n")}
    }
  }
  return(bestmod)
}

hide_edges = function(A, frac = 0.1) {
  tm = matrix(1:(nrow(A)^2), nrow = nrow(A))
  vals = tm[upper.tri(tm)]
  vals = sample(vals, size = floor(length(vals) * frac))
  A[vals] = NA
  return(symmetrize_mat(A))
}

CV_SBM = function(A, qs, Nfits = 50, Nobs, CV_folds = 10, verbose = TRUE) {
  liks_mat = matrix(0, nrow = CV_folds, ncol = length(qs))
  for(j in 1:CV_folds) {
    if (verbose) {cat("\nCV fold number", j, date(), "\n")}
    subA = hide_edges(A)
    for(Q in seq_along(qs)) {
      if (verbose) { cat(qs[Q], "-") }
      bestfit = search_best_SBM(A = subA, full_A = A, q = qs[Q], Nfits = Nfits, Nobs = Nobs, verbose = FALSE)
      liks_mat[j,Q] = SBM_likelihood_fit(A = A, fl = bestfit, hidden = TRUE, partial_A = subA)
    }
  }
  return(liks_mat)
}









