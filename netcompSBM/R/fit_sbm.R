##@S Currently stores ALL sbm fitting code. This will be split up later. 

## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (fit_SBM_v2)
#' Run the EM algorithm to fit a SBM
#' 
#' @param adjm Input adjacency matrix. This can contain NAs, values > 1, as a total count over a number of networks.
#' @param Nobs Number of network observations
#' @param Nclass Number of classes in the SBM
#' @param Niter Max number of EM steps
#' @param verbose Lots of output?
#' @param stop_thres Stopping threshold (will stop before max Niter, if change in probabilty estimate matrix is smaller than this value)
#' @param method allows 'bp', 'mf', 'bp_scaled'
#' 
#' @return List of SBM parameters
#' 
#' @export
#' 
fit_SBM_v2 = function(adjm, Nobs = 1, Nclass, Niter = 100, verbose = 1, stop_thres = 0.000001, method = "bp") {
  require(igraph)
  
  if (FALSE) {
    ## Testing parameters
    adjm = graph; Nclass = 2; Niter = 10; Nobs = 1; verbose = 1; stop_thres = 0.0001
  }
  
  ## Initialize variables (with random starts)
  ## TODO: [Upgrade] Allow input of randomness (or of starting points)
  N = nrow(adjm) ## This is the number of nodes
  nodeps = rep(1/Nclass, length = Nclass)
  edgeps = symmetrize_mat(matrix(sum(adjm, na.rm = TRUE) / (Nobs * N * (N-1) / 2) * runif(n = Nclass * Nclass, min = 0.1, max = 0.9), nrow = Nclass))
  
  if (method == "mf") { 
    H = matrix(0, nrow = N, ncol = Nclass)
    PHI = matrix(runif(Nclass*N, min = 0.1, max = 0.9), nrow = N)
    PHI = PHI / rowSums(PHI)
    
    results = EM_SBM_mf(adjm = adjm, Nobs = Nobs, nodeps = nodeps, edgeps = edgeps, H = H, PHI = PHI, 
                        Niter = Niter, stop_thres = stop_thres, verbose = verbose)
  }
  if (method == "bp") {
    PHI = matrix(0, nrow = Nclass, ncol = N)
    AUX = rep(0, times = Nclass)  
    
    results = EM_SBM_bp(adjm = adjm, Nobs = Nobs, nodeps = nodeps, edgeps = edgeps, PHI = PHI, AUX = AUX, 
                        Niter = Niter, stop_thres = stop_thres, verbose = verbose)
  }
  if (method == "bp_rescaled") {
    PHI = matrix(0, nrow = Nclass, ncol = N)
    AUX = rep(0, times = Nclass)  
    
    results = EM_SBM_bp_rescaled(adjm = adjm, Nobs = Nobs, nodeps = nodeps, edgeps = edgeps, PHI = PHI, AUX = AUX, 
                                 Niter = Niter, stop_thres = stop_thres, verbose = verbose)
  }
  
  return(results)
}




## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_mf)
#' <What does this function do>
#' 
#' @param adjm temp
#' @param Nobs temp
#' @param nodeps temp
#' @param edgeps temp
#' @param H temp
#' @param PHI temp
#' @param Niter temp
#' @param stop_thres temp
#' @param verbose temp
#' 
#' @return temp
#' 
#' @export
#' 
EM_SBM_mf = function(adjm, Nobs, nodeps, edgeps, H, PHI, Niter, stop_thres, verbose) {
  
  N = nrow(adjm)
  Nclass = nrow(edgeps)
  
  for(I in 1:Niter) {
    ## E step
    H = H * 0
    for(r in 1:Nclass) { for(s in 1:Nclass) {
      h1 = adjm * log(edgeps[r,s]) + (Nobs - adjm) * log(1 - edgeps[r,s])
      H[,r] = H[,r] + sapply(1:N, function(x) {sum(h1[x,-x] * PHI[-x,s], na.rm = TRUE)})
    }}
    H = H + abs(max(H)) #rescale to prevent exponentiation errors
    peH = t(nodeps * t(exp(H)))
    PHI = peH / rowSums(peH)
    
    ## M step
    nodeps_new = apply(PHI, 2, sum) / N
    edgeps_new = edgeps * 0
    for(r in 1:Nclass) {
      for(s in r:Nclass) {
        Psq = PHI[,r,drop = FALSE] %*% t(PHI[,s,drop = FALSE])
        num = adjm * Psq
        den = matrix(as.numeric(!is.na(adjm)), nrow = N) * Nobs * Psq
        edgeps_new[r,s] = max(min(sum(num[lower.tri(num)], na.rm = TRUE) / sum(den[lower.tri(den)]), 0.999), 0.001) ## bound the prob by 0.999 and 0.001 to prevent weird stuff?
      }
    }
    edgeps_new = symmetrize_mat(edgeps_new)    
    
    ## Check for convergence
    ## Compute change in parameter estimates
    delta = sum(abs(edgeps_new - edgeps)) + sum(abs(nodeps_new - nodeps))
    if (verbose > 0) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    if (delta < stop_thres) { break }
  }
  
  return(list(nodeps = nodeps, edgeps = edgeps, PHI = PHI, 
              classes = sapply(1:N, function(x) { order(PHI[x,], decreasing = TRUE)[1] }), nsteps = I))
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_bp)
#' <What does this function do>
#' 
#' @param adjm temp
#' @param Nobs temp
#' @param nodeps temp
#' @param edgeps temp
#' @param PHI temp
#' @param AUX temp
#' @param Niter temp
#' @param stop_thres temp
#' @param verbose temp
#' 
#' @return temp
#' 
#' @export
#' 
EM_SBM_bp = function(adjm, Nobs, nodeps, edgeps, PHI, AUX, Niter, stop_thres, verbose) {
  ## Helper functions for BP
  compute_aggmsg = function(id = NULL) {
    # assumes MSG_L and PHI exists
    # for node i: 
    if (is.null(id)) { inds = 1:N } else { inds = id }
    un_norm = lapply(inds, function(i) {
      if (length(EL[[i]] > 0)) {
        return(nodeps * exp(-AUX) * sapply(lapply(1:Nclass, function(x) {apply(cab[x,] * MSG_L[[i]], 2, sum)}), prod)) ## equation 26
      } else { nodeps * exp(-AUX) }
    })
    return(do.call(what = cbind, args = lapply(un_norm, function(x) { x / sum(x)})))
  }
  compute_aux = function() {
    ## assumes AUX & all other stuff exists; does global assignmet. 
    res = sapply(1:Nclass, function(x){ sum(edgeps[x,] * PHI) } ) ## eqn 27
    return(res - min(res)) ## make this simpler to exponentiate
  }
  randomize_msg_order = function() {
    # generates a random message order
    all_msgs = do.call(rbind, lapply(1:N, function(x) { if (length(EL[[x]] > 0)) {cbind(x, seq_along(EL[[x]]))} else { NULL }}))
    all_msgs = all_msgs[sample(1:nrow(all_msgs)),]
    return(all_msgs)
  }
  extract_msg_pair = function(ii,jj) {
    mij = MSG_L[[jj]][,which(EL[[jj]] == ii)]
    mji = MSG_L[[ii]][,which(EL[[ii]] == jj)]
    contrib_mat = t(t(mij)) %*% t(mji) * cab
    return(list(Z = sum(contrib_mat), contrib_mat = contrib_mat))
  }
  find_edgelist = function(EL) { 
    temp = sapply(1:N, function(x) { EL[[x]][EL[[x]] > x]})
    edgesmat = do.call(rbind, sapply(1:N, function(x) { if (length(temp[[x]]) > 0) cbind(x, temp[[x]])}))
    return(edgesmat)
  }
  
  N = nrow(adjm)
  Nclass = nrow(edgeps)
  G = graph.adjacency(adjm, mode = "undirected")
  EL = get.adjlist(G)
  
  max_estep_tries = max(floor(sqrt(N) / 2), 5)
  cab = edgeps * N
  
  ## Do iterations
  for(I in 1:Niter) {
    # Per step BP initialization
    MSG_L = lapply(EL, FUN = function(x) { ## stores the in messages. 
      if(length(x) > 0) {
        z = matrix(runif(length(x) * Nclass), nrow = Nclass) 
        z = sapply(1:ncol(z), function(y) {z[,y] / sum(z[,y])})
        return(z)
      } else { return(NULL) }
    })
    PHI = compute_aggmsg()
    AUX = compute_aux()
    print(AUX) 
    
    # E step
    tries = 0
    conv = stop_thres + 10
    while((conv > stop_thres) & (tries < max_estep_tries)) {
      conv = 0; tries = tries + 1
      msg_order = randomize_msg_order()
      for(j in 1:nrow(msg_order)) {
        msg_start = EL[[msg_order[[j,1]]]][msg_order[j,2]]
        msg_end = msg_order[j,1]
        rem_nbr_inds = which(!(EL[[msg_start]] %in% msg_end))
        
        new_msgl = nodeps * exp(-AUX)
        if (length(rem_nbr_inds) > 0) {
          new_msgl = new_msgl * sapply(lapply(1:Nclass, function(x) {apply(cab[x,] * MSG_L[[msg_start]][,rem_nbr_inds, drop = FALSE], 2, sum)}), prod)
        } # else multiply by 1; eg do nothing
        new_msgl = new_msgl / sum(new_msgl)
        conv = conv + sum(abs(new_msgl - MSG_L[[msg_end]][,msg_order[j,2]]))
        
        ## UPDATE PHI
        MSG_L[[msg_end]][,msg_order[j,2]] = new_msgl
        old_msg = PHI[,msg_end]
        PHI[,msg_end] = compute_aggmsg(id = msg_end)
        
        ## UPDATE AUX FIELD
        old_contrib = sapply(1:Nclass, function(x) {sum(1/N * cab[x,] * old_msg)})
        new_contrib = sapply(1:Nclass, function(x) {sum(1/N * cab[x,] * PHI[,msg_end])})
        AUX = AUX - old_contrib + new_contrib
      }
    }
    
    
    # M step
    nodeps_new = apply(PHI, 1, sum) / N
    
    ## Compute normalization factor
    edgelist = find_edgelist(EL)
    edgecont = lapply(seq_len(nrow(edgelist)), function(x) { extract_msg_pair(edgelist[x,1], edgelist[x,2])})
    
    cab_new = matrix(0, nrow = Nclass, ncol = Nclass)
    for(m in seq_len(nrow(edgelist))) {
      cab_new = cab_new + (edgecont[[m]]$contrib_mat / edgecont[[m]]$Z)
    }
    cab_new = cab_new + t(cab_new)
    cab_new = cab_new / crossprod(t(nodeps_new)) / N
    edgeps_new = cab_new / N    
    
    
    ## Compute change in parameter estimates
    delta = sum(abs(edgeps_new - edgeps)) + sum(abs(nodeps_new - nodeps))
    if (verbose > 0) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[,x], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    cab = cab_new
    
    ## Stop if threshold is met
    if (delta < stop_thres) { break }
  }
  return(list(nodeps = nodeps, edgeps = edgeps, PHI = t(PHI), 
              classes = sapply(1:N, function(x) { order(PHI[,x], decreasing = TRUE)[1] }), nsteps = I))
}


## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (EM_SBM_bp_rescaled)
#' <What does this function do>
#' 
#' @param adjm temp
#' @param Nobs temp
#' @param nodeps temp
#' @param edgeps temp
#' @param PHI temp
#' @param AUX temp
#' @param Niter temp
#' @param stop_thres temp
#' @param verbose temp
#' 
#' @return temp
#' 
#' @export
#' 
EM_SBM_bp_rescaled = function(adjm, Nobs, nodeps, edgeps, PHI, AUX, Niter, stop_thres, verbose) {
  ## This version uses the probabilities intead of the cab matrix. 
  
  ## Helper functions for BP
  compute_aggmsg = function(id = NULL) {
    # assumes MSG_L and PHI exists
    # for node i: 
    if (is.null(id)) { inds = 1:N } else { inds = id }
    un_norm = lapply(inds, function(i) {
      if (length(EL[[i]] > 0)) {
        return(nodeps * exp(-AUX) * sapply(lapply(1:Nclass, function(x) {apply(edgeps[x,] * MSG_L[[i]], 2, sum)}), prod)) ## equation 26
      } else { nodeps * exp(-AUX) }
    })
    return(do.call(what = cbind, args = lapply(un_norm, function(x) { x / sum(x)})))
  }
  compute_aux = function() {
    ## assumes AUX & all other stuff exists; does global assignmet. 
    res = sapply(1:Nclass, function(x){ sum(edgeps[x,] * PHI) } ) ## eqn 27
    return(res - min(res)) ## make this simpler to exponentiate
  }
  randomize_msg_order = function() {
    # generates a random message order
    all_msgs = do.call(rbind, lapply(1:N, function(x) { if (length(EL[[x]] > 0)) {cbind(x, seq_along(EL[[x]]))} else { NULL }}))
    all_msgs = all_msgs[sample(1:nrow(all_msgs)),]
    return(all_msgs)
  }
  extract_msg_pair = function(ii,jj) {
    mij = MSG_L[[jj]][,which(EL[[jj]] == ii)]
    mji = MSG_L[[ii]][,which(EL[[ii]] == jj)]
    contrib_mat = t(t(mij)) %*% t(mji) * edgeps
    return(list(Z = sum(contrib_mat), contrib_mat = contrib_mat))
  }
  find_edgelist = function(EL) { 
    temp = sapply(1:N, function(x) { EL[[x]][EL[[x]] > x]})
    edgesmat = do.call(rbind, sapply(1:N, function(x) { if (length(temp[[x]]) > 0) cbind(x, temp[[x]])}))
    return(edgesmat)
  }
  
  N = nrow(adjm)
  Nclass = nrow(edgeps)
  G = graph.adjacency(adjm, mode = "undirected")
  EL = get.adjlist(G)
  
  max_estep_tries = max(floor(sqrt(N) / 2), 5)
  
  ## Do iterations
  for(I in 1:Niter) {
    # Per step BP initialization
    MSG_L = lapply(EL, FUN = function(x) { ## stores the in messages. 
      if(length(x) > 0) {
        z = matrix(runif(length(x) * Nclass), nrow = Nclass) 
        z = sapply(1:ncol(z), function(y) {z[,y] / sum(z[,y])})
        return(z)
      } else { return(NULL) }
    })
    PHI = compute_aggmsg()
    AUX = compute_aux()
    print(AUX)
    # E step
    tries = 0
    conv = stop_thres + 10
    while((conv > stop_thres) & (tries < max_estep_tries)) {
      if (verbose > 5) { print(MSG_L[[1]]) }
      conv = 0; tries = tries + 1
      msg_order = randomize_msg_order()
      for(j in 1:nrow(msg_order)) {
        msg_start = EL[[msg_order[[j,1]]]][msg_order[j,2]]
        msg_end = msg_order[j,1]
        rem_nbr_inds = which(!(EL[[msg_start]] %in% msg_end))
        
        new_msgl = nodeps * exp(-AUX)
        if (length(rem_nbr_inds) > 0) {
          new_msgl = new_msgl * sapply(lapply(1:Nclass, function(x) {apply(edgeps[x,] * MSG_L[[msg_start]][,rem_nbr_inds, drop = FALSE], 2, sum)}), prod)
        } # else multiply by 1; eg do nothing
        new_msgl = new_msgl / sum(new_msgl)
        conv = conv + sum(abs(new_msgl - MSG_L[[msg_end]][,msg_order[j,2]]))
        
        ## UPDATE PHI
        MSG_L[[msg_end]][,msg_order[j,2]] = new_msgl
        old_msg = PHI[,msg_end]
        PHI[,msg_end] = compute_aggmsg(id = msg_end)
        
        ## UPDATE AUX FIELD
        old_contrib = sapply(1:Nclass, function(x) {sum(edgeps[x,] * old_msg)})
        new_contrib = sapply(1:Nclass, function(x) {sum(edgeps[x,] * PHI[,msg_end])})
        AUX = AUX - old_contrib + new_contrib
      }
    }
    
    # M step
    nodeps_new = apply(PHI, 1, sum) / N
    
    ## Compute normalization factor
    edgelist = find_edgelist(EL)
    edgecont = lapply(seq_len(nrow(edgelist)), function(x) { extract_msg_pair(edgelist[x,1], edgelist[x,2])})
    
    edgeps_new = matrix(0, nrow = Nclass, ncol = Nclass)
    for(m in seq_len(nrow(edgelist))) {
      edgeps_new = edgeps_new + (edgecont[[m]]$contrib_mat / edgecont[[m]]$Z)
    }
    edgeps_new = edgeps_new + t(edgeps_new)
    edgeps_new = edgeps_new / crossprod(t(nodeps_new)) / N
    edgeps_new = edgeps_new / N    
    
    
    ## Compute change in parameter estimates
    delta = sum(abs(edgeps_new - edgeps)) + sum(abs(nodeps_new - nodeps))
    if (verbose > 0) { cat("Iteration ", I, " ----- change in edgeps and nodeps = ", delta, "\n", sep = "") }
    if (verbose > 1) { cat("\t\t Log-likelihood: ", compute_sbm_loglik(class_assign = sapply(1:N, function(x) { order(PHI[,x], decreasing = TRUE)[1] }), adjm = adjm, nodeps = nodeps_new, edgeps = edgeps_new, Nobs = Nobs), "\n") }
    
    ## Update old parameters
    edgeps = edgeps_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    if (is.nan(delta)) { break } # stop if delta is strange; return current result. 
    if (delta < stop_thres) { break } # stop if delta reaches threshold. 
  }
  return(list(nodeps = nodeps, edgeps = edgeps, PHI = t(PHI), 
              classes = sapply(1:N, function(x) { order(PHI[,x], decreasing = TRUE)[1] }), nsteps = I))
}



## TODO: [Documentation-AUTO] Check/fix Roxygen2 Documentation (compute_sbm_loglik)
#' <What does this function do>
#' 
#' @param class_assign temp
#' @param adjm temp
#' @param nodeps temp
#' @param edgeps temp
#' @param Nobs temp
#' 
#' @return temp
#' 
#' @export
#' 
compute_sbm_loglik = function(class_assign, adjm, nodeps, edgeps, Nobs = 1) {
  ## loglik = sum of logprob for classes + sum of log edge probs
  N = length(class_assign)
  resll = sum(log(nodeps[class_assign]))
  caprob = matrix(NA, nrow = N, ncol = N)
  for(i in 1:N) { for (j in 1:N) {
    caprob[i,j] = edgeps[class_assign[i], class_assign[j]]
  }}
  temp = adjm * log(caprob) + (Nobs - adjm) * log(1 - caprob)
  temp2 = temp[lower.tri(temp, diag = FALSE)]
  resll = resll + sum(temp2[is.finite(temp2)], na.rm = TRUE)
  return(resll)
}


