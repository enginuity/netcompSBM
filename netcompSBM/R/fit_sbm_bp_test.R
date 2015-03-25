
#' @param adjm Input adjacency matrix. This can contain NAs, values > 1, as a total count over a number of networks.
#' @param Nobs Number of network observations
#' @param Nclass Number of classes in the SBM
#' @param Niter Max number of EM steps
#' @param verbose Lots of output?
#' @param stop_thres Stopping threshold (will stop before max Niter, if change in probabilty estimate matrix is smaller than this value)
fit_SBM_BP = function(adjm, Nobs = 1, Nclass, Niter = 100, verbose = 1, stop_thres = 0.000001) {
  if(FALSE) {
    adjm = graph
    Nobs = 1
    Nclass = 5
    stop_thres = 0.0001
  }  
  ## helper functions
  compute_aggmsg = function(id = NULL) {
    # assumes MSG_L and MSG exists
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
    return(sapply(1:Nclass, function(x){ sum(cab[x,] * MSG)/N } )) ## extra divison by N??) ## equation 27
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
  
  ## SETUP
  N = nrow(adjm)
  nodeps = rep(1/Nclass, length = Nclass)
  edgeps = matrix(sum(adjm, na.rm = TRUE) / (N * (N-1) / 2) * runif(Nclass * Nclass, min = 0.1, max = 0.9), nrow = Nclass)
  edgeps = symmetrize_mat(edgeps)
  
  cab = edgeps * N
  
  ## Get edges in adjm. 
  G = graph.adjacency(adjm, mode = "undirected")
  EL = get.adjlist(G)
  
  ## Do inner loop
  for(I in 1:Niter) {
    
    # Do setup for the E-step. 
    MSG_L = lapply(EL, FUN = function(x) { ## stores the in messages. 
      if(length(x) > 0) {
        z = matrix(runif(length(x) * Nclass), nrow = Nclass) 
        z = sapply(1:ncol(z), function(y) {z[,y] / sum(z[,y])})
        return(z)
      } else { return(NULL) }
    })
    
    
    MSG = matrix(0, nrow = Nclass, ncol = N) # stores nodal overall value
    AUX = rep(0, times = Nclass) # auxiliary field
    MSG = compute_aggmsg()
    AUX = compute_aux()
    
    max_estep_tries = max(floor(sqrt(N) / 2), 5) # max algo steps 
    tries = 0
    conv = stop_thres + 10
    
    ## Run BP algorithm. 
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
        
        ## UPDATE MSG
        MSG_L[[msg_end]][,msg_order[j,2]] = new_msgl
        old_msg = MSG[,msg_end]
        MSG[,msg_end] = compute_aggmsg(id = msg_end)
        
        ## UPDATE AUX FIELD
        old_contrib = sapply(1:Nclass, function(x) {sum(1/N * cab[x,] * old_msg)})
        new_contrib = sapply(1:Nclass, function(x) {sum(1/N * cab[x,] * MSG[,msg_end])})
        AUX = AUX - old_contrib + new_contrib
      }
    }
    
    ## M step
    nodeps_new = apply(MSG, 1, sum) / N
    
    ## Compute normalization factor
    edgelist = find_edgelist(EL)
    edgecont = lapply(seq_len(nrow(edgelist)), function(x) { extract_msg_pair(edgelist[x,1], edgelist[x,2])})
    
    edgeps_new = matrix(0, nrow = Nclass, ncol = Nclass)
    for(m in seq_len(nrow(edgelist))) {
      edgeps_new = edgeps_new + (edgecont[[m]]$contrib_mat / edgecont[[m]]$Z)
    }
    edgeps_new = (edgeps_new + t(edgeps_new))
    cab_new = edgeps_new / (t(t(nodeps_new)) %*% t(nodeps_new)) / N
    
    ## Compute change in edge probability matrix
    delta = sum(abs(cab_new - cab))
    if (verbose > 0) { cat("Iteration ", I, " ----- change in epmat = ", delta, "\n", sep = "") }
    
    ## Update old parameters
    cab = cab_new
    nodeps = nodeps_new
    
    ## Stop if threshold is met
    if (delta < stop_thres) { break }
    #if (any(nodeps > 0.9)) { stop("Local Solution") }
  }
  
  return(list(nodeprobs = nodeps, edgeprobs = cab, classes = sapply(1:N, function(x) {
    order(MSG[,x], decreasing = TRUE)[1]
  })))
}
