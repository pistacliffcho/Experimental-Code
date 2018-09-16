gmmMonoEM = setRefClass("gmmMonoEM", 
                        fields = c("data", "Omega", "miss_count",
                                   "ss_start",
                                   "raw_weights", 
                                   "comp_probs", "obs_dens", "obs_probs", 
                                   "k", "nRow", "nCol", "comp_pars"))

gmmMonoEM$methods(
  compute_comp_mle = function(comp_i){
    this_w = raw_weights * obs_probs[,comp_i]
    this_mle = computeMLE(data, this_w, Omega, ss_start)
    comp_pars[[comp_i]] <<- this_mle
  }
)

gmmMonoEM$methods(
  compute_comp_dens = function(comp_i){
    this_comp_pars = comp_pars[[comp_i]]
    mu = this_comp_pars$mu
    chol = this_comp_pars$chol
    log_dens = mono_dmvnorm(data, mu, chol, miss_count) * raw_weights
    obs_dens[,comp_i] <<- exp(log_dens)
  }
)

gmmMonoEM$methods(
  row_dens2probs = function(){
    rescale_dens = row_mult(obs_dens, comp_probs)
    rescale_dens = col_mult(rescale_dens, 1 / rowSums(rescale_dens) )
    obs_probs <<- rescale_dens
  }
)

gmmMonoEM$methods(
  update_p = function(){
    p_vec = numeric(k)
    for(i in 1:k){
      p_vec[i] = sum( obs_probs[,i] * raw_weights / sum(raw_weights) )
    }
    comp_probs <<- p_vec
  }
)

gmmMonoEM$methods(
  llk = function(){
    probs = rowSums(row_mult(obs_dens, comp_probs) )
    ans = sum(raw_weights * log(probs))
    return(ans)
  }
)

gmmMonoEM$methods(
  E_step = function(){
    for(i in 1:k){ compute_comp_dens(i) }
    row_dens2probs()
    update_p()
  }
)

gmmMonoEM$methods(
  M_step = function(){
    for(i in 1:k){ compute_comp_mle(i) }
  }
)

make_init_obs_probs = function(data, k){
  nRow = nrow(data)
  ans = matrix(runif(nRow * k), nrow = nRow, ncol = k)
  ans = ans/rowSums(ans)
  return(ans)
}

alt_chol = function(chol_mat){
  invChol = solve(chol_mat)
  diag(invChol) = diag(invChol) + 0.01
  chol_mat = solve(invChol)
  return(chol_mat)
}

make_gmmMonoEM = function(data, k = 2, w = NULL){
  ans = gmmMonoEM()
  preppedData = dataPrep(data)
  ans$Omega = preppedData$Omega
  ans$data = preppedData$data
  ans$miss_count = rowSums(is.na(ans$data))
  ans$k = k
  ans$obs_probs = make_init_obs_probs(data, k)
  ans$obs_dens = ans$obs_probs + NA
  ans$nRow = nrow(data)
  ans$nCol = ncol(data)
  ans$comp_pars = list()
  if(is.null(w)){ w = rep(1, ans$nRow) }
  ans$raw_weights = w
  
  ans$update_p()
  
  c_means = colMeans(data, na.rm = T)
  c2_means = colMeans(data^2, na.rm = T)
  sd_0 = sqrt( c2_means - c_means^2 )
  ans$ss_start = diag(sd_0, length(sd_0) )
  
  
  return(ans)
}
