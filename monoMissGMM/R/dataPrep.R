dataPrep = function(data){
  data_mat = as.matrix(data)
  isMissing = is.na(data_mat)
  missCount = rowSums(isMissing)
  
  nCol = ncol(data_mat)
  # Set containing missingness pattern
  Omega = list()
  for(i in seq_len(nCol)){
    Omega[[i]] = which(missCount == (i-1) ) - 1
  }
  
  ans = list(data = data_mat, 
             Omega = Omega)
  return(ans)
}