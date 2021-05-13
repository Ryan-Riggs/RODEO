###Error Functions.

##RRMSE: Relative root mean square error. RRMSE^2 = (rBias^2) + (SDRR^2)
RRMSE = function(m, o){
  error = as.numeric(m - o)
  relative_residual = error/as.numeric(o)
  sqr_residual = (relative_residual^2)
  mean_sqr = mean(sqr_residual, na.rm = TRUE)
  rrmse = sqrt(mean_sqr)
  return(rrmse*100)
}


##NRMSE: Normalized root mean square error.
NRMSE = function(m, o){
  error = as.numeric(m - o)
  mean_obs = mean(as.numeric(o), na.rm = TRUE)
  relative_residual = error/as.numeric(mean_obs)
  sqr_residual = (relative_residual^2)
  mean_sqr = mean(sqr_residual, na.rm = TRUE)
  nrmse = sqrt(mean_sqr)
  return(nrmse*100)
}


##rBias
rBias = function(m, o){
  error = as.numeric(m - o)
  mean_obs = mean(as.numeric(o), na.rm = TRUE)
  out = error/mean_obs
  rBias = mean(out, na.rm = TRUE)
  return(rBias*100)
}

#R value
Rvalue = function(m, o){
  pearson = cor.test(o, m, method = "pearson")
  r = pearson$estimate
  return(r)
}
