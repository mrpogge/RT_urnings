################################################################################
# This script contains the functions used multiple times through the results
################################################################################

#calculating hitting times
calculate_ht = function(v, lower, upper){
  lower_true = v>lower
  upper_true = v<upper
  bounded = lower_true & upper_true
  #first time when bounded is true
  if(sum(bounded) > 0){
    return(min(which(bounded)))
  } else {
    return(500)
  }
}

#calculating hitting times
calculate_ht_MSE = function(v, upper){
  upper_true = v<upper
  #first time when bounded is true
  if(any(upper_true)){
    return(min(which(upper_true)))
  } else {
    return(500)
  }
}

transform_logit = function(v){
  return(exp(v)/(1+exp(v)))
}

transform_expit = function(v){
  return(log(v/(1-v)))
}

twoPL=function(theta,delta,a){
  theta=log(theta/(1-theta))
  delta=log(delta/(1-delta))
  P=1/(1+exp(-a*(theta-delta)))
  P[delta==theta]=0.5
  return(P)
}

coverage = function(v, lower, upper){
  return(mean(v>lower & v<upper))
}