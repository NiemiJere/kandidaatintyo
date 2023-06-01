# Valitaan alkutila 

set.seed(123)

alpha<-c(0.5, -0.3)
beta<-c(0.02, -0.15, 0.5)

# Generoidaan residuaalit
k<-1000

create_initial_arma<-function(initial_values, a, b, res, len){
  time_series<-matrix(NA,nrow=k,ncol=1)
  time_series[1, 1]<-initial_values[1]
  time_series[2, 1]<-initial_values[2]
  time_series[3, 1]<-initial_values[3]
  
  for(i in 4:len){
    time_series[i, 1]<-a[1] * time_series[i - 1] + a[2] * time_series[i - 2] + 
      b[1] * res[i - 1] + b[2] * res[i - 2] + b[3] * res[i - 3] +
      res[i]
  }
  return(time_series)
}

shift<-function(values, shift){
  result<-values-shift
  return(result)
}

generate_gamma_distributed_residuals<-function(amount, shape, rate){
  gamma_random_numbers<-rgamma(amount, shape, rate)
  return (shift(gamma_random_numbers, shape/rate))
}

# Gamma-jakauman parametrit
shape<-2
rate<-0.5

# Tasajakauman parametrit
maximum = 10
minimum = -10

normally_distributed_residuals<-rnorm(k)
gamma_distributed_residuals<-generate_gamma_distributed_residuals(k, shape, rate)
uniformly_distributed_residuals<-runif(k, min=minumum, max=maximum)

'data_normal_distribution<-arima.sim(model=list(ar=alpha, ma=beta), n=k, innov=normally_distributed_residuals)
data_gamma_distribution<-arima.sim(model=list(ar=alpha, ma=beta), n=k, innov=gamma_distributed_residuals)
data_uniform_distribution<-arima.sim(model=list(ar=alpha, ma=beta), n=k, innov=uniformly_distributed_residuals)'

data_normal_distribution<-create_initial_arma(rnorm(3), alpha, beta, normally_distributed_residuals, k)
data_gamma_distribution<-create_initial_arma(generate_gamma_distributed_residuals(3, shape, rate), alpha, beta, gamma_distributed_residuals, k)
data_uniformly_distribution<-create_initial_arma(runif(3, min=minumum, max=maximum), alpha, beta, uniformly_distributed_residuals, k)

data<-data_normal_distribution

# Bootstrap-funktio
create_arma_parameters<-function(rounds, series, arma_coef_amount){
  ts_min_length<-10
  ARMA<-matrix(NA,nrow=rounds,ncol=arma_coef_amount)
  for(i in 1:rounds){
    series_length<-sample(ts_min_length:length(series), 1)
    starting_point<-sample(1:(length(series)-series_length), 1)
    bootstrap_series<-series[starting_point:(starting_point+series_length)]
    bootstrap_fit<-arima(bootstrap_series,method="ML",order=c(2, 0, 3))
    for(j in 1:arma_coef_amount){
      ARMA[i,j]<-bootstrap_fit$coef[j]
    }
  }
  return (ARMA)
}

sort_params<-function(series){
  sorted_matrix<-matrix(NA,nrow=nrow(series),ncol=ncol(series))
  for(i in 1:ncol(series)){
    sorted_matrix[,i]<-sort(series[, i])
  }
  return (sorted_matrix)
}

get_confidence_interval<-function(series, interval){
  low_bound_index<-round(((1-interval)/2)*length(series))+1
  up_bound_index<-round((1-((1-interval)/2))*length(series))
  return (c(series[low_bound_index], series[up_bound_index]))
}

calculate_bounds<-function(series, confidence_interval){
  bounds_per_param<-matrix(NA,nrow=ncol(series),ncol=2)
  for(i in 1:ncol(series)){
    bounds_per_param[i,]<-get_confidence_interval(series[,i], confidence_interval)
  }
  return (bounds_per_param)
}

parameter_is_in_bounds<-function(bounds, param){
  return (param>=bounds[1] & param<=bounds[2])
}

# Luodaan ARMA(2,3)-mallille matriisi hyödyntäen create_arma_parameters-funktiota
bootstrap_rounds<-1000

# Koska käsittelemme tässä työssä ARMA(2,3)-mallia, parametrejä on yhteensä 5
ar_params<-2
ma_params<-3
coefs<-ar_params+ma_params

# Valitaan aluksi luottamusväliksi 90% luottamusväli
confidence_interval<-0.9

params<-create_arma_parameters(bootstrap_rounds, data, coefs)
sorted_params<-sort_params(params)
calculated_bounds<-calculate_bounds(sorted_params, confidence_interval)
print(calculated_bounds)

# Ajetaan yllä oleva 100 kertaa ja katsotaan, monta kertaa kukin parametri on luottamusvälin sisällä
iteration_rounds<-100

times_within_bounds<-matrix(0,nrow=coefs,ncol=1)
original_params<-c(alpha, beta)

for(i in 1:iteration_rounds){
  succeeded<-FALSE
  repeat{
    tryCatch({
      params<-create_arma_parameters(bootstrap_rounds, data, coefs)
      succeeded<-TRUE
    }, error=function(e) {
      cat("An error occurred:", conditionMessage(e), "\n")
    })
    if(succeeded){
      break
    }
  }
  sorted_params<-sort_params(params)
  calculated_bounds<-calculate_bounds(sorted_params, confidence_interval)
  for(j in 1:nrow(calculated_bounds)){
    if(parameter_is_in_bounds(calculated_bounds[j,], original_params[j])){
      times_within_bounds[j, 1]<-times_within_bounds[j, 1]+1
    }
  }
  print(paste(round((i/iteration_rounds)*100),'% completed'))
}

print("Parameters with normally distributed residuals")
print(times_within_bounds)








