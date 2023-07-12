# Tällä koodilla tutkitaan ARMA-mallin takaisinotannan 
# (bootstrapping) luottamusvälejä. Tutkimuksessa luodaan eripituisia 
# aikasarjoja erilaisilla residuaalien variansseilla eri residuaalien 
# jakaumille. Tutkittavat jakaumat ovat normaalijakauma, hännällinen 
# gammajakauma ja tasajakauma.
# 
# © Jere Niemi, Aalto-yliopisto

# Asetetaan mallille toistettavuus
set.seed(123)

# Tässä tunnetun ARMA(2,3)-aikasarjan parametrit.
alpha <- c(0.5, -0.3)
beta <- c(0.02, -0.15, 0.5)

# Aikasarjan pituudet
time_series_length <- c(100, 1000, 5000)

# Gammajakauman parametrit
shape <- 2
rate <- 0.5

# Tasajakauman parametrit
maximum <- 5
minimum <- -5

# Luottamusväli
confidence_interval <- 0.9

# Asetetaan takaisinotantojen määrä 
bootstrap_rounds <- 1000

# Esitellään dataa työn kannalta "standardiaikasarajalla", 
# jonka pituus on 1000
showcase_ts_length <- 1000

# Funktio luo aikasarjan annettujen parametrien mukaisesti
create_arma_ts <- function(initial_values, a, b, res, len){
  time_series <- matrix(NA, nrow=len, ncol=1)
  time_series[1, 1] <- initial_values[1]
  time_series[2, 1] <- initial_values[2]
  time_series[3, 1] <- initial_values[3]
  
  for(i in 4:len){
    time_series[i, 1] <- 
      a[1] * time_series[i - 1] + a[2] * time_series[i - 2] + 
      b[1] * res[i - 1] + b[2] * res[i - 2] + b[3] * res[i - 3] +
      res[i]
  }
  return(time_series)
} 

# Apufunktio, joka siirtää jakaumaa siten, että jakauman
# keskiarvo on 0.
shift <- function(values, shift_amount){
  result <- values-shift_amount
  return(result)
}

# Funktio luo gammajakauman residuaalit ja siirtää ne siten, 
# että palautettavan jakauman keskiarvo on 0.
generate_gamma_distributed_residuals <- function(
  amount,
  shape,
  rate
){
  gamma_random_numbers <- rgamma(amount, shape, rate)
  return(shift(gamma_random_numbers, shape/rate))
}

# ARMA-mallin takaisinotantafunktio
estimate_arma_parameters <- function(rounds, series, ar, ma){
  ts_min_length <- 10
  arma_coef_amount <- ar + ma
  ARMA <- matrix(NA,nrow=rounds,ncol=arma_coef_amount)
  for(i in 1:rounds){
    series_length <- sample(ts_min_length:length(series), 1)
    starting_point <- sample(1:(length(series) - series_length), 1)
    bootstrap_series <- series[
      starting_point:(starting_point + series_length)
    ]
    bootstrap_fit <- arima(
      bootstrap_series, 
      method="ML", 
      order=c(ar, 0, ma)
    )
    for(j in 1:arma_coef_amount){
      ARMA[i,j] <- bootstrap_fit$coef[j]
    }
  }
  return(ARMA)
}

# Funktio järjestää matriisin alkiot suuruusjärjestykseen
sort_params <- function(series){
  sorted_matrix <- matrix(NA, nrow=nrow(series), ncol=ncol(series))
  for(i in 1:ncol(series)){
    sorted_matrix[,i] <- sort(series[, i])
  }
  return(sorted_matrix)
}

# Funktio palauttaa halutun luottamusvälin arvot
get_confidence_interval <- function(series, interval){
  low_bound_index <- round(((1 - interval)/2) * length(series)) + 1
  up_bound_index<-round((1 - ((1 - interval)/2)) * length(series))
  return(c(series[low_bound_index], series[up_bound_index]))
}

# Funktio laskee luottamusvälit kaikille ARMA-mallin parametreille
calculate_bounds <- function(series, confidence_interval){
  bounds_per_param <- matrix(NA, nrow=ncol(series), ncol=2)
  for(i in 1:ncol(series)){
    bounds_per_param[i,] <- get_confidence_interval(
      series[,i],
      confidence_interval
    )
  }
  return(bounds_per_param)
}

# Apufunktio, joka kertoo, onko oikea parametri 
# luottamusvälin sisällä
parameter_is_in_bounds <- function(bounds, param){
  return(param >= bounds[1] & param <= bounds[2])
}

calculate_confidence_interval_length <- function(bound){
  return(bound[2] - bound[1])
}

is_zero_in_bound <- function(bound){
  return(bound[2] >= 0 & bound[1] < 0)
}

# Esitellään dataa

normally_distributed_residuals <- rnorm(showcase_ts_length)

gamma_distributed_residuals <- 
  generate_gamma_distributed_residuals(
    showcase_ts_length,
    shape,
    rate
  )

uniformly_distributed_residuals <- runif(
  showcase_ts_length,
  min=minimum,
  max=maximum
)

data_normal_distribution <- create_arma_ts(
  rnorm(3),
  alpha,
  beta,
  normally_distributed_residuals,
  showcase_ts_length
)

data_gamma_distribution <- create_arma_ts(
  generate_gamma_distributed_residuals(3, shape, rate),
  alpha,
  beta,
  gamma_distributed_residuals,
  showcase_ts_length
)

data_uniform_distribution <- create_arma_ts(
  runif(3, min=minimum, max=maximum),
  alpha,
  beta,
  uniformly_distributed_residuals,
  showcase_ts_length
)

plot(
  data_normal_distribution,
  type="l",
  main="ARMA(2,3)-aikasarja normaalijakautuneilla residuaaleilla",
  ylab = 'Arvo',
  xlab = 'Indeksi'
)

plot(
  data_gamma_distribution,
  type="l",
  main="ARMA(2,3)-aikasarja gammajakautuneilla residuaaleilla",
  ylab = 'Arvo',
  xlab = 'Indeksi'
)

plot(
  data_uniform_distribution,
  type="l",
  main="ARMA(2,3)-aikasarja tasajakautuneilla residuaaleilla",
  ylab = 'Arvo',
  xlab = 'Indeksi'
)

hist(
  normally_distributed_residuals,
  main="Normaalijakautuneet residuaalit (mean=0, var=1)",
  ylab = 'Määrä',
  xlab = 'Arvo'
)

hist(
  gamma_distributed_residuals,
  main="Gammajakautuneet residuaalit (mean=0, var=8)",
  ylab = 'Määrä',
  xlab = 'Arvo'
)

hist(
  uniformly_distributed_residuals,
  main="Tasajakautuneet residuaalit (mean=0, var=8.3)",
  ylab = 'Määrä',
  xlab = 'Arvo'
)

# params <- estimate_arma_parameters(bootstrap_rounds, data, length(alpha), length(beta))
# sorted_params <- sort_params(params)
# calculated_bounds <- calculate_bounds(sorted_params, confidence_interval)
# print(calculated_bounds)

# Ajetaan iteraatio 100 kertaa ja katsotaan, monta kertaa kukin 
# parametri on luottamusvälin sisällä
iteration_rounds <- 100

times_within_bounds_normal_distribution <- matrix(
  0,
  nrow=length(alpha)+length(beta),
  ncol=1
)
times_within_bounds_gamma_distribution <- matrix(
  0,
  nrow=length(alpha)+length(beta),
  ncol=1
)
times_within_bounds_uniform_distribution <- matrix(
  0,
  nrow=length(alpha)+length(beta),
  ncol=1
)

bounds_normal_distribution <- matrix(
  0,
  nrow=iteration_rounds,
  ncol=2*(length(alpha) + length(beta))
)

bounds_gamma_distribution <- matrix(
  0,
  nrow=iteration_rounds,
  ncol=2*(length(alpha) + length(beta))
)

bounds_uniform_distribution <- matrix(
  0, nrow=iteration_rounds,
  ncol=2*(length(alpha) + length(beta))
)

bound_lengths_normal_distribution <- matrix(
  0,
  nrow=iteration_rounds,
  ncol=length(alpha) + length(beta)
)

bound_lengths_gamma_distribution <- matrix(
  0,
  nrow=iteration_rounds,
  ncol=length(alpha) + length(beta)
)

bound_lengths_uniform_distribution <- matrix(
  0,
  nrow=iteration_rounds,
  ncol=length(alpha) + length(beta)
)

errors_normal_distribution <- 0
errors_gamma_distribution <- 0
errors_uniform_distribution <- 0

zeros_within_bounds_normal_distribution <- matrix(
  0,
  nrow=length(alpha) + length(beta),
  ncol=1
)

zeros_within_bounds_gamma_distribution <- matrix(
  0,
  nrow=length(alpha) + length(beta),
  ncol=1
)
zeros_within_bounds_uniform_distribution <- matrix(
  0,
  nrow=length(alpha) + length(beta),
  ncol=1
)

original_params <- c(alpha, beta)

print(
  paste('Script is about to run with conf =', confidence_interval)
)

print(
  paste('Script is about to run with ts_length =', showcase_ts_length)
)

for(i in 1:iteration_rounds){
  
  normally_distributed_residuals <- rnorm(showcase_ts_length)
  gamma_distributed_residuals <- 
    generate_gamma_distributed_residuals(showcase_ts_length, shape, rate)
  uniformly_distributed_residuals <- 
    runif(showcase_ts_length, min=minimum, max=maximum)
  
  # normally_distributed_residuals <- rnorm(showcase_ts_length, 0, sqrt(2))
  # gamma_distributed_residuals <- generate_gamma_distributed_residuals(showcase_ts_length, 2*shape, rate)
  # uniformly_distributed_residuals <- runif(showcase_ts_length, min=sqrt(2)*minimum, max=sqrt(2)*maximum)
  
  # normally_distributed_residuals <- rnorm(showcase_ts_length, 0, sqrt(2)/2)
  # gamma_distributed_residuals <- generate_gamma_distributed_residuals(showcase_ts_length, 0.5*shape, rate)
  # uniformly_distributed_residuals <- runif(showcase_ts_length, min=(sqrt(2)/2)*minimum, max=(sqrt(2)/2)*maximum)
  
  data_normal_distribution <- create_arma_ts(
    rnorm(3),
    alpha,
    beta,
    normally_distributed_residuals,
    showcase_ts_length
  )
  data_gamma_distribution <- create_arma_ts(
    generate_gamma_distributed_residuals(3, shape, rate),
    alpha,
    beta,
    gamma_distributed_residuals,
    showcase_ts_length
  )
  data_uniform_distribution <- create_arma_ts(
    runif(3, min=minimum, max=maximum),
    alpha, beta,
    uniformly_distributed_residuals,
    showcase_ts_length
  )
  
  # data_normal_distribution <- create_arma_ts(rnorm(3, 0, sqrt(2)/2), alpha, beta, normally_distributed_residuals, showcase_ts_length)
  # data_gamma_distribution <- create_arma_ts(generate_gamma_distributed_residuals(3, 0.5*shape, rate), alpha, beta, gamma_distributed_residuals, showcase_ts_length)
  # data_uniform_distribution <- create_arma_ts(runif(3, min=(sqrt(2)/2)*minimum, max=(sqrt(2)/2)*maximum), alpha, beta, uniformly_distributed_residuals, showcase_ts_length)
  
  succeeded <- FALSE
  maxTries <- 1000
  counter <- 0
  
  # Normaalijakauma
  repeat{
    tryCatch({
      counter <- counter + 1
      if (counter == maxTries) {
        print('Maximum amount of iterations reached. 
              No stationary result found.')
      }
      params <- estimate_arma_parameters(
        bootstrap_rounds,
        data_normal_distribution,
        length(alpha),
        length(beta)
      )
      succeeded <- TRUE
    }, error=function(e) {
      errors_normal_distribution <- errors_normal_distribution + 1
      cat("An error occurred:", conditionMessage(e), "\n")
    })
    if(succeeded | counter >= maxTries){
      break
    }
  }
  sorted_params <- sort_params(params)
  calculated_bounds <- 
    calculate_bounds(sorted_params, confidence_interval)
  for(j in 1:nrow(calculated_bounds)){
    if(parameter_is_in_bounds(calculated_bounds[j,], original_params[j])){
      times_within_bounds_normal_distribution[j, 1] <- 
        times_within_bounds_normal_distribution[j, 1] + 1
    }
    if(is_zero_in_bound(calculated_bounds[j,])){
      zeros_within_bounds_normal_distribution[j, 1] <- 
        zeros_within_bounds_normal_distribution[j, 1] + 1
    }
    bound_lengths_normal_distribution[i, j] <- 
      calculate_confidence_interval_length(calculated_bounds[j,])
    bounds_normal_distribution[i, 2*j - 1] <- calculated_bounds[j, 1]
    bounds_normal_distribution[i, 2*j] <- calculated_bounds[j, 2]
  }
  
  counter <- 0

  # Gammajakauma
  repeat{
    tryCatch({
      counter <- counter + 1
      if (counter == maxTries) {
        print('Maximum amount of iterations reached. 
              No stationary result found.')
      }
      params <- estimate_arma_parameters(
        bootstrap_rounds,
        data_gamma_distribution,
        length(alpha),
        length(beta)
      )
      succeeded <- TRUE
    }, error=function(e) {
      cat("An error occurred:", conditionMessage(e), "\n")
      errors_gamma_distribution <- errors_gamma_distribution + 1
    })
    if(succeeded | counter >= maxTries){
      break
    }
  }
  sorted_params <- sort_params(params)
  calculated_bounds <- calculate_bounds(sorted_params, confidence_interval)
  for(j in 1:nrow(calculated_bounds)){
    if(parameter_is_in_bounds(calculated_bounds[j,], original_params[j])){
      times_within_bounds_gamma_distribution[j, 1] <- 
        times_within_bounds_gamma_distribution[j, 1] + 1
    }
    if(is_zero_in_bound(calculated_bounds[j,])){
      zeros_within_bounds_gamma_distribution[j, 1] <- 
        zeros_within_bounds_gamma_distribution[j, 1] + 1
    }
    bound_lengths_gamma_distribution[i, j] <- 
      calculate_confidence_interval_length(calculated_bounds[j,])
    bounds_gamma_distribution[i, 2*j - 1] <- calculated_bounds[j, 1]
    bounds_gamma_distribution[i, 2*j] <- calculated_bounds[j, 2]
  }

  counter <- 0

  # Tasajakauma
  repeat{
    tryCatch({
      counter <- counter + 1
      if (counter == maxTries) {
        print('Maximum amount of iterations reached. 
              No stationary result found.')
      }
      params <- estimate_arma_parameters(
        bootstrap_rounds,
        data_uniform_distribution,
        length(alpha), length(beta)
      )
      succeeded <- TRUE
    }, error=function(e) {
      cat("An error occurred:", conditionMessage(e), "\n")
      errors_uniform_distribution <- errors_uniform_distribution + 1
    })
    if(succeeded | counter >= maxTries){
      break
    }
  }
  sorted_params <- sort_params(params)
  calculated_bounds <- calculate_bounds(sorted_params, confidence_interval)
  for(j in 1:nrow(calculated_bounds)){
    if(parameter_is_in_bounds(calculated_bounds[j,], original_params[j])){
      times_within_bounds_uniform_distribution[j, 1] <- 
        times_within_bounds_uniform_distribution[j, 1] + 1
    }
    if(is_zero_in_bound(calculated_bounds[j,])){
      zeros_within_bounds_uniform_distribution[j, 1] <- 
        zeros_within_bounds_uniform_distribution[j, 1] + 1
    }
    bound_lengths_uniform_distribution[i, j] <- 
      calculate_confidence_interval_length(calculated_bounds[j,])
    bounds_uniform_distribution[i, 2*j - 1] <- calculated_bounds[j, 1]
    bounds_uniform_distribution[i, 2*j] <- calculated_bounds[j, 2]
  }
  
  print(paste(round((i/iteration_rounds) * 100),'% completed'))
  
}

# Tulosten printtaamista

print(bounds_normal_distribution)
print(bound_lengths_normal_distribution)

for(i in 1:ncol(bound_lengths_normal_distribution)){
  print(paste('index:', i))
  print(max(bound_lengths_normal_distribution[,i]))
  print(min(bound_lengths_normal_distribution[,i]))
  print(mean(bound_lengths_normal_distribution[,i]))
  print(var(bound_lengths_normal_distribution[,i]))
}

for(i in 1:ncol(bounds_normal_distribution)){
  print(paste('index:', i))
  print(mean(bounds_normal_distribution[,i]))
}

for(i in 1:ncol(bounds_gamma_distribution)){
  print(paste('index:', i))
  print(mean(bounds_gamma_distribution[,i]))
}

for(i in 1:ncol(bounds_uniform_distribution)){
  print(paste('index:', i))
  print(mean(bounds_uniform_distribution[,i]))
}

for(i in 1:ncol(bound_lengths_gamma_distribution)){
  print(paste('index:', i))
  print(max(bound_lengths_gamma_distribution[,i]))
  print(min(bound_lengths_gamma_distribution[,i]))
  print(mean(bound_lengths_gamma_distribution[,i]))
  print(var(bound_lengths_gamma_distribution[,i]))
}

for(i in 1:ncol(bound_lengths_uniform_distribution)){
  print(paste('index:', i))
  print(max(bound_lengths_uniform_distribution[,i]))
  print(min(bound_lengths_uniform_distribution[,i]))
  print(mean(bound_lengths_uniform_distribution[,i]))
  print(var(bound_lengths_uniform_distribution[,i]))
}

print(zeros_within_bounds_normal_distribution)

print(bounds_gamma_distribution)
print(bound_lengths_gamma_distribution)
print(zeros_within_bounds_gamma_distribution)

print(bounds_uniform_distribution)
print(bound_lengths_uniform_distribution)
print(zeros_within_bounds_uniform_distribution)

print(times_within_bounds_normal_distribution)
print(times_within_bounds_gamma_distribution)
print(times_within_bounds_uniform_distribution)

print(errors_normal_distribution)
print(errors_gamma_distribution)
print(errors_uniform_distribution)

intervals <- c(80, 82, 84, 86, 88, 90, 92, 94, 96, 98)

n_mat <- matrix(0, nrow=5, ncol=10)

'n_mat[1,] <- c(85, 91, 94, 96, 98, 99, 99, 100, 100, 100)
n_mat[2,] <- c(80, 89, 91, 92, 93, 96, 98, 100, 100, 100)
n_mat[3,] <- c(86, 91, 93, 98, 98, 99, 100, 100, 100, 100)
n_mat[4,] <- c(82, 89, 93, 93, 94, 95, 98, 100, 100, 100)
n_mat[5,] <- c(80, 90, 91, 93, 95, 97, 98, 99, 100, 100)'

'n_mat[1,] <- c(87, 86, 90, 93, 95, 97, 98, 99, 100, 100)
n_mat[2,] <- c(78, 87, 89, 92, 94, 96, 97, 98, 100, 100)
n_mat[3,] <- c(86, 89, 91, 93, 95, 97, 98, 100, 100, 100)
n_mat[4,] <- c(85, 82, 82, 84, 90, 94, 96, 100, 100, 100)
n_mat[5,] <- c(79, 87, 91, 93, 93, 96, 98, 98, 100, 100)'

n_mat[1,] <- c(89, 88, 91, 94, 95, 99, 100, 100, 100, 100)
n_mat[2,] <- c(81, 80, 85, 89, 90, 91, 94, 97, 99, 100)
n_mat[3,] <- c(85, 87, 93, 95, 97, 100, 100, 100, 100, 100)
n_mat[4,] <- c(81, 81, 85, 89, 90, 93, 95, 98, 98, 100)
n_mat[5,] <- c(78, 83, 88, 90, 91, 95, 98, 99, 99, 100)


normal_distribution_intervals_within_bounds <- array(c(
  c(85, 91, 94, 96, 98, 99, 99, 100, 100, 100),
  c(80, 89, 91, 92, 93, 96, 98, 100, 100, 100),
  c(86, 91, 93, 98, 98, 99, 100, 100, 100, 100),
  c(82, 89, 93, 93, 94, 95, 98, 100, 100, 100),
  c(80, 90, 91, 93, 95, 97, 98, 99, 100, 100)
), dim = c(5, 10))
    
gamma_distribution_intervals_within_bounds <- c(
  c(87, 86, 90, 93, 95, 97, 98, 99, 100, 100),
  c(78, 87, 89, 92, 94, 96, 97, 98, 100, 100),
  c(86, 89, 91, 93, 95, 97, 98, 100, 100, 100),
  c(85, 82, 82, 84, 90, 94, 96, 100, 100, 100),
  c(79, 87, 91, 93, 93, 96, 98, 98, 100, 100)
)

uniform_distribution_intervals_within_bounds <- c(
  c(89, 88, 91, 94, 95, 99, 100, 100, 100, 100),
  c(81, 80, 85, 89, 90, 91, 94, 97, 99, 100),
  c(85, 87, 93, 95, 97, 100, 100, 100, 100, 100),
  c(81, 81, 85, 89, 90, 93, 95, 98, 98, 100),
  c(78, 83, 88, 90, 91, 95, 98, 99, 99, 100)
)

print(normal_distribution_intervals_within_bounds)

plot(
  intervals,
  n_mat[1,],
  xlim = c(80, 100),
  ylim = c(75, 100),
  xlab = "Takaisinotantaluottamusväli",
  ylab = "Parametri luottamusvälin sisällä (kertaa)",
  main = "Takaisinotantaluottamusvälin pituuden vaikutus 
    oikean parametrin osumiseen luottamusvälille"
)

points(intervals, n_mat[1,], col="black")
lines(intervals, n_mat[1,], col="black",lty=2)

points(intervals, n_mat[2,], col="red", pch=0)
lines(intervals, n_mat[2,], col="red",lty=2)

points(intervals, n_mat[3,], col="green", pch=2)
lines(intervals, n_mat[3,], col="green",lty=2)

points(intervals, n_mat[4,], col="blue", pch=4)
lines(intervals, n_mat[4,], col="blue",lty=2)

points(intervals, n_mat[5,], col="orange", pch=8)
lines(intervals, n_mat[5,], col="orange",lty=2)

legend(
  "topleft",
  legend=c("ar1","ar2","ma1","ma2","ma3"),
  col=c("black","red","green", "blue", "orange"),
  pch=c(1,0,2,4,8), 
  lty=c(1,1,1,1,1),
  ncol=1
)
