library(ggplot2)
library(shinyWidgets)
library(data.table)

# Get Interpretatable Vals

getInterpretableVals <- function(p.g_g, 
                                 p.g_icu, 
                                 p.g_d, 
                                 p.icu_g,
                                 p.icu_icu, 
                                 p.icu_v,
                                 p.v_icu, 
                                 p.v_v){

  vec1 <- c(1 - p.g_g, -p.g_icu, 0)
  vec2 <- c(-p.icu_g, 1 - p.icu_icu, -p.icu_v)
  vec3 <- c(0, -p.v_icu, 1 - p.v_v)

  ls <- rbind(vec1, vec2, vec3)
  rs <- rbind(1, 1, 1)
  
  sol <- t(solve(ls, rs))
  colnames(sol) <- c('m.g', 'm.icu', 'm.v')
  
  return(sol[1,])
  
  
}

# Split gamma for I -> H and I -> D
findGammas <- function(per.hosp, 
                       hosp.delay.time, 
                       illness.length){
  
  per.nonhosp <- 1 - per.hosp
  
  avg.time.i <- (per.hosp * hosp.delay.time) + (per.nonhosp * illness.length)
  
  gamma.i_d <- per.nonhosp * (1 / avg.time.i)
  gamma.i_h <- per.hosp * (1 / avg.time.i)
  
  gamma.tot <- gamma.i_d + gamma.i_h
  
  gamma.list <- list(
    'gamma.d' = gamma.i_d,
    'gamma.h' = gamma.i_h, 
    'gamma.tot' = gamma.tot
  )
  
  return(gamma.list)
}

# Transition Matrix
createTransition <- function(p.g_g, 
                             p.g_icu, 
                             p.g_d, 
                             p.icu_g,
                             p.icu_icu, 
                             p.icu_v,
                             p.v_icu, 
                             p.v_v, 
                             p.v_m,
                             p.g_m = 0,
                             p.icu_m = 0){

  tmat <- matrix(c(p.g_g, p.g_icu, 0, p.g_d, p.g_m,
                   p.icu_g, p.icu_icu, p.icu_v, 0, p.icu_m,
                   0, p.v_icu, p.v_v, 0, p.v_m,
                   0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), nrow = 5, byrow = TRUE)
  

  return(tmat)
}

# Run Markov Simulation 
runMarkov <- function(new.hosp.vec, trans.mat){
  
  # h, icu, vent, d, m 
  df.final <- data.frame()
  
  step.vec <- c(0, 0, 0, 0, 0)
  
  for (hosp.val in new.hosp.vec){
    step.vec <- as.vector(step.vec %*% trans.mat)
    step.vec[1] <- step.vec[1] + hosp.val 
    df.final <- rbind(df.final, step.vec)
  }
  
  colnames(df.final) <- c('h', 'icu', 'vent', 'd', 'death')
  return(df.final)
}

# Runs SIR, combines with Markov
SIR <- function(S0, I0, R0, beta.vector, gamma_d, gamma_h, trans.mat, num.days) {
  
  # initialize S, I, R 
  S <- I <- R <- hosp <- rep(NA_real_, num.days)
  S[1] <- S0
  I[1] <- I0
  R[1] <- R0
  N = S[1] + I[1] + R[1]
  hosp[1] <- 0
  
  gamma <- gamma_d + gamma_h
  
  # run SIR model 
  for (tt in 1:(num.days - 1)) {
    beta <- beta.vector[tt]
    S[tt + 1] <-  -beta * S[tt] * I[tt] / N                 + S[tt]
    I[tt + 1] <-   beta * S[tt] * I[tt] / N - gamma * I[tt]  + I[tt]
    R[tt + 1] <-                          gamma_d * I[tt]  + R[tt]
    hosp[tt + 1] <- gamma_h * I[tt]
  }
  
  df.return <- data.frame(days = 1:num.days, S, I, R, hosp)
  markov.df <- runMarkov(hosp, trans.mat)
  markov.df$days <- 1:num.days 
  
  df.return <- merge(df.return, markov.df, by = 'days')
  df.return$R <- df.return$R + df.return$d
  df.return$d <- NULL
  df.return$h <- NULL
  
  return(df.return)
  
}


# finds current estimates of the number of active infections, 
# number recovered, and number 
find.curr.estimates = function(S0, beta.vector, gamma_d, gamma_h, trans.mat,
                               num.days, num.actual, start.inf = 1){
  
  # starting number of susceptible people
  start.susc <- S0 - start.inf
  start.res <- 0 
  
  SIR.df = SIR(start.susc, start.inf, start.res, beta.vector, gamma_d, 
               gamma_h, trans.mat, num.days)

  # find the difference between hospitalized column and the currently hospitalized number
  SIR.df$diff_proj <- abs(SIR.df$hosp - num.actual)
    
  # hacky 
  # the # hospitalized will be achieved twice according to model 
  # first as the hospitalizations go up, and second as the hospitalizations go down 
  # we want to find the day
  hosp.numbers <- SIR.df$hosp
  hosp.change <- hosp.numbers[2:length(hosp.numbers)] - hosp.numbers[1:length(hosp.numbers) - 1]
  hosp.change <- c(0, hosp.change)
  SIR.df$hosp.change <- hosp.change
  
  curr.day.df <- SIR.df[SIR.df$hosp.change > 0,]
  curr.day.df <- curr.day.df[curr.day.df$diff_proj == min(curr.day.df$diff_proj, na.rm = TRUE),] 
  
  curr.day <- as.integer(curr.day.df$day)
  infection.estimate <- as.integer(curr.day.df$I)
  susceptible.estimate <- as.integer(curr.day.df$S)
  recovered.estimate <- as.integer(curr.day.df$R)
  
  curr.day.list <- list(
    curr.day = curr.day,
    infection.estimate = infection.estimate, 
    susceptible.estimate = susceptible.estimate, 
    recovered.estimate = recovered.estimate
  )
  
  return(curr.day.list)
}


# gets beta based on doubling time, gamma, and N
getBetaFromDoubling <- function(doubling.time, gamma) {
  g <- 2^(1/doubling.time) - 1
  beta <- g + gamma
  return(beta)
}

getBetaFromRe <- function(Re, gamma) {
  beta <- Re * gamma
  return(beta)
}

