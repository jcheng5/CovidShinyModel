library(ggplot2)
library(shinyWidgets)
library(data.table)

# TODO: figure out how to add function documentation

# runs SIR model simulation, and then calculates hospitalizations, 
# icu numbers, and ventilator numbers 

createTransition <- function(p.g_g, 
                             p.g_icu, 
                             p.g_d, 
                             p.icu_g,
                             p.icu_icu, 
                             p.icu_v,
                             p.v_icu, 
                             p.v_v, 
                             p.v_m){
  tmat <- matrix(c(p.g_g, p.g_icu, 0, p.g_d, 0,
                   p.icu_g, p.icu_icu, p.icu_v, 0, 0,
                   0, p.v_icu, p.v_v, 0, p.v_m,
                   0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), byrow = TRUE)
  
  return(tmat)
}


runMarkov <- function(new.hosp.vec, trans.mat){
  
  # h, icu, vent, d, m 
  df.final <- data.frame()
  
  step.vec <- c(0, 0, 0, 0, 0)
  
  for (hosp.val in new.hosp.vec){
    step.vec <- as.vector(step.vec %*% trans.mat)
    step.vec[1] <- step.vec[1] + hosp.val 
    df.final <- rbind(df.final, step.vec)
  }
  
  colnames(df.final) <- c('h', 'icu', 'vent', 'd', 'm')
  return(df.final)
}


SIR <- function(S0, I0, R0, beta.vector, gamma_r, gamma_h, trans.mat, num.days) {
  
  # initialize S, I, R 
  S <- I <- R <- H <- rep(NA_real_, num.days)
  S[1] <- S0
  I[1] <- I0
  R[1] <- R0
  hosp[1] <- 0
  
  gamma <- gamma_r + gamma_h
  
  # run SIR model 
  for (tt in 1:(num.days - 1)) {
    beta <- beta.vector[tt]
    S[tt + 1] <-  -beta * S[tt] * I[tt]                  + S[tt]
    I[tt + 1] <-   beta * S[tt] * I[tt] - gamma * I[tt]  + I[tt]
    R[tt + 1] <-                          gamma_r * I[tt]  + R[tt]
    hosp[tt + 1] <- gamma_h * I[tt]
  }
  
  df.return <- data.frame(days = 1:num.days, S, I, R, hosp)
  markov.df <- runMarkov(H, trans.mat)
  markov.df$days <- 1:num.days 
  
  df.return <- merge(df.return, markov.df, by = 'days')
  df.return$R <- df.return$R + df.return$d
  df.return$d <- NULL
  df.return$h <- NULL
  
  return(df.return)
  
}


# finds current estimates of the number of active infections, 
# number recovered, and number 
find.curr.estimates = function(S0, beta.vector, gamma_r, gamma_h, trans.mat,
                               num.days, num.actual, start.inf = 1){
  
  # starting number of susceptible people
  start.susc <- S0 - start.inf
  start.res <- 0 
  
  SIR.df = SIR(start.susc, start.inf, start.res, beta.vector, gamma_r, 
               gamma_h, trans.mat, num.days)

  # find the difference between hospitalized column and the currently hospitalized number
  SIR.df$diff_proj <- abs(SIR.df$hosp - num_actual)
    
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

# gets doubling time based on R0 and gamma 
doubleTime <- function(R0, gamma) {
  1 / log2(R0 * gamma - gamma + 1)
}

# gets beta based on doubling time, gamma, and S0
getBeta <- function(doubling.time, gamma, S0) {
  g <- 2^(1/doubling.time) - 1
  beta <- (g + gamma) / S0
  return(beta)
}

