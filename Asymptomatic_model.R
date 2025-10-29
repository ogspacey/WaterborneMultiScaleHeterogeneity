### Asymptomatic model
### Oliver G. Spacey
### 2025.10.27

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load required packages
library(deSolve)   # differential equations
library(ggplot2)   # data visualisation
library(tidyverse) # data wrangling
library(reshape2)  # data wrangling
library(nleqslv)  # solve nonlinear equations

# Initial Model -----------------------------------------------------------
# Create function to run asymptomatic model
asmod <- function(t, y, parameters){
  # Pull state variables from y vector
  S1 = y[1] # Susceptible class without prior infection or vaccination
  IS = y[2] # Symptomatic infectious class
  R1 = y[3] # Recovered class following primary infection or vaccination
  S2 = y[4] # Susceptible class following previous infection or vaccination
  IA = y[5] # Asymptomatic/mild infectious class
  R2 = y[6] # Recovered class following asymptomatic infection
  W  = y[7] # Bacterial concentration in the water
  
  # Pull parameter values from the input vector 
  mu     = parameters["mu"]     # Birth and death rate of each compartment
  N      = parameters["N"]      # Total population size
  betaW  = parameters["betaW"]  # Indirect transmission rate
  betaDS = parameters["betaDS"] # Direct transmission rate from symptomatic individuals
  betaDA = parameters["betaDA"] # Direct transmission rate from asymptomatic individuals
  u      = parameters["u"]      # Vaccination rate
  pS     = parameters["pS"]     # Proportion of initial infections which are symptomatic
  gammaS = parameters["gammaS"] # Recovery rate for symptomatic infections
  gammaA = parameters["gammaA"] # Recovery rate for symptomatic infections
  sigma1 = parameters["sigma1"] # Rate of loss of immunity for recovered class 1
  sigma2 = parameters["sigma2"] # Rate of loss of immunity for recovered class 2
  alphaS = parameters["alphaS"] # Shedding rate from sympomatic individuals
  alphaA = parameters["alphaA"] # Shedding rate from asympomatic individuals
  muW    = parameters["muW"]    # Rate of growth/loss of bacteria in the water
  
  # Define equations
  dS1 = mu * N - (betaW * S1 * W + betaDS * S1 * IS + betaDA * S1 * IA) / N - (u + mu) * S1
  dIS = pS * (betaW * S1 * W + betaDS * S1 * IS + betaDA * S1 * IA) / N - (gammaS + mu) * IS
  dR1 = u * S1 + gammaS * IS - (sigma1 + mu) * R1
  dS2 = sigma1 * R1 + sigma2 * R2 - (betaW * S2 * W + betaDS * S2 * IS + betaDA * S2 * IA) / N - mu * S2
  dIA = (1 - pS) * (betaW * S1 * W + betaDS * S1 * IS + betaDA * S1 * IA) / N + (betaW * S2 * W + betaDS * S2 * IS + betaDA * S2 * IA) / N - (gammaA + mu) * IA
  dR2 = gammaA * IA - (sigma2 + mu) * R2
  dW  = alphaS * IS + alphaA * IA + muW * W
  res = c(dS1, dIS, dR1, dS2, dIA, dR2, dW)
  
  # Return list of gradients
  list(res)
}

# Set times to run simulation for
times = seq(0, 365*5, by = 1)

# Set parameters for initial model
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0.00001,          # Vaccination rate
          pS     = 0.5,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Set population size
N <- 10000

# Set initial conditions - epidemic conditions
start <- c(S1 = 0.99 * N,  # Susceptible class without prior infection or vaccination
           IS = 0.005 * N, # Symptomatic infectious class
           R1 = 0, # Recovered class following primary infection or vaccination
           S2 = 0, # Susceptible class following previous infection or vaccination
           IA = 0.005 * N, # Asymptomatic/mild infectious class
           R2 = 0, # Recovered class following asymptomatic infection
           W  = 0) # Bacterial concentration in the water

# Simulate dynamics
out <- ode(y = start, times = times, func = asmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$S1, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0, 10000), cex.lab = 1.5)
lines(x = out$time, y = out$IS, col = "red") 
lines(x = out$time, y = out$R1, col = "blue")
lines(x = out$time, y = out$S2, col = "grey")
lines(x = out$time, y = out$IA, col = "pink") 
lines(x = out$time, y = out$R2, col = "lightblue")
legend(x = 0.8 * max(out$time), y = 0.9 * N,                                
       legend = c("S1", "IS", "R1", "S2", "IA", "R2"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend

# Plot just symptomatic infected individuals
plot(x = out$time, y = out$IS, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", col = 'red', cex.lab = 1.5)

# Plot reservoir
plot(x = out$time, y = log(out$W), col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "log(Bacterial concentration in water)")

# Check constant population size
out$Ntot <- rowSums(out[, c("S1", "IS", "R1", "S2", "IA", "R2")])
plot(out$time, out$Ntot, type='l')

# Calculate R0 --------------------------------------------
# Create function to calculate R0
calc_R0 <- function(paras) {
  # Extract parameters
  mu      <- as.numeric(paras["mu"])
  u       <- as.numeric(paras["u"])
  betaDS  <- as.numeric(paras["betaDS"])
  betaDA  <- as.numeric(paras["betaDA"])
  betaW   <- as.numeric(paras["betaW"])
  gammaS  <- as.numeric(paras["gammaS"])
  gammaA  <- as.numeric(paras["gammaA"])
  muW     <- as.numeric(paras["muW"])
  alphaS  <- as.numeric(paras["alphaS"])
  alphaA  <- as.numeric(paras["alphaA"])
  pS      <- as.numeric(paras["pS"])
  
  # Susceptible fraction at disease-free equilibrium
  s <- mu / (u + mu)
  
  # Direct transmission component
  R0_direct <- s * (pS * betaDS / (gammaS + mu) +
                    (1 - pS) * betaDA / (gammaA + mu))
  
  # Environmental transmission component
  R0_env <- s * (betaW / (-muW)) * (
    pS * alphaS / (gammaS + mu) +
      (1 - pS) * alphaA / (gammaA + mu)
  )
  
  # Total R0
  R0 <- R0_direct + R0_env
  
  return(R0)
}

# Calculate for current model
calc_R0(paras)

# Break down contributions from each class of susceptibles to Re
calc_Re_breakdown <- function(out, paras) {
  with(as.list(paras), {
    s <- mu / (u + mu)
    
    # Per-susceptible contribution from S1
    direct_S1 <- pS * betaDS / (gammaS + mu) + (1-pS) * betaDA / (gammaA + mu)
    env_S1    <- (betaW / (-muW)) * (pS * alphaS / (gammaS + mu) + (1-pS) * alphaA / (gammaA + mu))
    Re_S1 <- (out$S1 / N) * (direct_S1 + env_S1) / s
    
    # Per-susceptible contribution from S2
    direct_S2 <- betaDS / (gammaS + mu) + betaDA / (gammaA + mu)
    env_S2    <- (betaW / (-muW)) * (alphaS / (gammaS + mu) + alphaA / (gammaA + mu))
    Re_S2 <- (out$S2 / N) * (direct_S2 + env_S2) / s
    
    data.frame(time = out$time, Re_S1 = Re_S1, Re_S2 = Re_S2, Re_total = Re_S1 + Re_S2)
  })
}

# Plot contributions
Re_break <- calc_Re_breakdown(out, paras)

plot(Re_break$time, Re_break$Re_total, type="l", col="black", lwd=2,
     xlab="Time (days)", ylab="R_e(t)", main="R_e(t) contributions by susceptible class",
     ylim = c(0,3))
lines(Re_break$time, Re_break$Re_S1, col="blue", lwd=2)
lines(Re_break$time, Re_break$Re_S2, col="green", lwd=2)
abline(h=1, lty=2, col="red")
legend("topright", legend=c("Total Re","S1 contribution","S2 contribution"),
       col=c("black","blue","green"), lwd=2, bty="n")


# Long term dynamics ------------------------------------------------------
# Run for long time
# Set times to run simulation for
times = seq(0, 365*300, by = 1)

# Set parameters for initial model
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0.00001,          # Vaccination rate
          pS     = 0.5,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Set population size
N <- 10000

# Set initial conditions - epidemic conditions
start <- c(S1 = 0.99 * N,  # Susceptible class without prior infection or vaccination
           IS = 0.005 * N, # Symptomatic infectious class
           R1 = 0, # Recovered class following primary infection or vaccination
           S2 = 0, # Susceptible class following previous infection or vaccination
           IA = 0.005 * N, # Asymptomatic/mild infectious class
           R2 = 0, # Recovered class following asymptomatic infection
           W  = 0) # Bacterial concentration in the water

# Simulate dynamics
out <- ode(y = start, times = times, func = asmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$S1, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0, 10000), cex.lab = 1.5)
lines(x = out$time, y = out$IS, col = "red") 
lines(x = out$time, y = out$R1, col = "blue")
lines(x = out$time, y = out$S2, col = "grey")
lines(x = out$time, y = out$IA, col = "pink") 
lines(x = out$time, y = out$R2, col = "lightblue")
legend(x = 0.8 * max(out$time), y = 0.9 * N,                                
       legend = c("S1", "IS", "R1", "S2", "IA", "R2"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend

# Plot just symptomatic infected individuals
plot(x = out$time, y = out$IS, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", col = 'red', cex.lab = 1.5)

# Plot reservoir
plot(x = out$time, y = log(out$W), col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "log(Bacterial concentration in water)")

# Calculate equilibria ----------------------------------------------------
# Function to solve for equilibrium
equilibrium_system <- function(x, paras){
  S1 <- x[1]; S2 <- x[2]; IS <- x[3]; IA <- x[4]; R1 <- x[5]; R2 <- x[6]; W <- x[7]
  with(as.list(paras), {
    lambda <- (betaDS * IS + betaDA * IA + betaW * W) / N
    
    F <- numeric(7) # Set all values to 0
    F[1] <- mu * N - lambda * S1 - (u + mu) * S1
    F[2] <- sigma1 * R1 + sigma2 * R2 - lambda * S2 - mu * S2
    F[3] <- pS * lambda * S1 - (gammaS + mu) * IS
    F[4] <- (1 - pS) * lambda * S1 + lambda * S2 - (gammaA + mu) * IA
    F[5] <- u * S1 + gammaS * IS - (sigma1 + mu) * R1
    F[6] <- gammaA * IA - (sigma2 + mu) * R2
    F[7] <- alphaS * IS + alphaA * IA + muW * W
    return(F)
  })
}

# Take initial guess - based on visual
x0 <- c(S1 = paras["N"]*0.01, S2 = paras["N"]*0.5, IS = paras["N"]*0.01, IA = paras["N"]*0.09,
        R1 = paras["N"]*0.09, R2 = 0.3, W = 100000)

# Solve using nleqslv
sol <- nleqslv(x0, equilibrium_system, paras = paras)

# Extract solution for endemic equilibrium
endemic_eq <- data.frame(
  S1 = sol$x[1],
  S2 = sol$x[2],
  IS = sol$x[3],
  IA = sol$x[4],
  R1 = sol$x[5],
  R2 = sol$x[6],
  W  = sol$x[7]
)
print(endemic_eq)

# Verify using ode run
print(tail(out, 1))


# Calculate equilibrium asymptomatic ratio
as_rat <- endemic_eq$IA / endemic_eq$IS 


# Only seems to converge where input is near final values - use long-run output
endemic_eq <- print(tail(out, 1))
endemic_eq$IA / endemic_eq$IS 

# Estimate impact of changing pS probability on asymptomatic ratio
# Set pS probability to run through
pSs <- seq(0.01, 1, by = 0.02)

# Run through probabilities and calculate equilibrium asymptomatic ratio
as_rats <- list()
for(i in 1:length(pSs)){
  paras["pS"] <- pSs[i]
  # Solve using nleqslv
  # sol <- nleqslv(x0, equilibrium_system, paras = paras)
  # endemic_eq <- data.frame(
  #   S1 = sol$x[1],
  #   S2 = sol$x[2],
  #   IS = sol$x[3],
  #   IA = sol$x[4],
  #   R1 = sol$x[5],
  #   R2 = sol$x[6],
  #   W  = sol$x[7]
  # )
  out <- ode(y = start, times = times, func = asmod, parms = paras)
  endemic_eq <- tail(out, 1)
  # Calculate equilibrium asymptomatic ratio
  as_rats[i] <- endemic_eq[,'IA'] / endemic_eq[,'IS']
  print(i)
}

# Plot endemic ratio against initial symptomatic ratio
plot(x = pSs, y = as_rats, type = 'l', xlab = "Probability of being symptomatic upon first infection", ylab = "Ratio of asymptomatic cases at equilibrium")

# Linear interpolation of pS for ratio of 820 (Hedge et al., 2024)
# Target ratio
ratio_target <- 820

# Interpolate x given y
result <- approx(x = as_rats, y = pSs, xout = ratio_target)
pS_interp <- result$y
pS_interp

# Impact of transmission heterogeneity ------------------------------
# Set times to run simulation  - 20 should be sufficient to calculate equilibrium
times = seq(0, 365*20, by = 1)


# Infectiousness of symptomatic and asymptomatic individuals
# Reset baseline parameters
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Vary infectiousness of symptomatic and asymptomatic individuals
betaDSs <- seq(0, 0.4, length.out = 50)
betaDAs <- seq(0, 0.4, length.out = 50)
  
# Set empty output arrays
betaR0s <- array(NA, dim = c(length(betaDSs), length(betaDAs)))
betaEq_infs <- array(NA, dim = c(length(betaDSs), length(betaDAs)))
betaEq_rats <- array(NA, dim = c(length(betaDSs), length(betaDAs)))
             
# Run through parameters and calculate outputs
for(i in 1:length(betaDSs)){
  paras["betaDS"] <- betaDSs[i]
  for(j in 1:length(betaDAs)){
    paras["betaDA"] <- betaDAs[j]
    
    # Calculate R0
    betaR0s[i, j] <- calc_R0(paras)
    
    # Estimate equilibrium number of infecteds - use end ODE output as more reliable
    out <- ode(y = start, times = times, func = asmod, parms = paras)
    equil <- tail(out, 1)
    betaEq_infs[i, j] <- equil[, 'IS']
    
    # Estimate equilibrium asymptomatic ratio
    betaEq_rats[i, j] <- equil[, 'IA'] / equil[, 'IS']
    
    print(i)
  }
}

# Plot results
# Wrangle to data frame of R0s
R0_df <- melt(betaR0s)
names(R0_df) <- c("betaDS", "betaDA", "R0")
# Map to actual x/y variable values
R0_df$betaDS <- betaDSs[R0_df$betaDS]
R0_df$betaDA <- betaDAs[R0_df$betaDA]

# Plot heat map for R0s
ggplot(R0_df, aes(x = betaDS, y = betaDA, fill = R0, z = R0)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(
    x = "Symptomatic transmission rate",
    y = "Asymptomatic transmission rate",
    fill = "R0"
  ) +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        legend.title = element_text(size = 16))

# Wrangle to data frame of Eq_infs
Eq_inf_df <- melt(betaEq_infs)
names(Eq_inf_df) <- c("betaDS", "betaDA", "Eq_inf")
# Map to actual x/y variable values
Eq_inf_df$betaDS <- betaDSs[Eq_inf_df$betaDS]
Eq_inf_df$betaDA <- betaDAs[Eq_inf_df$betaDA]

# Plot heat map for Eq_infs 
ggplot(Eq_inf_df, aes(x = betaDS, y = betaDA, fill = Eq_inf)) +
  geom_tile() +
  labs(x = "Symptomatic transmission rate",
       y = "Asymptomatic transmission rate",
       fill = "Eq_inf") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        legend.title = element_text(size = 16))

# Wrangle to data frame of Eq_rats
Eq_rat_df <- melt(betaEq_rats)
names(Eq_rat_df) <- c("betaDS", "betaDA", "Eq_rat")
# Map to actual x/y variable values
Eq_rat_df$betaDS <- betaDSs[Eq_rat_df$betaDS]
Eq_rat_df$betaDA <- betaDAs[Eq_rat_df$betaDA]

# Plot heat map for Eq_rats 
ggplot(Eq_rat_df, aes(x = betaDS, y = betaDA, fill = Eq_rat, z = Eq_rat)) +
  geom_tile() +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(x = "Symptomatic transmission rate",
       y = "Asymptomatic transmission rate",
       fill = "Eq_rat") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        legend.title = element_text(size = 16))

# Duration of infection of symptomatic and asymptomatic individuals
# Reset baseline parameters
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Vary duration of infection of symptomatic and asymptomatic individuals - set as duration of infection (days)
durSs <- seq(2, 20, length.out = 100)
durAs <- seq(2, 20, length.out = 100)

# Set empty output arrays
gammaR0s <- array(NA, dim = c(length(durSs), length(durAs)))
gammaEq_infs <- array(NA, dim = c(length(durSs), length(durAs)))
gammaEq_rats <- array(NA, dim = c(length(durSs), length(durAs)))

# Run through parameters and calculate outputs
for(i in 1:length(durSs)){
  paras["gammaS"] <- 1 / durSs[i] # Set gamma as 1/duration
  
  for(j in 1:length(durAs)){
    paras["gammaA"] <- 1 / durAs[j]
    
    # Calculate R0
    gammaR0s[i, j] <- calc_R0(paras)
    
    # Estimate equilibrium number of infecteds - use end ODE output as more reliable
    out <- ode(y = start, times = times, func = asmod, parms = paras)
    equil <- tail(out, 1)
    gammaEq_infs[i, j] <- equil[, 'IS']
    
    # Estimate equilibrium asymptomatic ratio
    gammaEq_rats[i, j] <- equil[, 'IA'] / equil[, 'IS']
  }
}

# Plot results
# Wrangle to data frame of R0s
R0_df <- melt(gammaR0s)
names(R0_df) <- c("durS", "durA", "R0")
# Map to actual x/y variable values
R0_df$durS <- durSs[R0_df$durS]
R0_df$durA <- durAs[R0_df$durA]

# Plot heat map for R0s
ggplot(R0_df, aes(x = durS, y = durA, fill = R0, z = R0)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(
    x = "Duration of symptomatic infection (days)",
    y = "Duration of asymptomatic infection (days)",
    fill = "R0"
  ) +
  theme_minimal()

# Wrangle to data frame of Eq_infs
Eq_inf_df <- melt(gammaEq_infs)
names(Eq_inf_df) <- c("durS", "durA", "Eq_inf")
# Map to actual x/y variable values
Eq_inf_df$durS <- durSs[Eq_inf_df$durS]
Eq_inf_df$durA <- durAs[Eq_inf_df$durA]

# Plot heat map for Eq_infs 
ggplot(Eq_inf_df, aes(x = durS, y = durA, fill = Eq_inf)) +
  geom_tile() +
  labs(x = "Duration of symptomatic infection (days)",
       y = "Duration of asymptomatic infection (days)",
       fill = "Eq_inf") +
  theme_minimal()

# Wrangle to data frame of Eq_rats
Eq_rat_df <- melt(gammaEq_rats)
names(Eq_rat_df) <- c("durS", "durA", "Eq_rat")
# Map to actual x/y variable values
Eq_rat_df$durS <- durSs[Eq_rat_df$durS]
Eq_rat_df$durA <- durAs[Eq_rat_df$durA]

# Plot heat map for Eq_rats 
ggplot(Eq_rat_df, aes(x = durS, y = durA, fill = Eq_rat, z = Eq_rat)) +
  geom_tile() +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(x = "Duration of symptomatic infection (days)",
       y = "Duration of asymptomatic infection (days)",
       fill = "Eq_rat") +
  theme_minimal()


# Shedding rates of symptomatic and asymptomatic individuals
# Reset baseline parameters
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Vary shedding rates of symptomatic and asymptomatic individuals
alphaSs <- seq(0, 100, by = 1)
alphaAs <- seq(0, 100, by = 1)

# Set empty output arrays
alphaR0s     <- array(NA, dim = c(length(alphaSs), length(alphaAs)))
alphaEq_infs <- array(NA, dim = c(length(alphaSs), length(alphaAs)))
alphaEq_rats <- array(NA, dim = c(length(alphaSs), length(alphaAs)))

# Run through parameters and calculate outputs
for(i in 1:length(alphaSs)){
  paras["alphaS"] <- alphaSs[i]
  
  for(j in 1:length(alphaAs)){
    paras["alphaA"] <- alphaAs[j]
    
    # Calculate R0
    alphaR0s[i, j] <- calc_R0(paras)
    
    # Estimate equilibrium number of infecteds - use end ODE output as more reliable
    out <- ode(y = start, times = times, func = asmod, parms = paras)
    equil <- tail(out, 1)
    alphaEq_infs[i, j] <- equil[, 'IS']
    
    # Estimate equilibrium asymptomatic ratio
    alphaEq_rats[i, j] <- equil[, 'IA'] / equil[, 'IS']
  }
}

# Plot results
# Wrangle to data frame of R0s
R0_df <- melt(alphaR0s)
names(R0_df) <- c("alphaS", "alphaA", "R0")
# Map to actual x/y variable values
R0_df$alphaS <- alphaSs[R0_df$alphaS]
R0_df$alphaA <- alphaAs[R0_df$alphaA]

# Plot heat map for R0s
ggplot(R0_df, aes(x = alphaS, y = alphaA, fill = R0, z = R0)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(
    x = "Symptomatic shedding rate",
    y = "Asymptomatic shedding rate",
    fill = "R0"
  ) +
  theme_minimal()

# Wrangle to data frame of Eq_infs
Eq_inf_df <- melt(alphaEq_infs)
names(Eq_inf_df) <- c("alphaS", "alphaA", "Eq_inf")
# Map to actual x/y variable values
Eq_inf_df$alphaS <- alphaSs[Eq_inf_df$alphaS]
Eq_inf_df$alphaA <- alphaAs[Eq_inf_df$alphaA]

# Plot heat map for Eq_infs 
ggplot(Eq_inf_df, aes(x = alphaS, y = alphaA, fill = Eq_inf)) +
  geom_tile() +
  labs(x = "Symptomatic shedding rate",
       y = "Asymptomatic shedding rate",
       fill = "Eq_inf") +
  theme_minimal()

# Wrangle to data frame of Eq_rats
Eq_rat_df <- melt(alphaEq_rats)
names(Eq_rat_df) <- c("alphaS", "alphaA", "Eq_rat")
# Map to actual x/y variable values
Eq_rat_df$alphaS <- alphaSs[Eq_rat_df$alphaS]
Eq_rat_df$alphaA <- alphaAs[Eq_rat_df$alphaA]

# Plot heat map for Eq_rats 
ggplot(Eq_rat_df, aes(x = alphaS, y = alphaA, fill = Eq_rat, z = Eq_rat)) +
  geom_tile() +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(x = "Symptomatic shedding rate",
       y = "Asymptomatic shedding rate",
       fill = "Eq_rat") +
  theme_minimal()


# Loss of immunity of symptomatic and asymptomatic individuals
# Reset baseline parameters
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Vary loss of immunity of symptomatic and asymptomatic individuals - duration of immunity (days)
sigma1s <- seq(30, 365*4, length.out = 100)
sigma2s <- seq(30, 365*4, length.out = 100)

# Set empty output arrays
sigmaR0s <- array(NA, dim = c(length(sigma1s), length(sigma2s)))
sigmaEq_infs <- array(NA, dim = c(length(sigma1s), length(sigma2s)))
sigmaEq_rats <- array(NA, dim = c(length(sigma1s), length(sigma2s)))

# Run through parameters and calculate outputs
for(i in 1:length(sigma1s)){
  paras["sigma1"] <- 1 / sigma1s[i]
  
  for(j in 1:length(sigma2s)){
    paras["sigma2"] <- 1 / sigma2s[j]
    
    # Calculate R0
    sigmaR0s[i, j] <- calc_R0(paras)
    
    # Estimate equilibrium number of infecteds - use end ODE output as more reliable
    out <- ode(y = start, times = times, func = asmod, parms = paras)
    equil <- tail(out, 1)
    sigmaEq_infs[i, j] <- equil[, 'IS']
    
    # Estimate equilibrium asymptomatic ratio
    sigmaEq_rats[i, j] <- equil[, 'IA'] / equil[, 'IS']
  }
}

# Plot results
# Wrangle to data frame of R0s
R0_df <- melt(sigmaR0s)
names(R0_df) <- c("sigma1", "sigma2", "R0")
# Map to actual x/y variable values
R0_df$sigma1 <- sigma1s[R0_df$sigma1]
R0_df$sigma2 <- sigma2s[R0_df$sigma2]

# Plot heat map for R0s
ggplot(R0_df, aes(x = sigma1, y = sigma2, fill = R0, z = R0)) +
  geom_tile(aes(fill = R0)) +
  geom_contour(color = "black", linewidth = 0.5) +
  labs(x = "Duration of immunity from primary infection/vaccination",
       y = "Duration of immunity from secondary infection",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of Eq_infs
Eq_inf_df <- melt(sigmaEq_infs)
names(Eq_inf_df) <- c("sigma1", "sigma2", "Eq_inf")
# Map to actual x/y variable values
Eq_inf_df$sigma1 <- sigma1s[Eq_inf_df$sigma1]
Eq_inf_df$sigma2 <- sigma2s[Eq_inf_df$sigma2]

# Plot heat map for Eq_infs 
ggplot(Eq_inf_df, aes(x = sigma1, y = sigma2, fill = Eq_inf)) +
  geom_tile() +
  labs(x = "Duration of immunity from primary infection/vaccination",
       y = "Duration of immunity from secondary infection",
       fill = "Eq_inf") +
  theme_minimal()

# Wrangle to data frame of Eq_rats
Eq_rat_df <- melt(sigmaEq_rats)
names(Eq_rat_df) <- c("sigma1", "sigma2", "Eq_rat")
# Map to actual x/y variable values
Eq_rat_df$sigma1 <- sigma1s[Eq_rat_df$sigma1]
Eq_rat_df$sigma2 <- sigma2s[Eq_rat_df$sigma2]

# Plot heat map for Eq_rats 
ggplot(Eq_rat_df, aes(x = sigma1, y = sigma2, fill = Eq_rat)) +
  geom_tile() +
  labs(x = "Duration of immunity from primary infection/vaccination",
       y = "Duration of immunity from secondary infection",
       fill = "Eq_rat") +
  theme_minimal()


# Impact of transmission heterogeneity with intervention ------------------
# Set baseline impact of vaccination on outputs
# Reset baseline parameters
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Set levels of vaccination
vacc_levels <- seq(0, 0.001, by = 0.00001)

# Set empty output arrays
vaccR0s     <- array(NA, dim = c(length(vacc_levels)))
vaccEq_infs <- array(NA, dim = c(length(vacc_levels)))
vaccEq_rats <- array(NA, dim = c(length(vacc_levels)))

# Run through vaccination levels and calculate outputs
for(i in 1:length(vacc_levels)){
  paras["u"] <- vacc_levels[i]
  
  # Calculate R0
  vaccR0s[i] <- calc_R0(paras)
  
  # Estimate equilibrium number of infecteds - use end ODE output as more reliable
  out <- ode(y = start, times = times, func = asmod, parms = paras)
  equil <- tail(out, 1)
  vaccEq_infs[i] <- equil[, 'IS']
  
  # Estimate equilibrium asymptomatic ratio
  vaccEq_rats[i] <- equil[, 'IA'] / equil[, 'IS']
}

# Plot R0 against vaccinated rate
plot(x = vacc_levels, y = vaccR0s, type = 'l',
     xlab = "Vaccination level per day", ylab = "R0")

# Plot equilibrium infecteds against vaccinated rate
plot(x = vacc_levels, y = vaccEq_infs, type = 'l',
     xlab = "Vaccination level per day", ylab = "Proportion infected at equilibrium")

# Plot equilibrium ratio against vaccinated rate
plot(x = vacc_levels, y = vaccEq_rats, type = 'l',
     xlab = "Vaccination level per day", ylab = "Asymptomatic ratio at equilibrium")


# Repeat under different conditions
# Change betaDs - relative transmission rates of asymptomatic and symptomatic individuals
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW  = 0.02,        # Indirect transmission rate
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,  # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Run through proportions of primary infections that are symptomatic
betaDSs <- seq(0, 0.4, by = 0.04)

# Create empty array to store outputs
vaccR0s <- array(NA, dim = c(length(betaDSs), length(vacc_levels)))

# Run through vaccination levels and calculate outputs
for(i in 1:length(betaDSs)){
  paras["betaDS"] <- betaDSs[i]
  paras["betaDA"] <- 0.4 - betaDSs[i]
  
  for(j in 1: length(vacc_levels)){
    paras["u"] <- vacc_levels[j]
  
    # Calculate R0
    vaccR0s[i,j] <- calc_R0(paras)
    
    # Estimate equilibrium number of infecteds - use end ODE output as more reliable
    out <- ode(y = start, times = times, func = asmod, parms = paras)
    equil <- tail(out, 1)
    vaccEq_infs[i] <- equil[, 'IS']
    
    # Estimate equilibrium asymptomatic ratio
    vaccEq_rats[i] <- equil[, 'IA'] / equil[, 'IS']
    print(i)
  }
}

# Wrangle to data frame of R0s under vaccination
vaccR0_df <- melt(vaccR0s)
names(vaccR0_df) <- c("betaDS", "vacc_level", "R0")
# Map to actual x/y variable values
vaccR0_df$betaDS <- betaDSs[vaccR0_df$betaDS]
vaccR0_df$betaDS <- (0.4 - vaccR0_df$betaDS) * 100 / 0.4
vaccR0_df$vacc_level <- vacc_levels[vaccR0_df$vacc_level]
vaccR0_df$betaDS <- as.factor(vaccR0_df$betaDS)

# Plot R0 against vaccinated rate
# Create a red-to-blue palette with as many colors as levels in betaDS
cols <- colorRampPalette(c("blue", "red"))(length(levels(vaccR0_df$betaDS)))

ggplot(data = vaccR0_df, aes(x = vacc_level, y = R0, color = betaDS)) +
  geom_line() +
  labs(
    x = "Vaccination level per day",
    col = "Symptomatic transmission rate \nas a proportion of total \ntransmission"
  ) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(axis.title = element_text(size = 16))


# Seasonality -------------------------------------------------------------
# Introduce seasonality on betaW

# Create function to run asymptomatic seasonal model
assmod <- function(t, y, parameters){
  # Pull state variables from y vector
  S1 = y[1] # Susceptible class without prior infection or vaccination
  IS = y[2] # Symptomatic infectious class
  R1 = y[3] # Recovered class following primary infection or vaccination
  S2 = y[4] # Susceptible class following previous infection or vaccination
  IA = y[5] # Asymptomatic/mild infectious class
  R2 = y[6] # Recovered class following asymptomatic infection
  W  = y[7] # Bacterial concentration in the water
  
  # Pull parameter values from the input vector 
  mu     = parameters["mu"]     # Birth and death rate of each compartment
  N      = parameters["N"]      # Total population size
  betaW0 = parameters["betaW0"] # Baseline indirect transmission rate
  A      = parameters["A"]      # Seasonal amplitude
  phi    = parameters["phi"]    # Peak timing
  T      = parameters["T"]      # Seasonality period
  betaDS = parameters["betaDS"] # Direct transmission rate from symptomatic individuals
  betaDA = parameters["betaDA"] # Direct transmission rate from asymptomatic individuals
  u      = parameters["u"]      # Vaccination rate
  pS     = parameters["pS"]     # Proportion of initial infections which are symptomatic
  gammaS = parameters["gammaS"] # Recovery rate for symptomatic infections
  gammaA = parameters["gammaA"] # Recovery rate for symptomatic infections
  sigma1 = parameters["sigma1"] # Rate of loss of immunity for recovered class 1
  sigma2 = parameters["sigma2"] # Rate of loss of immunity for recovered class 2
  alphaS = parameters["alphaS"] # Shedding rate from sympomatic individuals
  alphaA = parameters["alphaA"] # Shedding rate from asympomatic individuals
  muW    = parameters["muW"]    # Rate of growth/loss of bacteria in the water
  
  # Seasonal betaW
  betaW_t <- betaW0 * (1 + A * sin(2 * pi * (t - phi) / T))
  
  # Define equations
  dS1 = mu * N - (betaW_t * S1 * W + betaDS * S1 * IS + betaDA * S1 * IA) / N - (u + mu) * S1
  dIS = pS * (betaW_t * S1 * W + betaDS * S1 * IS + betaDA * S1 * IA) / N - (gammaS + mu) * IS
  dR1 = u * S1 + gammaS * IS - (sigma1 + mu) * R1
  dS2 = sigma1 * R1 + sigma2 * R2 - (betaW_t * S2 * W + betaDS * S2 * IS + betaDA * S2 * IA) / N - mu * S2
  dIA = (1 - pS) * (betaW_t * S1 * W + betaDS * S1 * IS + betaDA * S1 * IA) / N + (betaW_t * S2 * W + betaDS * S2 * IS + betaDA * S2 * IA) / N - (gammaA + mu) * IA
  dR2 = gammaA * IA - (sigma2 + mu) * R2
  dW  = alphaS * IS + alphaA * IA + muW * W
  res = c(dS1, dIS, dR1, dS2, dIA, dR2, dW)
  
  # Return list of gradients
  list(res)
}

# Set times to run simulation for
times = seq(0, 365*10, by = 1)

# set parameters for initial model
paras = c(mu     = 1/(45*365), # Birth and death rate of each compartment
          N      = 10000,      # Total population size
          betaW0 = 0.02,       # Indirect transmission rate
          A      = 1,          # Seasonal amplitude
          phi    = 100,        # Peak timing
          T      = 182.5,      # Seasonality period
          betaDS = 0.2,        # Direct transmission rate from symptomatic individuals
          betaDA = 0.2,        # Direct transmission rate from asymptomatic individuals
          u      = 0,          # Vaccination rate
          pS     = pS_interp,        # Proportion of initial infections which are symptomatic
          gammaS = 1/5,        # Recovery rate for symptomatic infections
          gammaA = 1/10,       # Recovery rate for asymptomatic infections
          sigma1 = 1/(3*365),  # Rate of loss of immunity for recovered class 1
          sigma2 = 1/50,       # Rate of loss of immunity for recovered class 2
          alphaS = 100,         # Shedding rate from symptomatic individuals
          alphaA = 10,         # Shedding rate from asymptomatic individuals
          muW    = -10)        # Rate of growth/loss of bacteria in the water

# Set population size
N <- 10000

# Set initial conditions - epidemic conditions
start <- c(S1 = 0.99 * N,  # Susceptible class without prior infection or vaccination
           IS = 0.005 * N, # Symptomatic infectious class
           R1 = 0, # Recovered class following primary infection or vaccination
           S2 = 0, # Susceptible class following previous infection or vaccination
           IA = 0.005 * N, # Asymptomatic/mild infectious class
           R2 = 0, # Recovered class following asymptomatic infection
           W  = 0) # Bacterial concentration in the water

# Simulate dynamics
out <- ode(y = start, times = times, func = assmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$S1, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0, 10000), cex.lab = 1.5)
lines(x = out$time, y = out$IS, col = "red") 
lines(x = out$time, y = out$R1, col = "blue")
lines(x = out$time, y = out$S2, col = "grey")
lines(x = out$time, y = out$IA, col = "pink") 
lines(x = out$time, y = out$R2, col = "lightblue")
legend(x = 0.8 * max(out$time), y = 0.9 * N,                                
       legend = c("S1", "IS", "R1", "S2", "IA", "R2"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend



# Plot reservoir
plot(x = out$time, y = log(out$W), col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "log(Bacterial concentration in water)")

# Check constant population size
out$Ntot <- rowSums(out[, c("S1", "IS", "R1", "S2", "IA", "R2")])
plot(out$time, out$Ntot, type='l')


# Now calculate Re with seasonality
calc_Re_seasonal_breakdown <- function(out, paras) {
  with(as.list(paras), {
    s <- mu / (u + mu)
    
    # Pre-allocate vectors
    Re_S1 <- numeric(length(out$time))
    Re_S2 <- numeric(length(out$time))
    
    for (i in seq_along(out$time)) {
      t <- out$time[i]
      
      # Seasonal betaW
      betaW_t <- betaW0 * (1 + A * sin(2*pi*(t - phi)/T))
      
      # Contributions per susceptible
      direct_S1 <- pS * betaDS / (gammaS + mu) + (1 - pS) * betaDA / (gammaA + mu)
      env_S1    <- (betaW_t / (-muW)) * (pS * alphaS / (gammaS + mu) + (1 - pS) * alphaA / (gammaA + mu))
      Re_S1[i]  <- (out$S1[i] / N) * (direct_S1 + env_S1) / s
      
      direct_S2 <- betaDS / (gammaS + mu) + betaDA / (gammaA + mu)
      env_S2    <- (betaW_t / (-muW)) * (alphaS / (gammaS + mu) + alphaA / (gammaA + mu))
      Re_S2[i]  <- (out$S2[i] / N) * (direct_S2 + env_S2) / s
    }
    
    data.frame(time = out$time,
               Re_S1 = Re_S1,
               Re_S2 = Re_S2,
               Re_total = Re_S1 + Re_S2)
  })
}


# Compute seasonal Re breakdown
Re_seasonal_break <- calc_Re_seasonal_breakdown(out, paras)

# Plot
plot(Re_seasonal_break$time, Re_seasonal_break$Re_total, type="l", col="black", lwd=2,
     xlab="Time (days)", ylab = expression(R[e](t)),
     main="Seasonal R_e(t) decomposed by susceptible class",
     ylim = c(0,2))
lines(Re_seasonal_break$time, Re_seasonal_break$Re_S1, col="blue", lwd=2)
lines(Re_seasonal_break$time, Re_seasonal_break$Re_S2, col="green", lwd=2)
abline(h=1, lty=2, col="red")
legend("topright",
       legend=c("Total R_e","Primary S1","Secondary S2"),
       col=c("black","blue","green"), lwd=2, bty="n")




