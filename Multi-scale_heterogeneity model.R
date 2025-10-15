# Multi-scale heterogeneity model
# Oliver G. Spacey
# 2025.10.14

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load required packages
library(deSolve)   # differential equations
library(ggplot2)   # data visualisation
library(reshape2)  # data wrangling
library(scales)    # data wrangling
library(tidyverse) # data wrangling

# Initial Model -----------------------------------------------------------
# Create function to run multi-scale model
msmod <- function(t, y, parameters){
  # Pull state variables from y vector
  SA  = y[1] # Susceptibles in patch A
  IAH = y[2] # High-infectiousness infectious individuals in patch A
  IAL = y[3] # Low-infectiousness infectious individuals in patch A
  RA  = y[4] # Recovered-and-immune individuals in patch A
  SB  = y[5] # Susceptibles in patch B
  IBH = y[6] # High-infectiousness infectious individuals in patch B
  IBL = y[7] # Low-infectiousness infectious individuals in patch B
  RB  = y[8] # Recovered-and-immune individuals in patch B
  W   = y[9] # Bacterial concentration in the environmental reservoir
  
  # Pull parameter values from the input vector 
  mu      = parameters["mu"]      # Birth and death rate of each compartment
  betaWA  = parameters["betaWA"]  # Indirect transmission rate to patch A
  betaWB  = parameters["betaWB"]  # Indirect transmission rate to patch B
  kappa   = parameters["kappa"]   # Saturating coefficient - infectious dose sufficient for 50% infection   
  pAH     = parameters["pAH"]     # Proportion of high-infectiousness individuals in patch A
  pAL     = parameters["pAL"]           # Proportion of low-infectiousness individuals in patch A
  pBH     = parameters["pBH"]     # Proportion of high-infectiousness individuals in patch B
  pBL     = parameters["pBL"]          # Proportion of low-infectiousness individuals in patch B
  betaDAH = parameters["betaDAH"] # Direct transmission rate for high-infectiousness individuals in patch A
  betaDAL = parameters["betaDAL"] # Direct transmission rate for low-infectiousness individuals in patch A
  betaDBH = parameters["betaDBH"] # Direct transmission rate for high-infectiousness individuals in patch B
  betaDBL = parameters["betaDBL"] # Direct transmission rate for low-infectiousness individuals in patch B
  zeta    = parameters["zeta"]    # Inter-patch mixing coefficient
  gamma   = parameters["gamma"]   # Recovery rate
  aAH     = parameters["aAH"]     # Shedding rate from high-infectiousness individuals in patch A
  aAL     = parameters["aAL"]     # Shedding rate from low-infectiousness individuals in patch A
  aBH     = parameters["aBH"]     # Shedding rate from high-infectiousness individuals in patch B
  aBL     = parameters["aBL"]     # Shedding rate from low-infectiousness individuals in patch B
  xi      = parameters["xi"]      # Rate of loss of bacteria from the environmental reservoir
  
  # Define equations
  dSA  = mu - mu * SA - betaWA * SA * (W / (W + kappa)) - (betaDAH * IAH + betaDAL * IAL) * SA - zeta * (betaDBH * IBH + betaDBL * IBL) * SA
  dIAH = pAH * (betaWA * SA * (W / (W + kappa)) + (betaDAH * IAH + betaDAL * IAL) * SA + zeta * (betaDBH * IBH + betaDBL * IBL) * SA) - (mu + gamma) * IAH
  dIAL = pAL * (betaWA * SA * (W / (W + kappa)) + (betaDAH * IAH + betaDAL * IAL) * SA + zeta * (betaDBH * IBH + betaDBL * IBL) * SA) - (mu + gamma) * IAL
  dRA  = gamma * (IAH + IAL) - mu * RA
  dSB  = mu - mu * SB - betaWB * SB * (W / (W + kappa)) - (betaDBH * IBH + betaDBL * IBL) * SB - zeta * (betaDAH * IAH + betaDAL * IAL) * SB
  dIBH = pBH * (betaWB * SB * (W / (W + kappa)) + (betaDBH * IBH + betaDBL * IBL) * SB + zeta * (betaDAH * IAH + betaDAL * IAL) * SB) - (mu + gamma) * IBH
  dIBL = pBL * (betaWB * SB * (W / (W + kappa)) + (betaDBH * IBH + betaDBL * IBL) * SB + zeta * (betaDAH * IAH + betaDAL * IAL) * SB) - (mu + gamma) * IBL
  dRB  = gamma * (IBH + IBL) - mu * RB
  dW   = xi * (aAH * IAH + aAL * IAL + aBH * IBH + aBL * IBL - W)
  res  = c(dSA, dIAH, dIAL, dRA, dSB, dIBH, dIBL, dRB, dW)
  
  # Return list of gradients
  list(res)
}

# Set times to run simulation for
times = seq(0, 365*1000, by = 1)

# Set population sizes and original parameters inspired by Mukandivire et al.
N_A <- 10000 # Patch A population size
N_B <- 10000 # Patch B population size
pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
pAL <- 0.5
pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
pBL <- 0.5
b_D_AH <- 0.00003 # Direct transmission rate for high-infectiousness individuals in patch A
b_D_AL <- 0.00003 # Direct transmission rate for low-infectiousness individuals in patch A
b_D_BH <- 0.00003 # Direct transmission rate for high-infectiousness individuals in patch B
b_D_BL <- 0.00003 # Direct transmission rate for low-infectiousness individuals in patch B
b_W_A <- 0.1 # Indirect transmission rate for patch A
b_W_B <- 0.1 # Indirect transmission rate for patch B
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 10 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 10 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A * pAH + alpha_AL * N_A * pAL + alpha_BH * N_B * pBH + alpha_BL * N_B * pBL # Total shedding from all individuals, for scaling
xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
  
# Calculate non-dimensionalised parameters
# Initial homogeneous model where ΣaβW/κ = ΣβD 
paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
          betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
          betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
          pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = 0.5,
          pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = 0.5,
          betaDAH = b_D_AH * N_A * pAH,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A * pAL,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B * pBH,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B * pBL,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.05, # Inter-patch mixing is 5% strength
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A * pAH / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A * pAL / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B * pBH / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B * pBL / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate total indirect transmission - ΣaβW/κ 
tot_w_trnm <- as.numeric((paras['aAH'] * paras['betaWA'] + paras['aAL'] * paras['betaWA'] + paras['aBH'] * paras['betaWB'] + paras['aBL'] * paras['betaWB']) / paras['kappa'])

# Calculate total indirect transmission - ΣβD 
tot_d_trnm <- as.numeric(paras['betaDAH'] + paras['betaDAL'] + paras['betaDBH'] + paras['betaDBL'])

# Set initial conditions - epidemic conditions
start <- c(SA  = 0.98,  # Susceptibles in patch A
           IAH = 0.01, # High-infectiousness infectious individuals in patch A
           IAL = 0.01, # Low-infectiousness infectious individuals in patch A
           RA  = 0,     # Recovered-and-immune individuals in patch A
           SB  = 0.98,  # Susceptibles in patch B
           IBH = 0.01, # High-infectiousness infectious individuals in patch B
           IBL = 0.01, # Low-infectiousness infectious individuals in patch B
           RB  = 0,     # Recovered-and-immune individuals in patch B
           W   = 0)  # Bacterial concentration in the environmental reservoir

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top of data frame
head(round(out, 3))

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")


# Basic Reproduction Number ------------------------------------------------------------

# Set parameters as previous parameters and compartments at disease-free equilibrium to evaluate R0
dfeparameters <- c(SA  = 1, # Susceptibles in patch A
               IAH = 0, # High-infectiousness infectious individuals in patch A
               IAL = 0, # Low-infectiousness infectious individuals in patch A
               RA  = 0, # Recovered-and-immune individuals in patch A
               SB  = 1, # Susceptibles in patch B
               IBH = 0, # High-infectiousness infectious individuals in patch B
               IBL = 0, # Low-infectiousness infectious individuals in patch B
               RB  = 0, # Recovered-and-immune individuals in patch B
               W   = 0, # Bacterial concentration in the environmental reservoir
               paras)   # Add other parameters

# Convert dfe parameters to list
dfeparameters <- as.list(dfeparameters)

get.R0 <- function(dfeparas){
  # Define F matrix - new infections to infectious compartments
  F1 = quote(pAH * betaWA * SA * (W / (W + kappa)) + pAH * betaDAH * SA * IAH + pAH * betaDAL * IAL * SA + zeta * pAH * betaDBH * SA * IBH + zeta * pAH * betaDBL * SA * IBL)
  F2 = quote(pAL * betaWA * SA * (W / (W + kappa)) + pAL * betaDAH * SA * IAH + pAL * betaDAL * IAL * SA + zeta * pAL * betaDBH * SA * IBH + zeta * pAL * betaDBL * SA * IBL)
  F3 = quote(pBH * betaWB * SB * (W / (W + kappa)) + pBH * betaDBH * SB * IBH + pBH * betaDBL * IBL * SB + zeta * pBH * betaDAH * SB * IAH + zeta * pBH * betaDAL * SB * IAL)
  F4 = quote(pBL * betaWB * SB * (W / (W + kappa)) + pBL * betaDBH * SB * IBH + pBL * betaDBL * IBL * SB + zeta * pBL * betaDAH * SB * IAH + zeta * pBL * betaDAL * SB * IAL)
  F5 = 0
  
  # Define V matrix - losses and transfers from infectious compartments
  V1 = quote((gamma + mu) * IAH)
  V2 = quote((gamma + mu) * IAL)
  V3 = quote((gamma + mu) * IBH)
  V4 = quote((gamma + mu) * IBL)
  V5 = quote(xi * W - xi * aAH * IAH - xi * aAL * IAL - xi * aBH * IBH - xi * aBL * IBL)
  
  # Generate partial derivaties for the Jacobian matrices
  f11 = D(F1, "IAH")
  f12 = D(F1, "IAL")
  f13 = D(F1, "IBH")
  f14 = D(F1, "IBL")
  f15 = D(F1, "W")
  f21 = D(F2, "IAH")
  f22 = D(F2, "IAL")
  f23 = D(F2, "IBH")
  f24 = D(F2, "IBL")
  f25 = D(F2, "W")
  f31 = D(F3, "IAH")
  f32 = D(F3, "IAL")
  f33 = D(F3, "IBH")
  f34 = D(F3, "IBL")
  f35 = D(F3, "W")
  f41 = D(F4, "IAH")
  f42 = D(F4, "IAL")
  f43 = D(F4, "IBH")
  f44 = D(F4, "IBL")
  f45 = D(F4, "W")
  f51 = D(F5, "IAH")
  f52 = D(F5, "IAL")
  f53 = D(F5, "IBH")
  f54 = D(F5, "IBL")
  f55 = D(F5, "W")

  v11 = D(V1, "IAH")
  v12 = D(V1, "IAL")
  v13 = D(V1, "IBH")
  v14 = D(V1, "IBL")
  v15 = D(V1, "W")
  v21 = D(V2, "IAH")
  v22 = D(V2, "IAL")
  v23 = D(V2, "IBH")
  v24 = D(V2, "IBL")
  v25 = D(V2, "W")
  v31 = D(V3, "IAH")
  v32 = D(V3, "IAL")
  v33 = D(V3, "IBH")
  v34 = D(V3, "IBL")
  v35 = D(V3, "W")
  v41 = D(V4, "IAH")
  v42 = D(V4, "IAL")
  v43 = D(V4, "IBH")
  v44 = D(V4, "IBL")
  v45 = D(V4, "W")
  v51 = D(V5, "IAH")
  v52 = D(V5, "IAL")
  v53 = D(V5, "IBH")
  v54 = D(V5, "IBL")
  v55 = D(V5, "W")
  
  # Construct matrices
  f = with(dfeparas,
           matrix(c(eval(f11), eval(f12), eval(f13), eval(f14), eval(f15),
                    eval(f21), eval(f22), eval(f23), eval(f24), eval(f25),
                    eval(f31), eval(f32), eval(f33), eval(f34), eval(f35),
                    eval(f41), eval(f42), eval(f43), eval(f44), eval(f45),
                    eval(f51), eval(f52), eval(f53), eval(f54), eval(f55)),
                  nrow = 5, byrow = TRUE))
  v = with(dfeparas,
           matrix(c(eval(v11), eval(v12), eval(v13), eval(v14), eval(v15),
                    eval(v21), eval(v22), eval(v23), eval(v24), eval(v25),
                    eval(v31), eval(v32), eval(v33), eval(v34), eval(v35),
                    eval(v41), eval(v42), eval(v43), eval(v44), eval(v45),
                    eval(v51), eval(v52), eval(v53), eval(v54), eval(v55)),
                  nrow = 5, byrow = TRUE))
  
  # Calculate R0 = ρ(FV^(-1))
  return(max(Mod(eigen(f %*% solve(v))$values)))
}

# Get R0 for initial model
get.R0(dfeparameters)

# Long-term Dynamics ----------------------------------------------

# Run initial model for 50 years
times = seq(0, 365*50, by = 1)

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top of data frame
head(round(out, 3))

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")

# Patch heterogeneity only ------------------------------------------------
# Vary patch heterogeneity in indirect and direct transmission
# Analogy - two sites differ in their exposure to dirty water (e.g., better sanitation in one patch), and/or rates of direct transmission
# Expect similar results to Robertson et al

### REDO THESE GRAPHS WITH NEW STRUCTURE - NON-DIMENSIONALISED FIRST
### ADD IN EPIDEMIC LENGTH AND SIZE FROM START
### MAKE MORE ROBUST CALCULATION OF EPIDEMIC SIZE

# Vary ΣaβW and ΣβD for each patch
# Direct and indirect transmission equal ΣpAjβDAj + ΣpBjβDBj = ΣpAjβWAaAj/κ + ΣpBjβWBaBj/κ
# pAj = pBj = 0.5, ΣβDAj + ΣβDBj = ΣβWAaAj + ΣβWBaBj; vary only indirect transmission rate βW

# Generate empty array for R0 outputs
patch_het_outputs <- array(data = NA, dim = c(100, 100))

# Create function to loop through transmission rates and calulate R0
patch.het <- function(betaWs, betaDs){
  for(i in 1:length(betaWs)){
  # Set betaW
  betaW <- betaWs[i]
    for(j in 1:length(betaDs)){
    # Set betaD
    betaD <- betaDs[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
    pAL <- 0.5
    pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
    pBL <- 0.5
    b_D_AH <- 0.00003 # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- 0.00003 # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- 0.00003 # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- 0.00003 # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- 0.1 # Indirect transmission rate for patch A
    b_W_B <- 0.1 # Indirect transmission rate for patch B
    alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A * pAH + alpha_AL * N_A * pAL + alpha_BH * N_B * pBH + alpha_BL * N_B * pBL # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    # Initial homogeneous model where ΣaβW/κ = ΣβD 
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = b_D_AH * N_A * pAH,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = b_D_AL * N_A * pAL,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B * pBH,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B * pBL,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05, # Inter-patch mixing is 5% strength
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A * pAH / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
              aAL     = alpha_AL * N_A * pAL / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B * pBH / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B * pBL / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
              xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Set parameters at DFE
    dfeparameters <- c(SA  = 1, # Susceptibles in patch A
                       IAH = 0, # High-infectiousness infectious individuals in patch A
                       IAL = 0, # Low-infectiousness infectious individuals in patch A
                       RA  = 0, # Recovered-and-immune individuals in patch A
                       SB  = 1, # Susceptibles in patch B
                       IBH = 0, # High-infectiousness infectious individuals in patch B
                       IBL = 0, # Low-infectiousness infectious individuals in patch B
                       RB  = 0, # Recovered-and-immune individuals in patch B
                       W   = 0, # Bacterial concentration in the environmental reservoir
                       paras)   # Add other parameters
    
    # Convert dfe parameters to list
    dfeparams <- as.list(dfeparameters)
    
    # Calculate R0
    patch_het_outputs[i,j] <- get.R0(dfeparams)
    }
  }
  return(patch_het_outputs)
}

# Set range of indirect transmission rates - βW values
betaWs <- seq(0, 0.4*10^6, length.out = 100)

# Set range of direct transmission rates - βD values
betaDs <- seq(0, 0.2, length.out = 100)

# Calculate R0s
patch_het_outputs <- patch.het(betaWs, betaDs)

# Plot outputs
# Wrangle R0s
patch_het_R0_df <- melt(patch_het_outputs)
names(patch_het_R0_df) <- c("betaW", "betaD", "R0")
# Map to actual x/y variable values
patch_het_R0_df$betaW <- betaWs[patch_het_R0_df$betaW]/max(betaWs) * 100
patch_het_R0_df$betaD <- betaDs[patch_het_R0_df$betaD]/max(betaDs) * 100

# Plot heat map for R0s in varying degrees of - note switched axes for ease of interpretation
ggplot(patch_het_R0_df, aes(y = betaD, x = betaW, fill = R0)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0.75, 2.0)) +
  labs(y = "% direct transmission in patch A",
       x = "% indirect transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Direct transmission dominates
# Set range of indirect transmission rates - βW values
betaWs <- seq(0, 0.3*10^6, length.out = 100)

# Set range of direct transmission rates - βD values
betaDs <- seq(0, 0.3, length.out = 100)

# Calculate R0s
patch_het_outputs <- patch.het(betaWs, betaDs)

# Plot outputs
# Wrangle R0s
patch_het_R0_df <- melt(patch_het_outputs)
names(patch_het_R0_df) <- c("betaW", "betaD", "R0")
# Map to actual x/y variable values
patch_het_R0_df$betaW <- betaWs[patch_het_R0_df$betaW]/max(betaWs) * 100
patch_het_R0_df$betaD <- betaDs[patch_het_R0_df$betaD]/max(betaDs) * 100

# Plot heat map for R0s in varying degrees of - note switched axes for ease of interpretation
ggplot(patch_het_R0_df, aes(y = betaD, x = betaW, fill = R0)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0.8, 2.0)) +
  labs(y = "% direct transmission in patch A",
       x = "% indirect transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Indirect transmission dominates
# Set range of indirect transmission rates - βW values
betaWs <- seq(0, 0.5*10^6, length.out = 100)
# Set range of direct transmission rates - βD values
betaDs <- seq(0, 0.1, length.out = 100)

# Calculate R0s
patch_het_outputs <- patch.het(betaWs, betaDs)

# Plot outputs
# Wrangle R0s
patch_het_R0_df <- melt(patch_het_outputs)
names(patch_het_R0_df) <- c("betaW", "betaD", "R0")
# Map to actual x/y variable values
patch_het_R0_df$betaW <- betaWs[patch_het_R0_df$betaW]/max(betaWs) * 100
patch_het_R0_df$betaD <- betaDs[patch_het_R0_df$betaD]/max(betaDs) * 100

# Plot heat map for R0s in varying degrees of - note switched axes for ease of interpretation
ggplot(patch_het_R0_df, aes(y = betaD, x = betaW, fill = R0)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0.8, 2.0)) +
  labs(y = "% direct transmission in patch A",
       x = "% indirect transmission in patch A",
       fill = "R0") +
  theme_minimal()


# Individual-level heterogeneity only ------------------------------------------
# Vary pAHβDAH vs pALβDAL and pBHβDBH vs pBLβDBL such that ΣβDA = ΣβDB (still patch-level homogeneity), 
# and ΣpβD = 0.1 - proportion of superspreaders is inversely related to their infectiousness, if pAH = 0.5, βDAH = 0.2, but if pAH = 0.1, βDAH = 1
# This is related to the 80:20 rule
# Low infectiousness individuals retain the same transmissibility
# Analogy - within each site there is a varying degree of superspreading -  differing behavioural practices, some individuals preparing food etc

# Vary within-patch direct transmission heterogeneity

# Generate empty array for R0 outputs
direct_het_outputs <- array(data = NA, dim = c(100, 100))

# Create function to loop through transmission rates and calulate R0
direct.het <- function(pAHs, pBHs){
  for(i in 1:length(pAHs)){
    # Set pAH
    pAH <- pAHs[i]
    for(j in 1:length(pBHs)){
      # Set pBH
      pBH <- pBHs[j]
      
      # Set parameters
      paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
                betaWA  = 0.4*10^6,        # Indirect transmission rate to patch A
                betaWB  = 0.4*10^6,        # Indirect transmission rate to patch B
                kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
                pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
                pAL     = 1 - pAH,
                pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
                pBL     = 1 - pBH,
                betaDAH = 0.2*sqrt((1-pAH)/pAH),       # Direct transmission rate for high-infectiousness individuals in patch A
                betaDAL = 0.2*sqrt(pAH/(1-pAH)),       # Direct transmission rate for low-infectiousness individuals in patch A
                betaDBH = 0.2*sqrt((1-pBH)/pBH),       # Direct transmission rate for high-infectiousness individuals in patch B
                betaDBL = 0.2*sqrt(pBH/(1-pBH)),       # Direct transmission rate for low-infectiousness individuals in patch B
                zeta    = 0.5,       # Inter-patch mixing is 5% strength of intra-patch mixing
                gamma   = 1/4,        # Duration of infection of 4 days
                aAH     = 0.5*sqrt((1-pAH)/pAH),        # Shedding rate from high-infectiousness individuals in patch A
                aAL     = 0.5*sqrt(pAH/(1-pAH)),        # Shedding rate from low-infectiousness individuals in patch A
                aBH     = 0.5*sqrt((1-pBH)/pBH),        # Shedding rate from high-infectiousness individuals in patch B
                aBL     = 0.5*sqrt(pBH/(1-pBH)),        # Shedding rate from low-infectiousness individuals in patch B
                xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
      
      # Set parameters at DFE
      dfeparameters <- c(SA  = 1, # Susceptibles in patch A
                         IAH = 0, # High-infectiousness infectious individuals in patch A
                         IAL = 0, # Low-infectiousness infectious individuals in patch A
                         RA  = 0, # Recovered-and-immune individuals in patch A
                         SB  = 1, # Susceptibles in patch B
                         IBH = 0, # High-infectiousness infectious individuals in patch B
                         IBL = 0, # Low-infectiousness infectious individuals in patch B
                         RB  = 0, # Recovered-and-immune individuals in patch B
                         W   = 0, # Bacterial concentration in the environmental reservoir
                         paras)   # Add other parameters
      
      # Convert dfe parameters to list
      dfeparams <- as.list(dfeparameters)
      
      # Calculate R0
      direct_het_outputs[i,j] <- get.R0(dfeparams)
    }
  }
  return(direct_het_outputs)
}

# Set proportion in infectious class per patch
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Calculate R0s
direct_het_R0_df <- melt(direct.het(pAHs, pBHs))
names(direct_het_R0_df) <- c("pAH", "pBH", "R0")
# Map to actual x/y variable values
direct_het_R0_df$pAH <- pAHs[direct_het_R0_df$pAH]
direct_het_R0_df$pBH <- pBHs[direct_het_R0_df$pBH]

# Plot heat map for R0s in varying degrees of - note switched axes for ease of interpretation
ggplot(direct_het_R0_df, aes(y = pBH, x = pAH, fill = R0)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "proportion of superspreaders in patch B",
       x = "proportion of superspreaders in patch A",
       fill = "R0") +
  theme_minimal()

# Inverse relationship between shedding and transmission

# Generate empty array for R0 outputs
direct_het_outputs <- array(data = NA, dim = c(100, 100))

# Create function to loop through transmission rates and calulate R0
direct.het <- function(pAHs, pBHs){
  for(i in 1:length(pAHs)){
    # Set pAH
    pAH <- pAHs[i]
    for(j in 1:length(pBHs)){
      # Set pBH
      pBH <- pBHs[j]
      
      # Set parameters
      paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
                betaWA  = 0.4*10^6,        # Indirect transmission rate to patch A
                betaWB  = 0.4*10^6,        # Indirect transmission rate to patch B
                kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
                pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
                pAL     = 1 - pAH,
                pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
                pBL     = 1 - pBH,
                betaDAH = 0.2*sqrt(pAH/(1-pAH)),       # Direct transmission rate for high-infectiousness individuals in patch A
                betaDAL = 0.2*sqrt((1-pAH)/pAH),       # Direct transmission rate for low-infectiousness individuals in patch A
                betaDBH = 0.2*sqrt(pBH/(1-pBH)),       # Direct transmission rate for high-infectiousness individuals in patch B
                betaDBL = 0.2*sqrt((1-pBH)/pBH),       # Direct transmission rate for low-infectiousness individuals in patch B
                zeta    = 0.5,       # Inter-patch mixing is 5% strength of intra-patch mixing
                gamma   = 1/4,        # Duration of infection of 4 days
                aAH     = 0.5*sqrt((1-pAH)/pAH),        # Shedding rate from high-infectiousness individuals in patch A
                aAL     = 0.5*sqrt(pAH/(1-pAH)),        # Shedding rate from low-infectiousness individuals in patch A
                aBH     = 0.5*sqrt((1-pBH)/pBH),        # Shedding rate from high-infectiousness individuals in patch B
                aBL     = 0.5*sqrt(pBH/(1-pBH)),        # Shedding rate from low-infectiousness individuals in patch B
                xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
      
      # Set parameters at DFE
      dfeparameters <- c(SA  = 1, # Susceptibles in patch A
                         IAH = 0, # High-infectiousness infectious individuals in patch A
                         IAL = 0, # Low-infectiousness infectious individuals in patch A
                         RA  = 0, # Recovered-and-immune individuals in patch A
                         SB  = 1, # Susceptibles in patch B
                         IBH = 0, # High-infectiousness infectious individuals in patch B
                         IBL = 0, # Low-infectiousness infectious individuals in patch B
                         RB  = 0, # Recovered-and-immune individuals in patch B
                         W   = 0, # Bacterial concentration in the environmental reservoir
                         paras)   # Add other parameters
      
      # Convert dfe parameters to list
      dfeparams <- as.list(dfeparameters)
      
      # Calculate R0
      direct_het_outputs[i,j] <- get.R0(dfeparams)
    }
  }
  return(direct_het_outputs)
}

# Set proportion in infectious class per patch
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Calculate R0s
direct_het_R0_df <- melt(direct.het(pAHs, pBHs))
names(direct_het_R0_df) <- c("pAH", "pBH", "R0")
# Map to actual x/y variable values
direct_het_R0_df$pAH <- pAHs[direct_het_R0_df$pAH]
direct_het_R0_df$pBH <- pBHs[direct_het_R0_df$pBH]

# Plot heat map for R0s in varying degrees of - note switched axes for ease of interpretation
ggplot(direct_het_R0_df, aes(y = pBH, x = pAH, fill = R0)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "proportion of superspreaders in patch B",
       x = "proportion of superspreaders in patch A",
       fill = "R0") +
  theme_minimal()



# Covary βDAH and aAH positively - expect that more shedding promotes both direct and indirect transmission
# Given a linear positive relationship between βDAH and aAH, measure relationship between heterogeneity


# Both patch and individual-level heterogeneity ---------------------------
# Vary degree of patch-level and individual-level heterogeneity 
# Greater patch-level heterogeneity = ΣaβW and ΣβD more different in A vs B, may be positively or negatively related

# Greater individual-level heterogeneity = 

# Vary patch heterogeneity and degree of individual heterogeneity in each patch
# Generate empty array for R0 outputs
both_het_outputs <- array(data = NA, dim = c(100, 100))

# Setting same degree of heterogeneity in each patch does not impact overall dynamics

# Create function to loop through transmission rates and calulate R0
both.het <- function(betaWs, betaDs, pAH, pBH){
  for(i in 1:length(betaWs)){
    # Set betaW
    betaW <- betaWs[i]
    for(j in 1:length(betaDs)){
      # Set betaD
      betaD <- betaDs[j]
      
      # Set parameters
      paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
                betaWA  = betaW,        # Indirect transmission rate to patch A
                betaWB  = max(betaWs) - betaW,        # Indirect transmission rate to patch B
                kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
                pAH     = pAH,        # Equal proportion of high- and low-infectiousness individuals in patch A
                pAL     = 1-pAH,
                pBH     = pBH,        # Equal proportion of high- and low-infectiousness individuals in patch B
                pBL     = 1-pBH,
                betaDAH = betaD*(1-pAH)/pAH,       # Direct transmission rate for high-infectiousness individuals in patch A
                betaDAL = betaD*pAH/(1-pAH),       # Direct transmission rate for low-infectiousness individuals in patch A
                betaDBH = (max(betaDs) - betaD)*(1-pBH)/pBH,       # Direct transmission rate for high-infectiousness individuals in patch B
                betaDBL = (max(betaDs) - betaD)*pBH/(1-pBH),       # Direct transmission rate for low-infectiousness individuals in patch B
                zeta    = 0.05,       # Inter-patch mixing is 5% strength of intra-patch mixing
                gamma   = 1/4,        # Duration of infection of 4 days
                aAH     = 0.5*(1-pAH)/pAH,        # Shedding rate from high-infectiousness individuals in patch A
                aAL     = 0.5*pAH/(1-pAH),        # Shedding rate from low-infectiousness individuals in patch A
                aBH     = 0.5*(1-pBH)/pBH,        # Shedding rate from high-infectiousness individuals in patch B
                aBL     = 0.5*pBH/(1-pBH),        # Shedding rate from low-infectiousness individuals in patch B
                xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
      
      # Set parameters at DFE
      dfeparameters <- c(SA  = 1, # Susceptibles in patch A
                         IAH = 0, # High-infectiousness infectious individuals in patch A
                         IAL = 0, # Low-infectiousness infectious individuals in patch A
                         RA  = 0, # Recovered-and-immune individuals in patch A
                         SB  = 1, # Susceptibles in patch B
                         IBH = 0, # High-infectiousness infectious individuals in patch B
                         IBL = 0, # Low-infectiousness infectious individuals in patch B
                         RB  = 0, # Recovered-and-immune individuals in patch B
                         W   = 0, # Bacterial concentration in the environmental reservoir
                         paras)   # Add other parameters
      
      # Convert dfe parameters to list
      dfeparams <- as.list(dfeparameters)
      
      # Calculate R0
      both_het_outputs[i,j] <- get.R0(dfeparams)
    }
  }
  return(both_het_outputs)
}

# Set range of indirect transmission rates - βW values
betaWs <- seq(0, 0.5*10^6, length.out = 100)

# Set range of direct transmission rates - βD values
betaDs <- seq(0, 0.1, length.out = 100)

# Calculate R0s
both_het_R0_df <- melt(both.het(betaWs = betaWs, betaDs = betaDs, pAH = 0.05, pBH = 0.5))
names(both_het_R0_df) <- c("betaW", "betaD", "R0")
# Map to actual x/y variable values
both_het_R0_df$betaW <- betaWs[both_het_R0_df$betaW]/max(betaWs) * 100
both_het_R0_df$betaD <- betaDs[both_het_R0_df$betaD]/max(betaDs) * 100

# Plot heat map for R0s in varying degrees of - note switched axes for ease of interpretation
ggplot(both_het_R0_df, aes(y = betaD, x = betaW, fill = R0)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "% direct transmission in patch A",
       x = "% indirect transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Analogy - heterogeneity between sites and within sites 


# Investigate other outputs from heterogeneity ----------------------------
# Investigate impact on size and speed of epidemic

# Homogeneous epidemic
paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
          betaWA  = 0.4*10^6,        # Indirect transmission rate to patch A - scale by kappa
          betaWB  = 0.4*10^6,        # Indirect transmission rate to patch B - scale by kappa
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
          pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = 0.5,
          pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = 0.5,
          betaDAH = 0.2,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = 0.2,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = 0.2,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = 0.2,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.05,         # Inter-patch mixing is 5% strength of intra-patch mixing
          gamma   = 1/4,        # Duration of infection of 4 days
          aAH     = 0.5,        # Shedding rate from high-infectiousness individuals in patch A
          aAL     = 0.5,        # Shedding rate from low-infectiousness individuals in patch A
          aBH     = 0.5,        # Shedding rate from high-infectiousness individuals in patch B
          aBL     = 0.5,        # Shedding rate from low-infectiousness individuals in patch B
          xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir

# Set initial conditions
start <- c(SA  = 0.98,  # Susceptibles in patch A
           IAH = 0.01, # High-infectiousness infectious individuals in patch A
           IAL = 0.01, # Low-infectiousness infectious individuals in patch A
           RA  = 0,     # Recovered-and-immune individuals in patch A
           SB  = 0.98,  # Susceptibles in patch B
           IBH = 0.01, # High-infectiousness infectious individuals in patch B
           IBL = 0.01, # Low-infectiousness infectious individuals in patch B
           RB  = 0,     # Recovered-and-immune individuals in patch B
           W   = 0)  # Bacterial concentration in the environmental reservoir

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top of data frame
head(round(out, 3))

# Plot outputs - same for both patches and infectiousness classes,
plot(x = out$time, y = out$SA, ylab = "Fraction", xlab = "Time", 
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$W, col = "purple")

# Create summary table, combining compartments
out_summary <- out %>%
               mutate(S = SA + SB,
                      I = IAH + IAL + IBH + IBL,
                      R = RA + RB)


# Set parameters as previous parameters and compartments at disease-free equilibrium to evaluate R0
dfeparameters <- c(SA  = 1, # Susceptibles in patch A
                   IAH = 0, # High-infectiousness infectious individuals in patch A
                   IAL = 0, # Low-infectiousness infectious individuals in patch A
                   RA  = 0, # Recovered-and-immune individuals in patch A
                   SB  = 1, # Susceptibles in patch B
                   IBH = 0, # High-infectiousness infectious individuals in patch B
                   IBL = 0, # Low-infectiousness infectious individuals in patch B
                   RB  = 0, # Recovered-and-immune individuals in patch B
                   W   = 0, # Bacterial concentration in the environmental reservoir
                   paras)   # Add other parameters

# Convert dfe parameters to list
dfeparameters <- as.list(dfeparameters)

# Get R0 for initial model
get.R0(dfeparameters)

# Calculate peak epidemic size
peak <- max(out_summary$I)

# Calculate time that peak occurred
peak_time <- out_summary$time[which(out_summary$I == max(out_summary$I))]

# Output peak sizes and times under varying degrees of heterogeneity
# Create function to extract R0, epidemic peak and time
get.outputs <- function(params){
  
  # Set initial conditions
  start <- c(SA  = 0.99,  # Susceptibles in patch A
             IAH = 0.005, # High-infectiousness infectious individuals in patch A
             IAL = 0.005, # Low-infectiousness infectious individuals in patch A
             RA  = 0,     # Recovered-and-immune individuals in patch A
             SB  = 0.99,  # Susceptibles in patch B
             IBH = 0.005, # High-infectiousness infectious individuals in patch B
             IBL = 0.005, # Low-infectiousness infectious individuals in patch B
             RB  = 0,     # Recovered-and-immune individuals in patch B
             W   = 0)  # Bacterial concentration in the environmental reservoir
  
  # Set times to run simulation for
  times = seq(0, 365, by = 1)
  
  # Simulate dynamics
  out <- ode(y = start, times = times, func = msmod, parms = params)
  
  # Convert output to data frame
  out = as.data.frame(out)
  
  # Create summary table, combining compartments
  out_summary <- out %>%
    mutate(S = SA + SB,
           I = IAH + IAL + IBH + IBL,
           R = RA + RB)
  
  # Set parameters as previous parameters and compartments at disease-free equilibrium to evaluate R0
  dfeparameters <- c(SA  = 1, # Susceptibles in patch A
                     IAH = 0, # High-infectiousness infectious individuals in patch A
                     IAL = 0, # Low-infectiousness infectious individuals in patch A
                     RA  = 0, # Recovered-and-immune individuals in patch A
                     SB  = 1, # Susceptibles in patch B
                     IBH = 0, # High-infectiousness infectious individuals in patch B
                     IBL = 0, # Low-infectiousness infectious individuals in patch B
                     RB  = 0, # Recovered-and-immune individuals in patch B
                     W   = 0, # Bacterial concentration in the environmental reservoir
                     params)   # Add other parameters
  
  # Convert dfe parameters to list
  dfeparameters <- as.list(dfeparameters)
  
  # Get R0
  R0 <- get.R0(dfeparameters)
  
  # Calculate epidemic size
  size <- max(out_summary$R)
  
  # Calculate epidemic length
  length <- out_summary$time[which(out_summary$R == max(out_summary$R))]
  
  # Combine outputs into list
  outputs <- c(R0, size, length)
  
  # Return list
  return(outputs)
  
}

# Patch-level heterogeneity
# Vary βW and βD

# Generate empty arrays for outputs
patch_het_R0   <- array(data = NA, dim = c(50, 50))
patch_het_size <- array(data = NA, dim = c(50, 50))
patch_het_length <- array(data = NA, dim = c(50, 50))

# Set range of indirect transmission rates - βW values
betaWs <- seq(0, 0.8*10^6, length.out = 50)
# Set range of direct transmission rates - βD values
betaDs <- seq(0, 0.4, length.out = 50)

# Run through each combination and get outputs
for(i in 1:length(betaWs)){
  # Set betaW
  betaW <- betaWs[i]
  for(j in 1:length(betaDs)){
    # Set betaD
    betaD <- betaDs[j]
    
    # Set parameters
    paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
              betaWA  = betaW,        # Indirect transmission rate to patch A
              betaWB  = max(betaWs) - betaW,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = betaD,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = betaD,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = max(betaDs) - betaD,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = max(betaDs) - betaD,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05,       # Inter-patch mixing is 5% strength of intra-patch mixing
              gamma   = 1/4,        # Duration of infection of 4 days
              aAH     = 0.5,        # Shedding rate from high-infectiousness individuals in patch A
              aAL     = 0.5,        # Shedding rate from low-infectiousness individuals in patch A
              aBH     = 0.5,        # Shedding rate from high-infectiousness individuals in patch B
              aBL     = 0.5,        # Shedding rate from low-infectiousness individuals in patch B
              xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
    
    # Get outputs
    outputs <- get.outputs(params = paras)
    patch_het_R0[i,j] <- outputs[1]
    patch_het_size[i,j] <- outputs[2]
    patch_het_length[i,j] <- outputs[3]
  }
}
s
# Wrangle to data frame of R0s
patch_het_R0_df <- melt(patch_het_R0)
names(patch_het_R0_df) <- c("betaW", "betaD", "R0")
# Map to actual x/y variable values
patch_het_R0_df$betaW <- betaWs[patch_het_R0_df$betaW]/max(betaWs) * 100
patch_het_R0_df$betaD <- betaDs[patch_het_R0_df$betaD]/max(betaDs) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission
ggplot(patch_het_R0_df, aes(x = betaW, y = betaD, fill = R0)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(1.5, 2, 2.5, 3, 3.5)),
                       limits = c(1.5, 3.5)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_het_size_df <- melt(patch_het_size)
names(patch_het_size_df) <- c("betaW", "betaD", "size")
# Map to actual x/y variable values
patch_het_size_df$betaW <- betaWs[patch_het_size_df$betaW]/max(betaWs) * 100
patch_het_size_df$betaD <- betaDs[patch_het_size_df$betaD]/max(betaDs) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_het_size_df, aes(x = betaW, y = betaD, fill = size/2)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(0.5, 0.625, 0.75, 0.875, 1)),
                       limits = c(0.5, 1)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
patch_het_length_df <- melt(patch_het_length)
names(patch_het_length_df) <- c("betaW", "betaD", "epidemic_length")
# Map to actual x/y variable values
patch_het_length_df$betaW <- betaWs[patch_het_length_df$betaW]/max(betaWs) * 100
patch_het_length_df$betaD <- betaDs[patch_het_length_df$betaD]/max(betaDs) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_het_length_df, aes(x = betaW, y = betaD, fill = epidemic_length)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(75, 95, 115, 135, 155)),
                       limits = c(75, 155)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()

# Individual heterogeneity only
# Vary pH, βH and aH; keep overall constant

# Generate empty arrays for outputs
indiv_het_R0   <- array(data = NA, dim = c(100, 100))
indiv_het_size <- array(data = NA, dim = c(100, 100))
indiv_het_time <- array(data = NA, dim = c(100, 100))

# Set proportion in infectious class per patch
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Run through each combination and get outputs
for(i in 1:length(pAHs)){
  # Set pAH
  pAH <- pAHs[i]
  for(j in 1:length(pBHs)){
    # Set pBH
    pBH <- pBHs[j]
    
    # Set parameters
    paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
              betaWA  = 0.4*10^6,        # Indirect transmission rate to patch A
              betaWB  = 0.4*10^6,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
              pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
              pAL     = 1 - pAH,
              pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
              pBL     = 1 - pBH,
              betaDAH = 0.2*(1-pAH)/pAH,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = 0.2*pAH/(1-pAH),       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = 0.2*(1-pBH)/pBH,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = 0.2*pBH/(1-pBH),       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05,       # Inter-patch mixing is 5% strength of intra-patch mixing
              gamma   = 1/4,        # Duration of infection of 4 days
              aAH     = 0.5*(1-pAH)/pAH,        # Shedding rate from high-infectiousness individuals in patch A
              aAL     = 0.5*pAH/(1-pAH),        # Shedding rate from low-infectiousness individuals in patch A
              aBH     = 0.5*(1-pBH)/pBH,        # Shedding rate from high-infectiousness individuals in patch B
              aBL     = 0.5*pBH/(1-pBH),     # Shedding rate from low-infectiousness individuals in patch B
              xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
    
    # Get outputs
    outputs <- get.outputs(params = paras)
    indiv_het_R0[i,j] <- outputs[1]
    indiv_het_size[i,j] <- outputs[2]
    indiv_het_time[i,j] <- outputs[3]
  }
}

# Wrangle to data frame of R0s
indiv_het_R0_df <- melt(indiv_het_R0)
names(indiv_het_R0_df) <- c("pAH", "pBH", "R0")
# Map to actual x/y variable values
indiv_het_R0_df$pAH <- pAHs[indiv_het_R0_df$pAH]
indiv_het_R0_df$pBH <- pBHs[indiv_het_R0_df$pBH]

# Plot heat map for R0s for different proportions of superspreaders
ggplot(indiv_het_R0_df, aes(x = pAH, y = pBH, fill = R0)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(1.5, 2, 2.5, 3, 3.5)),
                       limits = c(1.5, 3.5)) +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
indiv_het_size_df <- melt(indiv_het_size)
names(indiv_het_size_df) <- c("pAH", "pBH", "epidemic_size")
# Map to actual x/y variable values
indiv_het_size_df$pAH <- pAHs[indiv_het_size_df$pAH]
indiv_het_size_df$pBH <- pBHs[indiv_het_size_df$pBH]

# Plot heat map for sizes for different proportions of superspreaders
ggplot(indiv_het_size_df, aes(x = pAH, y = pBH, fill = epidemic_size/2)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(0.84, 0.88, 0.92, 0.96, 1)),
                       limits = c(0.84, 1)) +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
indiv_het_time_df <- melt(indiv_het_time)
names(indiv_het_time_df) <- c("pAH", "pBH", "epidemic_length")
# Map to actual x/y variable values
indiv_het_time_df$pAH <- pAHs[indiv_het_time_df$pAH]
indiv_het_time_df$pBH <- pBHs[indiv_het_time_df$pBH]

# Plot heat map for epidemic lengths for different proportions of superspreaders
ggplot(indiv_het_time_df, aes(x = pAH, y = pBH, fill =  epidemic_length)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(75, 95, 115, 135, 155)),
                       limits = c(75, 155)) +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic length") +
  theme_minimal()

# Individual and patch heterogeneity
# 5% superspreaders in both patches; effect of differing degrees of patch-level heterogeneity

# Generate empty arrays for outputs
patch_het_5_R0   <- array(data = NA, dim = c(50, 50))
patch_het_5_size <- array(data = NA, dim = c(50, 50))
patch_het_5_length <- array(data = NA, dim = c(50, 50))

# Set range of indirect transmission rates - βW values
betaWs <- seq(0, 0.8*10^6, length.out = 50)
# Set range of direct transmission rates - βD values
betaDs <- seq(0, 0.4, length.out = 50)

# Run through each combination and get outputs
for(i in 1:length(betaWs)){
  # Set betaW
  betaW <- betaWs[i]
  for(j in 1:length(betaDs)){
    # Set betaD
    betaD <- betaDs[j]
    
    # Set parameters
    paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
              betaWA  = betaW,        # Indirect transmission rate to patch A
              betaWB  = max(betaWs) - betaW,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
              pAH     = 0.05,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.95,
              pBH     = 0.05,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.95,
              betaDAH = betaD*(0.95/0.05),       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = betaD*(0.05/0.95),       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = (max(betaDs) - betaD)*(0.95/0.05),       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = (max(betaDs) - betaD)*(0.05/0.95),       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05,       # Inter-patch mixing is 5% strength of intra-patch mixing
              gamma   = 1/4,        # Duration of infection of 4 days
              aAH     = 0.5*(0.95/0.05),        # Shedding rate from high-infectiousness individuals in patch A
              aAL     = 0.5*(0.05/0.95),        # Shedding rate from low-infectiousness individuals in patch A
              aBH     = 0.5*(0.95/0.05),        # Shedding rate from high-infectiousness individuals in patch B
              aBL     = 0.5*(0.05/0.95),        # Shedding rate from low-infectiousness individuals in patch B
              xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
    
    # Get outputs
    outputs <- get.outputs(params = paras)
    patch_het_5_R0[i,j] <- outputs[1]
    patch_het_5_size[i,j] <- outputs[2]
    patch_het_5_length[i,j] <- outputs[3]
  }
}

# Wrangle to data frame of R0s
patch_het_5_R0_df <- melt(patch_het_5_R0)
names(patch_het_5_R0_df) <- c("betaW", "betaD", "R0")
# Map to actual x/y variable values
patch_het_5_R0_df$betaW <- betaWs[patch_het_5_R0_df$betaW]/max(betaWs) * 100
patch_het_5_R0_df$betaD <- betaDs[patch_het_5_R0_df$betaD]/max(betaDs) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission
ggplot(patch_het_5_R0_df, aes(x = betaW, y = betaD, fill = R0)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(1.5, 2, 2.5, 3, 3.5)),
                       limits = c(1.5, 3.5)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_het_5_size_df <- melt(patch_het_5_size)
names(patch_het_5_size_df) <- c("betaW", "betaD", "size")
# Map to actual x/y variable values
patch_het_5_size_df$betaW <- betaWs[patch_het_5_size_df$betaW]/max(betaWs) * 100
patch_het_5_size_df$betaD <- betaDs[patch_het_5_size_df$betaD]/max(betaDs) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_het_5_size_df, aes(x = betaW, y = betaD, fill = size/2)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(0.5, 0.625, 0.75, 0.875, 1)),
                       limits = c(0.5, 1)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
patch_het_5_length_df <- melt(patch_het_5_length)
names(patch_het_5_length_df) <- c("betaW", "betaD", "epidemic_length")
# Map to actual x/y variable values
patch_het_5_length_df$betaW <- betaWs[patch_het_5_length_df$betaW]/max(betaWs) * 100
patch_het_5_length_df$betaD <- betaDs[patch_het_5_length_df$betaD]/max(betaDs) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_het_5_length_df, aes(x = betaW, y = betaD, fill = epidemic_length)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(75, 95, 115, 135, 155)),
                       limits = c(75, 155)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()

# More direct transmission in patch A; effect of differing individual transmission heterogeneity
# Generate empty arrays for outputs
indiv_het_2_R0   <- array(data = NA, dim = c(100, 100))
indiv_het_2_size <- array(data = NA, dim = c(100, 100))
indiv_het_2_time <- array(data = NA, dim = c(100, 100))

# Set proportion in infectious class per patch
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Run through each combination and get outputs
for(i in 1:length(pAHs)){
  # Set pAH
  pAH <- pAHs[i]
  for(j in 1:length(pBHs)){
    # Set pBH
    pBH <- pBHs[j]
    
    # Set parameters
    paras = c(mu      = 1/(50*365), # Mean 50-year lifespan
              betaWA  = 0.7*10^6,        # Indirect transmission rate to patch A
              betaWB  = 0.1*10^6,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection   
              pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
              pAL     = 1 - pAH,
              pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
              pBL     = 1 - pBH,
              betaDAH = 0.38*(1-pAH)/pAH,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = 0.38*pAH/(1-pAH),       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = 0.02*(1-pBH)/pBH,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = 0.02*pBH/(1-pBH),       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05,       # Inter-patch mixing is 5% strength of intra-patch mixing
              gamma   = 1/4,        # Duration of infection of 4 days
              aAH     = 0.5*(1-pAH)/pAH,        # Shedding rate from high-infectiousness individuals in patch A
              aAL     = 0.5*pAH/(1-pAH),        # Shedding rate from low-infectiousness individuals in patch A
              aBH     = 0.5*(1-pBH)/pBH,        # Shedding rate from high-infectiousness individuals in patch B
              aBL     = 0.5*pBH/(1-pBH),     # Shedding rate from low-infectiousness individuals in patch B
              xi      = 30/365)     # Rate of loss of bacteria from the environmental reservoir
    
    # Get outputs
    outputs <- get.outputs(params = paras)
    indiv_het_2_R0[i,j] <- outputs[1]
    indiv_het_2_size[i,j] <- outputs[2]
    indiv_het_2_time[i,j] <- outputs[3]
  }
}

# Wrangle to data frame of R0s
indiv_het_2_R0_df <- melt(indiv_het_2_R0)
names(indiv_het_2_R0_df) <- c("pAH", "pBH", "R0")
# Map to actual x/y variable values
indiv_het_2_R0_df$pAH <- pAHs[indiv_het_2_R0_df$pAH]
indiv_het_2_R0_df$pBH <- pBHs[indiv_het_2_R0_df$pBH]

# Plot heat map for R0s for different proportions of superspreaders
ggplot(indiv_het_2_R0_df, aes(x = pAH, y = pBH, fill = R0)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(1.5, 2, 2.5, 3, 3.5)),
                       limits = c(1.5, 3.5)) +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
indiv_het_2_size_df <- melt(indiv_het_2_size)
names(indiv_het_2_size_df) <- c("pAH", "pBH", "epidemic_size")
# Map to actual x/y variable values
indiv_het_2_size_df$pAH <- pAHs[indiv_het_2_size_df$pAH]
indiv_het_2_size_df$pBH <- pBHs[indiv_het_2_size_df$pBH]

# Plot heat map for sizes for different proportions of superspreaders
ggplot(indiv_het_2_size_df, aes(x = pAH, y = pBH, fill = epidemic_size/2)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(0.6, 0.7, 0.8, 0.9, 1)),
                       limits = c(0.6, 1)) +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
indiv_het_2_time_df <- melt(indiv_het_2_time)
names(indiv_het_2_time_df) <- c("pAH", "pBH", "epidemic_length")
# Map to actual x/y variable values
indiv_het_2_time_df$pAH <- pAHs[indiv_het_2_time_df$pAH]
indiv_het_2_time_df$pBH <- pBHs[indiv_het_2_time_df$pBH]

# Plot heat map for epidemic lengths for different proportions of superspreaders
ggplot(indiv_het_2_time_df, aes(x = pAH, y = pBH, fill =  epidemic_length)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(75, 95, 115, 135, 155)),
                       limits = c(75, 155)) +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic length") +
  theme_minimal()

# Chain Patch Model -------------------------------------------------------
# Set up single patch model



