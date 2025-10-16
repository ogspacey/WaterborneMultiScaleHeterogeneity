### Multi-scale heterogeneity model
### Oliver G. Spacey
### 2025.10.15

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
times = seq(0, 365*50, by = 1)

# Set population sizes and original parameters inspired by Mukandivire et al.
N_A <- 10000 # Patch A population size
N_B <- 10000 # Patch B population size
pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
pAL <- 0.5
pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
pBL <- 0.5
b_D_AH <- 0.000015 # Direct transmission rate for high-infectiousness individuals in patch A
b_D_AL <- 0.000015 # Direct transmission rate for low-infectiousness individuals in patch A
b_D_BH <- 0.000015 # Direct transmission rate for high-infectiousness individuals in patch B
b_D_BL <- 0.000015 # Direct transmission rate for low-infectiousness individuals in patch B
b_W_A <- 0.05 # Indirect transmission rate for patch A
b_W_B <- 0.05 # Indirect transmission rate for patch B
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 10 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 10 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
  
# Calculate non-dimensionalised parameters
# Initial homogeneous model where ΣaβW/κ = ΣβD 
paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
          betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
          betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
          pAH     = pAH,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = pAL,
          pBH     = pBH,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = pBL,
          betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0, # Inter-patch mixing cut off for now
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate total indirect transmission - ΣaβW/κ 
tot_w_trnm <- as.numeric((paras['aAH']*paras['betaWA'] + paras['aAL']*paras['betaWA'] + paras['aBH']*paras['betaWB'] + paras['aBL']*paras['betaWB']) / paras['kappa'])

# Calculate total direct transmission - ΣβD 
tot_d_trnm <- as.numeric(paras['betaDAH'] + paras['betaDAL'] + paras['betaDBH'] + paras['betaDBL'])

# Calculate indirect transmission per patch - ΣpaβW/κ
a_w_trnm <- as.numeric((paras['pAH']*paras['aAH']*paras['betaWA'] + paras['pAL']*paras['aAL']*paras['betaWA']) / paras['kappa'])
b_w_trnm <- as.numeric((paras['pBH']*paras['aBH']*paras['betaWB'] + paras['pBL']*paras['aBL']*paras['betaWB']) / paras['kappa'])

# Calculate direct transmission per patch - ΣpβD
a_d_trnm <- as.numeric(paras['pAH']*paras['betaDAH'] + paras['pAL']*paras['betaDAL'])
b_d_trnm <- as.numeric(paras['pBH']*paras['betaDBH'] + paras['pBL']*paras['betaDBL'])


# Set initial conditions - epidemic conditions
start <- c(SA  = 0.99,  # Susceptibles in patch A
           IAH = 0.005, # High-infectiousness infectious individuals in patch A
           IAL = 0.005, # Low-infectiousness infectious individuals in patch A
           RA  = 0,     # Recovered-and-immune individuals in patch A
           SB  = 0.99,  # Susceptibles in patch B
           IBH = 0.005, # High-infectiousness infectious individuals in patch B
           IBL = 0.005, # Low-infectiousness infectious individuals in patch B
           RB  = 0,     # Recovered-and-immune individuals in patch B
           W   = 0)  # Bacterial concentration in the environmental reservoir

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 0.8 * max(out$time), y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")

# Basic Reproduction Number and Other Outputs ------------------------------------------------------------

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

# Create function to extract R0, epidemic peak and time
get.outputs <- function(start, params){
  
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
  
  # Calculate epidemic size - both patches recovered population after first outbreak
  size <- max(out_summary$R)
  
  # Calculate epidemic length
  length <- out_summary$time[which(out_summary$R == max(out_summary$R))]
  
  # Combine outputs into list
  outputs <- c(R0, size, length)
  
  # Return list
  return(outputs)
}

# Long-term Dynamics ----------------------------------------------

# Run initial model for 500 years
times = seq(0, 365*500, by = 1)

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Plot outputs - same for both patches and infectiousness classes, so just plot one
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 0.8*max(out$time), y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")

# Patch heterogeneity only ------------------------------------------------
# Vary patch heterogeneity in indirect and direct transmission
# Analogy - two sites differ in their exposure to dirty water (e.g., better sanitation in one patch), and/or rates of direct transmission
# Expect similar results to Robertson et al

# Vary ΣaβW and ΣβD for each patch
# Direct and indirect transmission equal ΣpAjβDAj + ΣpBjβDBj = ΣpAjβWAaAj/κ + ΣpBjβWBaBj/κ
# pAj = pBj = 0.5, ΣβDAj + ΣβDBj = ΣβWAaAj + ΣβWBaBj; vary only indirect transmission rate βW

# Generate empty arrays for outputs
patch_het_R0   <- array(data = NA, dim = c(50, 50))
patch_het_size <- array(data = NA, dim = c(50, 50))
patch_het_length <- array(data = NA, dim = c(50, 50))

# Set range of indirect transmission rates - bW values
b_Ws <- seq(0, 0.1, length.out = 50)

# Set range of direct transmission rates - bD values
b_Ds <- seq(0, 0.00003, length.out = 50)

# Loop through transmission rates and calculate outputs
for(i in 1:length(b_Ws)){
  # Set betaW
  b_W <- b_Ws[i]
    for(j in 1:length(b_Ds)){
    # Set betaD
    b_D <- b_Ds[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
    pAL <- 0.5
    pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
    pBL <- 0.5
    b_D_AH <- b_D # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- b_D # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- max(b_Ds) - b_D # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- max(b_Ds) - b_D # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- b_W # Indirect transmission rate for patch A
    b_W_B <- max(b_Ws) - b_W # Indirect transmission rate for patch B
    alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0, # Inter-patch mixing is 5% strength
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Set starting vector
      start <- c(SA  = 0.99,  # Susceptibles in patch A
                 IAH = 0.005, # High-infectiousness infectious individuals in patch A
                 IAL = 0.005, # Low-infectiousness infectious individuals in patch A
                 RA  = 0,     # Recovered-and-immune individuals in patch A
                 SB  = 0.99,  # Susceptibles in patch B
                 IBH = 0.005, # High-infectiousness infectious individuals in patch B
                 IBL = 0.005, # Low-infectiousness infectious individuals in patch B
                 RB  = 0,     # Recovered-and-immune individuals in patch B
                 W   = 0)  # Bacterial concentration in the environmental reservoir
    
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    patch_het_R0[i,j] <- outputs[1]
    patch_het_size[i,j] <- outputs[2]
    patch_het_length[i,j] <- outputs[3]
    }
}

# Plot outputs
# Wrangle to data frame of R0s
patch_het_R0_df <- melt(patch_het_R0)
names(patch_het_R0_df) <- c("b_W", "b_D", "R0")
# Map to actual x/y variable values
patch_het_R0_df$b_W <- b_Ws[patch_het_R0_df$b_W]/max(b_Ws) * 100
patch_het_R0_df$b_D <- b_Ds[patch_het_R0_df$b_D]/max(b_Ds) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission - set gradients for comparison later
ggplot(patch_het_R0_df, aes(x = b_W, y = b_D, fill = R0)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(1.2, 1.9, 2.6, 3.3, 4)),
  #                      limits = c(1.2, 4)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_het_size_df <- melt(patch_het_size)
names(patch_het_size_df) <- c("b_W", "b_D", "size")
# Map to actual x/y variable values
patch_het_size_df$b_W <- b_Ws[patch_het_size_df$b_W]/max(b_Ws) * 100
patch_het_size_df$b_D <- b_Ds[patch_het_size_df$b_D]/max(b_Ds) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_het_size_df, aes(x = b_W, y = b_D, fill = size/2)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(0.5, 0.625, 0.75, 0.875, 1)),
  #                      limits = c(0.5, 1)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
patch_het_length_df <- melt(patch_het_length)
names(patch_het_length_df) <- c("b_W", "b_D", "epidemic_length")
# Map to actual x/y variable values
patch_het_length_df$b_W <- b_Ws[patch_het_length_df$b_W]/max(b_Ws) * 100
patch_het_length_df$b_D <- b_Ds[patch_het_length_df$b_D]/max(b_Ds) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_het_length_df, aes(x = b_W, y = b_D, fill = epidemic_length)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(80, 130, 180, 230, 280)),
  #                      limits = c(80, 280)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()

# Direct transmission dominates

# Generate empty arrays for outputs
patch_het_R0   <- array(data = NA, dim = c(50, 50))
patch_het_size <- array(data = NA, dim = c(50, 50))
patch_het_length <- array(data = NA, dim = c(50, 50))

# Set range of indirect transmission rates - bW values
b_Ws <- seq(0, 0.05, length.out = 50)

# Set range of direct transmission rates - bD values
b_Ds <- seq(0, 0.00006, length.out = 50)

# Loop through transmission rates and calculate outputs
for(i in 1:length(b_Ws)){
  # Set betaW
  b_W <- b_Ws[i]
  for(j in 1:length(b_Ds)){
    # Set betaD
    b_D <- b_Ds[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
    pAL <- 0.5
    pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
    pBL <- 0.5
    b_D_AH <- b_D # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- b_D # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- max(b_Ds) - b_D # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- max(b_Ds) - b_D # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- b_W # Indirect transmission rate for patch A
    b_W_B <- max(b_Ws) - b_W # Indirect transmission rate for patch B
    alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05, # Inter-patch mixing is 5% strength
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Set starting vector
    start <- c(SA  = 0.98,  # Susceptibles in patch A
               IAH = 0.01, # High-infectiousness infectious individuals in patch A
               IAL = 0.01, # Low-infectiousness infectious individuals in patch A
               RA  = 0,     # Recovered-and-immune individuals in patch A
               SB  = 0.98,  # Susceptibles in patch B
               IBH = 0.01, # High-infectiousness infectious individuals in patch B
               IBL = 0.01, # Low-infectiousness infectious individuals in patch B
               RB  = 0,     # Recovered-and-immune individuals in patch B
               W   = 0)  # Bacterial concentration in the environmental reservoir
    
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    patch_het_R0[i,j] <- outputs[1]
    patch_het_size[i,j] <- outputs[2]
    patch_het_length[i,j] <- outputs[3]
  }
}

# Plot outputs
# Wrangle to data frame of R0s
patch_het_R0_df <- melt(patch_het_R0)
names(patch_het_R0_df) <- c("b_W", "b_D", "R0")
# Map to actual x/y variable values
patch_het_R0_df$b_W <- b_Ws[patch_het_R0_df$b_W]/max(b_Ws) * 100
patch_het_R0_df$b_D <- b_Ds[patch_het_R0_df$b_D]/max(b_Ds) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission
ggplot(patch_het_R0_df, aes(x = b_W, y = b_D, fill = R0)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(1.2, 1.9, 2.6, 3.3, 4)),
  #                      limits = c(1.2, 4)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_het_size_df <- melt(patch_het_size)
names(patch_het_size_df) <- c("b_W", "b_D", "size")
# Map to actual x/y variable values
patch_het_size_df$b_W <- b_Ws[patch_het_size_df$b_W]/max(b_Ws) * 100
patch_het_size_df$b_D <- b_Ds[patch_het_size_df$b_D]/max(b_Ds) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_het_size_df, aes(x = b_W, y = b_D, fill = size/2)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(0.5, 0.625, 0.75, 0.875, 1)),
  #                      limits = c(0.5, 1)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
patch_het_length_df <- melt(patch_het_length)
names(patch_het_length_df) <- c("b_W", "b_D", "epidemic_length")
# Map to actual x/y variable values
patch_het_length_df$b_W <- b_Ws[patch_het_length_df$b_W]/max(b_Ws) * 100
patch_het_length_df$b_D <- b_Ds[patch_het_length_df$b_D]/max(b_Ds) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_het_length_df, aes(x = b_W, y = b_D, fill = epidemic_length)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(80, 130, 180, 230, 280)),
  #                      limits = c(80, 280)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()

# Indirect transmission dominates
# Generate empty arrays for outputs
patch_het_R0   <- array(data = NA, dim = c(50, 50))
patch_het_size <- array(data = NA, dim = c(50, 50))
patch_het_length <- array(data = NA, dim = c(50, 50))

# Set range of indirect transmission rates - bW values
b_Ws <- seq(0, 0.2, length.out = 50)

# Set range of direct transmission rates - bD values
b_Ds <- seq(0, 0.000015, length.out = 50)

# Loop through transmission rates and calculate outputs
for(i in 1:length(b_Ws)){
  # Set betaW
  b_W <- b_Ws[i]
  for(j in 1:length(b_Ds)){
    # Set betaD
    b_D <- b_Ds[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
    pAL <- 0.5
    pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
    pBL <- 0.5
    b_D_AH <- b_D # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- b_D # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- max(b_Ds) - b_D # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- max(b_Ds) - b_D # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- b_W # Indirect transmission rate for patch A
    b_W_B <- max(b_Ws) - b_W # Indirect transmission rate for patch B
    alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0.05, # Inter-patch mixing is 5% strength
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Set starting vector
    start <- c(SA  = 0.98,  # Susceptibles in patch A
               IAH = 0.01, # High-infectiousness infectious individuals in patch A
               IAL = 0.01, # Low-infectiousness infectious individuals in patch A
               RA  = 0,     # Recovered-and-immune individuals in patch A
               SB  = 0.98,  # Susceptibles in patch B
               IBH = 0.01, # High-infectiousness infectious individuals in patch B
               IBL = 0.01, # Low-infectiousness infectious individuals in patch B
               RB  = 0,     # Recovered-and-immune individuals in patch B
               W   = 0)  # Bacterial concentration in the environmental reservoir
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    patch_het_R0[i,j] <- outputs[1]
    patch_het_size[i,j] <- outputs[2]
    patch_het_length[i,j] <- outputs[3]
  }
}

# Plot outputs
# Wrangle to data frame of R0s
patch_het_R0_df <- melt(patch_het_R0)
names(patch_het_R0_df) <- c("b_W", "b_D", "R0")
# Map to actual x/y variable values
patch_het_R0_df$b_W <- b_Ws[patch_het_R0_df$b_W]/max(b_Ws) * 100
patch_het_R0_df$b_D <- b_Ds[patch_het_R0_df$b_D]/max(b_Ds) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission
ggplot(patch_het_R0_df, aes(x = b_W, y = b_D, fill = R0)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(1.2, 1.9, 2.6, 3.3, 4)),
  #                      limits = c(1.2, 4)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_het_size_df <- melt(patch_het_size)
names(patch_het_size_df) <- c("b_W", "b_D", "size")
# Map to actual x/y variable values
patch_het_size_df$b_W <- b_Ws[patch_het_size_df$b_W]/max(b_Ws) * 100
patch_het_size_df$b_D <- b_Ds[patch_het_size_df$b_D]/max(b_Ds) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_het_size_df, aes(x = b_W, y = b_D, fill = size/2)) +
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
names(patch_het_length_df) <- c("b_W", "b_D", "epidemic_length")
# Map to actual x/y variable values
patch_het_length_df$b_W <- b_Ws[patch_het_length_df$b_W]/max(b_Ws) * 100
patch_het_length_df$b_D <- b_Ds[patch_het_length_df$b_D]/max(b_Ds) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_het_length_df, aes(x = b_W, y = b_D, fill = epidemic_length)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
                       values = rescale(c(80, 130, 180, 230, 280)),
                       limits = c(80, 280)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()

# Individual-level heterogeneity only ------------------------------------------
# Vary pAHβDAH vs pALβDAL and pBHβDBH vs pBLβDBL such that ΣβDA = ΣβDB (still patch-level homogeneity), 
# and ΣpβD = 0.1 - proportion of superspreaders is inversely related to their infectiousness, if pAH = 0.5, βDAH = 0.2, but if pAH = 0.1, βDAH = 1
# This is related to the 80:20 rule
# Analogy - within each site there is a varying degree of superspreading -  differing behavioural practices, some individuals preparing food etc
# Covary βDAH and aAH positively - expect that more shedding promotes both direct and indirect transmission
# Vary within-patch direct transmission heterogeneity
# Shedding and transmissibility both higher for superspreaders

# Vary degree of superspreader

# Set proportions in infectiousness classes
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Generate empty arrays for outputs
indiv_het_R0   <- array(data = NA, dim = c(100, 100))
indiv_het_size <- array(data = NA, dim = c(100, 100))
indiv_het_length <- array(data = NA, dim = c(100, 100))

# Loop through transmission rates and calculate outputs
for(i in 1:length(pAHs)){
    # Set pAH
    pAH <- pAHs[i]
    for(j in 1:length(pBHs)){
      # Set pBH
      pBH <- pBHs[j]
      
      # Set population sizes and original parameters inspired by Mukandivire et al.
      N_A <- 10000 # Patch A population size
      N_B <- 10000 # Patch B population size
      pAH <- pAH # Proportion of high- and low-infectiousness individuals in patch A
      pAL <- 1 - pAH
      pBH <- pBH # Proportion of high- and low-infectiousness individuals in patch B
      pBL <- 1 - pBH
      b_D_AH <- 0.000015 * 0.5 / pAH # Direct transmission rate for high-infectiousness individuals in patch A
      b_D_AL <- 0.000015 * 0.5 / pAL # Direct transmission rate for low-infectiousness individuals in patch A
      b_D_BH <- 0.000015 * 0.5 / pBH # Direct transmission rate for high-infectiousness individuals in patch B
      b_D_BL <- 0.000015 * 0.5 / pBL # Direct transmission rate for low-infectiousness individuals in patch B
      b_W_A <- 0.05 # Indirect transmission rate for patch A
      b_W_B <- 0.05 # Indirect transmission rate for patch B
      alpha_AH <- 10 * 0.5 / pAH # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
      alpha_AL <- 10 * 0.5 / pAL # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
      alpha_BH <- 10 * 0.5 / pBH # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
      alpha_BL <- 10 * 0.5 / pBL # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
      tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
      xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
      
      # Calculate non-dimensionalised parameters
      paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
                betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
                betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
                kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
                pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
                pAL     = pAL,
                pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
                pBL     = pBL,
                betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A - scale by pAL/pAH to increase transmissibility for this class
                betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
                betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
                betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
                zeta    = 0, # No inter-patch mixing 
                gamma   = 1/5, # Duration of infection of 5 days
                aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding  - scale by pAL/pAH to increase transmissibility for this class
                aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
                aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
                aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
      
      # Get outputs
      outputs <- get.outputs(start = start, params = paras)
      indiv_het_R0[i,j] <- outputs[1]
      indiv_het_size[i,j] <- outputs[2]
      indiv_het_length[i,j] <- outputs[3]
      
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
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
indiv_het_length_df <- melt(indiv_het_length)
names(indiv_het_length_df) <- c("pAH", "pBH", "epidemic_length")
# Map to actual x/y variable values
indiv_het_length_df$pAH <- pAHs[indiv_het_length_df$pAH]
indiv_het_length_df$pBH <- pBHs[indiv_het_length_df$pBH]

# Plot heat map for epidemic lengths for different proportions of superspreaders
ggplot(indiv_het_length_df, aes(x = pAH, y = pBH, fill =  epidemic_length)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic length") +
  theme_minimal()


# Both patch and individual-level heterogeneity ---------------------------

# Vary degree of patch-level heterogeneity under conditions of individual-level heterogeneity 
# Even individual heterogeneity - 20% superspreaders in both patches

# Generate empty arrays for outputs
patch_indiv_even_het_R0   <- array(data = NA, dim = c(100, 100))
patch_indiv_even_het_size <- array(data = NA, dim = c(100, 100))
patch_indiv_even_het_length <- array(data = NA, dim = c(100, 100))

# Set range of indirect transmission rates - bW values
b_Ws <- seq(0, 0.1, length.out = 100)

# Set range of direct transmission rates - bD values
b_Ds <- seq(0, 0.00003, length.out = 100)

# Loop through transmission rates and calculate outputs
for(i in 1:length(b_Ws)){
  # Set betaW
  b_W <- b_Ws[i]
  for(j in 1:length(b_Ds)){
    # Set betaD
    b_D <- b_Ds[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- 0.2 # Biased proportion of high- and low-infectiousness individuals in patch A
    pAL <- 0.8
    pBH <- 0.2 # Biased proportion of high- and low-infectiousness individuals in patch B
    pBL <- 0.8
    b_D_AH <- b_D * 0.5 / pAH # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- b_D * 0.5 / pAL # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- (max(b_Ds) - b_D) * 0.5 / pBH # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- (max(b_Ds) - b_D) * 0.5 / pBL # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- b_W # Indirect transmission rate for patch A
    b_W_B <- max(b_Ws) - b_W # Indirect transmission rate for patch B
    alpha_AH <- 10 * 0.5 / pAH # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 * 0.5 / pAL # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 * 0.5 / pBH # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 * 0.5 / pBL # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0, # Inter-patch mixing is 5% strength
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Set starting vector
    start <- c(SA  = 0.99,  # Susceptibles in patch A
               IAH = 0.005, # High-infectiousness infectious individuals in patch A
               IAL = 0.005, # Low-infectiousness infectious individuals in patch A
               RA  = 0,     # Recovered-and-immune individuals in patch A
               SB  = 0.99,  # Susceptibles in patch B
               IBH = 0.005, # High-infectiousness infectious individuals in patch B
               IBL = 0.005, # Low-infectiousness infectious individuals in patch B
               RB  = 0,     # Recovered-and-immune individuals in patch B
               W   = 0)  # Bacterial concentration in the environmental reservoir
    
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    patch_indiv_even_het_R0[i,j] <- outputs[1]
    patch_indiv_even_het_size[i,j] <- outputs[2]
    patch_indiv_even_het_length[i,j] <- outputs[3]
  }
}

# Plot outputs
# Wrangle to data frame of R0s
patch_indiv_even_het_R0_df <- melt(patch_indiv_even_het_R0)
names(patch_indiv_even_het_R0_df) <- c("b_W", "b_D", "R0")
# Map to actual x/y variable values
patch_indiv_even_het_R0_df$b_W <- b_Ws[patch_indiv_even_het_R0_df$b_W]/max(b_Ws) * 100
patch_indiv_even_het_R0_df$b_D <- b_Ds[patch_indiv_even_het_R0_df$b_D]/max(b_Ds) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission - set gradients for comparison later
ggplot(patch_indiv_even_het_R0_df, aes(x = b_W, y = b_D, fill = R0)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(1.2, 1.9, 2.6, 3.3, 4)),
  #                      limits = c(1.2, 4)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_indiv_even_het_size_df <- melt(patch_indiv_even_het_size)
names(patch_indiv_even_het_size_df) <- c("b_W", "b_D", "size")
# Map to actual x/y variable values
patch_indiv_even_het_size_df$b_W <- b_Ws[patch_indiv_even_het_size_df$b_W]/max(b_Ws) * 100
patch_indiv_even_het_size_df$b_D <- b_Ds[patch_indiv_even_het_size_df$b_D]/max(b_Ds) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_indiv_even_het_size_df, aes(x = b_W, y = b_D, fill = size/2)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(0.5, 0.625, 0.75, 0.875, 1)),
  #                      limits = c(0.5, 1)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
patch_indiv_even_het_length_df <- melt(patch_indiv_even_het_length)
names(patch_indiv_even_het_length_df) <- c("b_W", "b_D", "epidemic_length")
# Map to actual x/y variable values
patch_indiv_even_het_length_df$b_W <- b_Ws[patch_indiv_even_het_length_df$b_W]/max(b_Ws) * 100
patch_indiv_even_het_length_df$b_D <- b_Ds[patch_indiv_even_het_length_df$b_D]/max(b_Ds) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_indiv_even_het_length_df, aes(x = b_W, y = b_D, fill = epidemic_length)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(80, 130, 180, 230, 280)),
  #                      limits = c(80, 280)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()

# Uneven individual heterogeneity - 5% superspreaders in patch A only, patch B homogeneous

# Generate empty arrays for outputs
patch_indiv_uneven_het_R0   <- array(data = NA, dim = c(100, 100))
patch_indiv_uneven_het_size <- array(data = NA, dim = c(100, 100))
patch_indiv_uneven_het_length <- array(data = NA, dim = c(100, 100))

# Set range of indirect transmission rates - bW values
b_Ws <- seq(0, 0.1, length.out = 50)

# Set range of direct transmission rates - bD values
b_Ds <- seq(0, 0.00003, length.out = 50)

# Loop through transmission rates and calculate outputs
for(i in 1:length(b_Ws)){
  # Set betaW
  b_W <- b_Ws[i]
  for(j in 1:length(b_Ds)){
    # Set betaD
    b_D <- b_Ds[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- 0.2 # Biased proportion of high- and low-infectiousness individuals in patch A
    pAL <- 0.8
    pBH <- 0.5 # Biased proportion of high- and low-infectiousness individuals in patch B
    pBL <- 0.5
    b_D_AH <- b_D * 0.5 / pAH # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- b_D * 0.5 / pAL # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- (max(b_Ds) - b_D) * 0.5 / pBH # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- (max(b_Ds) - b_D) * 0.5 / pBL # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- b_W # Indirect transmission rate for patch A
    b_W_B <- max(b_Ws) - b_W # Indirect transmission rate for patch B
    alpha_AH <- 10 * 0.5 / pAH # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 * 0.5 / pAL # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 * 0.5 / pBH # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 * 0.5 / pBL # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
              pAL     = 0.5,
              pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
              pBL     = 0.5,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0, # Inter-patch mixing is 5% strength
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Set starting vector
    start <- c(SA  = 0.99,  # Susceptibles in patch A
               IAH = 0.005, # High-infectiousness infectious individuals in patch A
               IAL = 0.005, # Low-infectiousness infectious individuals in patch A
               RA  = 0,     # Recovered-and-immune individuals in patch A
               SB  = 0.99,  # Susceptibles in patch B
               IBH = 0.005, # High-infectiousness infectious individuals in patch B
               IBL = 0.005, # Low-infectiousness infectious individuals in patch B
               RB  = 0,     # Recovered-and-immune individuals in patch B
               W   = 0)  # Bacterial concentration in the environmental reservoir
    
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    patch_indiv_uneven_het_R0[i,j] <- outputs[1]
    patch_indiv_uneven_het_size[i,j] <- outputs[2]
    patch_indiv_uneven_het_length[i,j] <- outputs[3]
  }
}

# Plot outputs
# Wrangle to data frame of R0s
patch_indiv_uneven_het_R0_df <- melt(patch_indiv_uneven_het_R0)
names(patch_indiv_uneven_het_R0_df) <- c("b_W", "b_D", "R0")
# Map to actual x/y variable values
patch_indiv_uneven_het_R0_df$b_W <- b_Ws[patch_indiv_uneven_het_R0_df$b_W]/max(b_Ws) * 100
patch_indiv_uneven_het_R0_df$b_D <- b_Ds[patch_indiv_uneven_het_R0_df$b_D]/max(b_Ds) * 100

# Plot heat map for R0s for different contributions of indirect and direct transmission - set gradients for comparison later
ggplot(patch_indiv_uneven_het_R0_df, aes(x = b_W, y = b_D, fill = R0)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(1.2, 1.9, 2.6, 3.3, 4)),
  #                      limits = c(1.2, 4)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
patch_indiv_uneven_het_size_df <- melt(patch_indiv_uneven_het_size)
names(patch_indiv_uneven_het_size_df) <- c("b_W", "b_D", "size")
# Map to actual x/y variable values
patch_indiv_uneven_het_size_df$b_W <- b_Ws[patch_indiv_uneven_het_size_df$b_W]/max(b_Ws) * 100
patch_indiv_uneven_het_size_df$b_D <- b_Ds[patch_indiv_uneven_het_size_df$b_D]/max(b_Ds) * 100

# Plot heat map for sizes for different contributions of indirect and direct transmission
ggplot(patch_indiv_uneven_het_size_df, aes(x = b_W, y = b_D, fill = size/2)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(0.5, 0.625, 0.75, 0.875, 1)),
  #                      limits = c(0.5, 1)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
patch_indiv_uneven_het_length_df <- melt(patch_indiv_uneven_het_length)
names(patch_indiv_uneven_het_length_df) <- c("b_W", "b_D", "epidemic_length")
# Map to actual x/y variable values
patch_indiv_uneven_het_length_df$b_W <- b_Ws[patch_indiv_uneven_het_length_df$b_W]/max(b_Ws) * 100
patch_indiv_uneven_het_length_df$b_D <- b_Ds[patch_indiv_uneven_het_length_df$b_D]/max(b_Ds) * 100

# Plot heat map for epidemic lengths for different contributions of indirect and direct transmission
ggplot(patch_indiv_uneven_het_length_df, aes(x = b_W, y = b_D, fill = epidemic_length)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("blue", "green", "yellow", "orange", "red"),
  #                      values = rescale(c(80, 130, 180, 230, 280)),
  #                      limits = c(80, 280)) +
  labs(x = "% indirect transmission in patch A",
       y = "% direct transmission in patch A",
       fill = "epidemic length") +
  theme_minimal()


# Vary degree of individual-level heterogeneity under conditions of patch-level heterogeneity
# More direct and indirect transmission in patch A - direct and indirect positively covary
# Set proportions in infectiousness classes
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Generate empty arrays for outputs
indiv_patch_pos_cov_het_R0   <- array(data = NA, dim = c(100, 100))
indiv_patch_pos_cov_het_size <- array(data = NA, dim = c(100, 100))
indiv_patch_pos_cov_het_length <- array(data = NA, dim = c(100, 100))

# Loop through transmission rates and calculate outputs
for(i in 1:length(pAHs)){
  # Set pAH
  pAH <- pAHs[i]
  for(j in 1:length(pBHs)){
    # Set pBH
    pBH <- pBHs[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- pAH # Proportion of high- and low-infectiousness individuals in patch A
    pAL <- 1 - pAH
    pBH <- pBH # Proportion of high- and low-infectiousness individuals in patch B
    pBL <- 1 - pBH
    b_D_AH <- 0.0000225 * 0.5 / pAH # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- 0.0000225 * 0.5 / pAL # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- 0.0000075 * 0.5 / pBH # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- 0.0000075 * 0.5 / pBL # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- 0.075 # Indirect transmission rate for patch A
    b_W_B <- 0.025 # Indirect transmission rate for patch B
    alpha_AH <- 10 * 0.5 / pAH # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 * 0.5 / pAL # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 * 0.5 / pBH # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 * 0.5 / pBL # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
              pAL     = pAL,
              pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
              pBL     = pBL,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A - scale by pAL/pAH to increase transmissibility for this class
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0, # No inter-patch mixing 
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding  - scale by pAL/pAH to increase transmissibility for this class
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    indiv_patch_pos_cov_het_R0[i,j] <- outputs[1]
    indiv_patch_pos_cov_het_size[i,j] <- outputs[2]
    indiv_patch_pos_cov_het_length[i,j] <- outputs[3]
    
  }
}

# Wrangle to data frame of R0s
indiv_patch_pos_cov_het_R0_df <- melt(indiv_patch_pos_cov_het_R0)
names(indiv_patch_pos_cov_het_R0_df) <- c("pAH", "pBH", "R0")
# Map to actual x/y variable values
indiv_patch_pos_cov_het_R0_df$pAH <- pAHs[indiv_patch_pos_cov_het_R0_df$pAH]
indiv_patch_pos_cov_het_R0_df$pBH <- pBHs[indiv_patch_pos_cov_het_R0_df$pBH]

# Plot heat map for R0s for different proportions of superspreaders
ggplot(indiv_patch_pos_cov_het_R0_df, aes(x = pAH, y = pBH, fill = R0)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
indiv_patch_pos_cov_het_size_df <- melt(indiv_patch_pos_cov_het_size)
names(indiv_patch_pos_cov_het_size_df) <- c("pAH", "pBH", "epidemic_size")
# Map to actual x/y variable values
indiv_patch_pos_cov_het_size_df$pAH <- pAHs[indiv_patch_pos_cov_het_size_df$pAH]
indiv_patch_pos_cov_het_size_df$pBH <- pBHs[indiv_patch_pos_cov_het_size_df$pBH]

# Plot heat map for sizes for different proportions of superspreaders
ggplot(indiv_patch_pos_cov_het_size_df, aes(x = pAH, y = pBH, fill = epidemic_size/2)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
indiv_patch_pos_cov_het_length_df <- melt(indiv_patch_pos_cov_het_length)
names(indiv_patch_pos_cov_het_length_df) <- c("pAH", "pBH", "epidemic_length")
# Map to actual x/y variable values
indiv_patch_pos_cov_het_length_df$pAH <- pAHs[indiv_patch_pos_cov_het_length_df$pAH]
indiv_patch_pos_cov_het_length_df$pBH <- pBHs[indiv_patch_pos_cov_het_length_df$pBH]

# Plot heat map for epidemic lengths for different proportions of superspreaders
ggplot(indiv_patch_pos_cov_het_length_df, aes(x = pAH, y = pBH, fill =  epidemic_length)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic length") +
  theme_minimal()

# More direct transmission in patch A, more indirect transmission in patch B - direct and indirect negatively covary
# Set proportions in infectiousness classes
pAHs <- seq(0.01, 0.5, length.out = 100)
pBHs <- seq(0.01, 0.5, length.out = 100)

# Generate empty arrays for outputs
indiv_patch_neg_cov_het_R0   <- array(data = NA, dim = c(100, 100))
indiv_patch_neg_cov_het_size <- array(data = NA, dim = c(100, 100))
indiv_patch_neg_cov_het_length <- array(data = NA, dim = c(100, 100))

# Loop through transmission rates and calculate outputs
for(i in 1:length(pAHs)){
  # Set pAH
  pAH <- pAHs[i]
  for(j in 1:length(pBHs)){
    # Set pBH
    pBH <- pBHs[j]
    
    # Set population sizes and original parameters inspired by Mukandivire et al.
    N_A <- 10000 # Patch A population size
    N_B <- 10000 # Patch B population size
    pAH <- pAH # Proportion of high- and low-infectiousness individuals in patch A
    pAL <- 1 - pAH
    pBH <- pBH # Proportion of high- and low-infectiousness individuals in patch B
    pBL <- 1 - pBH
    b_D_AH <- 0.0000225 * 0.5 / pAH # Direct transmission rate for high-infectiousness individuals in patch A
    b_D_AL <- 0.0000225 * 0.5 / pAL # Direct transmission rate for low-infectiousness individuals in patch A
    b_D_BH <- 0.0000075 * 0.5 / pBH # Direct transmission rate for high-infectiousness individuals in patch B
    b_D_BL <- 0.0000075 * 0.5 / pBL # Direct transmission rate for low-infectiousness individuals in patch B
    b_W_A <- 0.025 # Indirect transmission rate for patch A
    b_W_B <- 0.075 # Indirect transmission rate for patch B
    alpha_AH <- 10 * 0.5 / pAH # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_AL <- 10 * 0.5 / pAL # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
    alpha_BH <- 10 * 0.5 / pBH # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    alpha_BL <- 10 * 0.5 / pBL # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
    tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
    xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span
    
    # Calculate non-dimensionalised parameters
    paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
              betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A
              betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B
              kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
              pAH     = pAH,        # Proportion of high- and low-infectiousness individuals in patch A
              pAL     = pAL,
              pBH     = pBH,        # Proportion of high- and low-infectiousness individuals in patch B
              pBL     = pBL,
              betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A - scale by pAL/pAH to increase transmissibility for this class
              betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
              betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
              betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
              zeta    = 0, # No inter-patch mixing 
              gamma   = 1/5, # Duration of infection of 5 days
              aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding  - scale by pAL/pAH to increase transmissibility for this class
              aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
              aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
              aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
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
    
    # Get outputs
    outputs <- get.outputs(start = start, params = paras)
    indiv_patch_neg_cov_het_R0[i,j] <- outputs[1]
    indiv_patch_neg_cov_het_size[i,j] <- outputs[2]
    indiv_patch_neg_cov_het_length[i,j] <- outputs[3]
    
  }
}

# Wrangle to data frame of R0s
indiv_patch_neg_cov_het_R0_df <- melt(indiv_patch_neg_cov_het_R0)
names(indiv_patch_neg_cov_het_R0_df) <- c("pAH", "pBH", "R0")
# Map to actual x/y variable values
indiv_patch_neg_cov_het_R0_df$pAH <- pAHs[indiv_patch_neg_cov_het_R0_df$pAH]
indiv_patch_neg_cov_het_R0_df$pBH <- pBHs[indiv_patch_neg_cov_het_R0_df$pBH]

# Plot heat map for R0s for different proportions of superspreaders
ggplot(indiv_patch_neg_cov_het_R0_df, aes(x = pAH, y = pBH, fill = R0)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "R0") +
  theme_minimal()

# Wrangle to data frame of sizes
indiv_patch_neg_cov_het_size_df <- melt(indiv_patch_neg_cov_het_size)
names(indiv_patch_neg_cov_het_size_df) <- c("pAH", "pBH", "epidemic_size")
# Map to actual x/y variable values
indiv_patch_neg_cov_het_size_df$pAH <- pAHs[indiv_patch_neg_cov_het_size_df$pAH]
indiv_patch_neg_cov_het_size_df$pBH <- pBHs[indiv_patch_neg_cov_het_size_df$pBH]

# Plot heat map for sizes for different proportions of superspreaders
ggplot(indiv_patch_neg_cov_het_size_df, aes(x = pAH, y = pBH, fill = epidemic_size/2)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic size") +
  theme_minimal()

# Wrangle to data frame of epidemic lengths
indiv_patch_neg_cov_het_length_df <- melt(indiv_patch_neg_cov_het_length)
names(indiv_patch_neg_cov_het_length_df) <- c("pAH", "pBH", "epidemic_length")
# Map to actual x/y variable values
indiv_patch_neg_cov_het_length_df$pAH <- pAHs[indiv_patch_neg_cov_het_length_df$pAH]
indiv_patch_neg_cov_het_length_df$pBH <- pBHs[indiv_patch_neg_cov_het_length_df$pBH]

# Plot heat map for epidemic lengths for different proportions of superspreaders
ggplot(indiv_patch_neg_cov_het_length_df, aes(x = pAH, y = pBH, fill =  epidemic_length)) +
  geom_tile() +
  labs(x = "proportion superspreaders in patch A",
       y = "proportion superspreaders in patch B",
       fill = "epidemic length") +
  theme_minimal()

# Chain Patch Model -------------------------------------------------------
# Set up single patch model

# Set times to run simulation for - run for long time to find equil
times = seq(0, 365*500, by = 1)

# Set population sizes and original parameters inspired by Mukandivire et al.
N_A <- 10000 # Patch A population size
N_B <- 10000 # Patch B population size
pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
pAL <- 0.5
pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
pBL <- 0.5
b_D_AH <- 0.000015 # Direct transmission rate for high-infectiousness individuals in patch A
b_D_AL <- 0.000015 # Direct transmission rate for low-infectiousness individuals in patch A
b_D_BH <- 0.000015 # Direct transmission rate for high-infectiousness individuals in patch B
b_D_BL <- 0.000015 # Direct transmission rate for low-infectiousness individuals in patch B
b_W_A <- 0.1 # Indirect transmission rate for patch A
b_W_B <- 0 # Indirect transmission rate for patch B - patch B entirely disconnected
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 0 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 0 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
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
          betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.05, # Inter-patch mixing is 5% strength
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate total indirect transmission - ΣaβW/κ 
tot_w_trnm <- as.numeric((paras['aAH'] * paras['betaWA'] + paras['aAL'] * paras['betaWA'] + paras['aBH'] * paras['betaWB'] + paras['aBL'] * paras['betaWB']) / paras['kappa'])

# Calculate total indirect transmission - ΣβD 
tot_d_trnm <- as.numeric(paras['betaDAH'] + paras['betaDAL'] + paras['betaDBH'] + paras['betaDBL'])

# Set initial conditions - epidemic conditions
start <- c(SA  = 0.99,  # Susceptibles in patch A
           IAH = 0.005, # High-infectiousness infectious individuals in patch A
           IAL = 0.005, # Low-infectiousness infectious individuals in patch A
           RA  = 0,     # Recovered-and-immune individuals in patch A
           SB  = 0.99,  # Susceptibles in patch B
           IBH = 0.005, # High-infectiousness infectious individuals in patch B
           IBL = 0.005, # Low-infectiousness infectious individuals in patch B
           RB  = 0,     # Recovered-and-immune individuals in patch B
           W   = 0)  # Bacterial concentration in the environmental reservoir

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top and bottom of data frame
head(round(out, 3))
tail(round(out, 3))

# Plot outputs
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", 
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 80000, y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")

# Get equilibrium value
equil_SA <- tail(out, n = 1)$SA
equil_IAH <- tail(out, n = 1)$IAH
equil_IAL <- tail(out, n = 1)$IAL
equil_RA <- tail(out, n = 1)$RA
# Test whether equilibrium estimated properly
paras['gamma']*(equil_IAH + equil_IAL)/paras['mu']




# Run with patch A already at equilibrium
# Set initial conditions - epidemic conditions
start <- c(SA  = equil_SA,  # Susceptibles in patch A
           IAH = equil_IAH, # High-infectiousness infectious individuals in patch A
           IAL = equil_IAL, # Low-infectiousness infectious individuals in patch A
           RA  = equil_RA,     # Recovered-and-immune individuals in patch A
           SB  = 0.99,  # Susceptibles in patch B
           IBH = 0.005, # High-infectiousness infectious individuals in patch B
           IBL = 0.005, # Low-infectiousness infectious individuals in patch B
           RB  = 0,     # Recovered-and-immune individuals in patch B
           W   = 0)  # Bacterial concentration in the environmental reservoir

# Reset parameters - increase indirect transmission, decrease direct transmission
# Set population sizes and original parameters inspired by Mukandivire et al.
N_A <- 10000 # Patch A population size
N_B <- 10000 # Patch B population size
pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
pAL <- 0.5
pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
pBL <- 0.5
b_D_AH <- 0.0000075 # Direct transmission rate for high-infectiousness individuals in patch A
b_D_AL <- 0.0000075 # Direct transmission rate for low-infectiousness individuals in patch A
b_D_BH <- 0.0000075 # Direct transmission rate for high-infectiousness individuals in patch B
b_D_BL <- 0.0000075 # Direct transmission rate for low-infectiousness individuals in patch B
b_W_A <- 0.2 # Indirect transmission rate for patch A
b_W_B <- 0 # Indirect transmission rate for patch B - patch B entirely disconnected
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 0 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 0 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate non-dimensionalised parameters
paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
          betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
          betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
          pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = 0.5,
          pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = 0.5,
          betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.05, # Inter-patch mixing is 5% strength
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top and bottom of data frame
head(round(out, 3))
tail(out, 3)

# Plot outputs
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", main = "Increase indirect transmission, decrease direct transmission",
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 0.8*max(out$time), y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")


# Reset parameters - increase direct transmission, indirect transmission
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
b_W_B <- 0 # Indirect transmission rate for patch B - patch B entirely disconnected
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 0 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 0 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate non-dimensionalised parameters
paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
          betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
          betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
          pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = 0.5,
          pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = 0.5,
          betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.05, # Inter-patch mixing is 5% strength
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span


# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top and bottom of data frame
head(round(out, 3))
tail(round(out, 3))

# Plot outputs
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", main = "Increase direct transmission, decrease indirect transmission",
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 0.8*max(out$time), y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")

# Reset parameters - increase zeta
# Set population sizes and original parameters inspired by Mukandivire et al.
N_A <- 10000 # Patch A population size
N_B <- 10000 # Patch B population size
pAH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch A
pAL <- 0.5
pBH <- 0.5 # Equal proportion of high- and low-infectiousness individuals in patch B
pBL <- 0.5
b_D_AH <- 0.000015 # Direct transmission rate for high-infectiousness individuals in patch A
b_D_AL <- 0.000015 # Direct transmission rate for low-infectiousness individuals in patch A
b_D_BH <- 0.000015 # Direct transmission rate for high-infectiousness individuals in patch B
b_D_BL <- 0.000015 # Direct transmission rate for low-infectiousness individuals in patch B
b_W_A <- 0.2 # Indirect transmission rate for patch A
b_W_B <- 0 # Indirect transmission rate for patch B - patch B entirely disconnected
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 0 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 0 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate non-dimensionalised parameters
paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
          betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
          betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
          pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = 0.5,
          pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = 0.5,
          betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.2, # Inter-patch mixing is 5% strength
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span


# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top and bottom of data frame
head(round(out, 3))
tail(round(out, 3))

# Plot outputs
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", main = "Increase zeta",
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 0.8*max(out$time), y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")

# Reset parameters - more direct transmission heterogeneity
# Set population sizes and original parameters inspired by Mukandivire et al.
N_A <- 10000 # Patch A population size
N_B <- 10000 # Patch B population size
pAH <- 0.8 # Equal proportion of high- and low-infectiousness individuals in patch A
pAL <- 0.2
pBH <- 0.8 # Equal proportion of high- and low-infectiousness individuals in patch B
pBL <- 0.2
b_D_AH <- 0.000015 * 0.5 / pAH # Direct transmission rate for high-infectiousness individuals in patch A
b_D_AL <- 0.000015 * 0.5 / pAL # Direct transmission rate for low-infectiousness individuals in patch A
b_D_BH <- 0.000015 * 0.5 / pBH # Direct transmission rate for high-infectiousness individuals in patch B
b_D_BL <- 0.000015 * 0.5 / pBL # Direct transmission rate for low-infectiousness individuals in patch B
b_W_A <- 0.2 # Indirect transmission rate for patch A
b_W_B <- 0 # Indirect transmission rate for patch B - patch B entirely disconnected
alpha_AH <- 10 # Shedding rate for high-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_AL <- 10 # Shedding rate for low-infectiousness individuals in patch A; cells mL^-1 d^-1 person^-1
alpha_BH <- 0 # Shedding rate for high-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
alpha_BL <- 0 # Shedding rate for low-infectiousness individuals in patch B; cells mL^-1 d^-1 person^-1 
tot_shed <- alpha_AH * N_A + alpha_AL * N_A + alpha_BH * N_B + alpha_BL * N_B # Total shedding from all individuals, for scaling
xi <- 1/30  # Rate of loss of bacteria from the environmental reservoir - 30 day life span

# Calculate non-dimensionalised parameters
paras = c(mu      = 1/(43.5*365), # Mean 43.5-year lifespan
          betaWA  = b_W_A * tot_shed / xi,        # Indirect transmission rate to patch A - scale by kappa for saturating effect
          betaWB  = b_W_B * tot_shed / xi,        # Indirect transmission rate to patch B - scale by kappa for saturating effect
          kappa   = 10^6,       # Saturating coefficient - infectious dose sufficient for 50% infection cells/mL 
          pAH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch A
          pAL     = 0.5,
          pBH     = 0.5,        # Equal proportion of high- and low-infectiousness individuals in patch B
          pBL     = 0.5,
          betaDAH = b_D_AH * N_A,       # Direct transmission rate for high-infectiousness individuals in patch A
          betaDAL = b_D_AL * N_A,       # Direct transmission rate for low-infectiousness individuals in patch A
          betaDBH = b_D_BH * N_B,       # Direct transmission rate for high-infectiousness individuals in patch B
          betaDBL = b_D_BL * N_B,       # Direct transmission rate for low-infectiousness individuals in patch B
          zeta    = 0.2, # Inter-patch mixing is 5% strength
          gamma   = 1/5, # Duration of infection of 5 days
          aAH     = alpha_AH * N_A / tot_shed,  # Shedding rate from high-infectiousness individuals in patch A; scaled as proportion of total shedding
          aAL     = alpha_AL * N_A / tot_shed,  # Shedding rate from low-infectiousness individuals in patch A; scaled as proportion of total shedding
          aBH     = alpha_BH * N_B / tot_shed,  # Shedding rate from high-infectiousness individuals in patch B; scaled as proportion of total shedding
          aBL     = alpha_BL * N_B / tot_shed,  # Shedding rate from low-infectiousness individuals in patch B; scaled as proportion of total shedding
          xi      = xi)     # Rate of loss of bacteria from the environmental reservoir - 30 day life span


# Simulate dynamics
out <- ode(y = start, times = times, func = msmod, parms = paras)

# Convert output to data frame
out = as.data.frame(out)

# Check top and bottom of data frame
head(round(out, 3))
tail(round(out, 3))

# Plot outputs
plot(x = out$time, y = out$SA, 
     ylab = "Fraction of population", xlab = "Time (days)", main = "Direct transmission heterogeneity",
     type = "l", ylim = c(0,1))
lines(x = out$time, y = out$IAH + out$IAL, col = "red") 
lines(x = out$time, y = out$RA, col = "blue")
lines(x = out$time, y = out$SB, col = "grey")
lines(x = out$time, y = out$IBH + out$IBL, col = "pink") 
lines(x = out$time, y = out$RB, col = "lightblue")
legend(x = 0.8*max(out$time), y = 0.9,                                
       legend = c("SA", "IA", "RA", "SB", "IB", "RB"),
       col = c("black", "red", "blue", "grey", "pink", "lightblue"),
       lty = 1,                                 # line type
       cex = 1.5,                               # text size
       bty = "n")                               # no box around legend
# Plot reservoir
plot(x = out$time, y = out$W * tot_shed/xi, col = "purple", type = "l", 
     xlab = "Time (days)", ylab = "Bacterial concentration in water")


