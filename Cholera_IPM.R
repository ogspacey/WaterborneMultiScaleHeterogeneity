### Cholera IPM
### Oliver G. Spacey
### 2025.10.15

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load required packages
library(ggplot2)   # data visualisation

# Vital Rate Models -----------------------------------------------------------
# Generate sequence for state variable (log10(bacterial concentration))
xx <- seq(0, 10, by = 0.1)

# Infection probability - logistic relationship between log(infectious dose) and infection probability
# Centre on 10^4 - Kothary and Babu, 2001
inf_prob <- function(x) {
  beta0 <- -4.5
  beta1 <- 1
  
  logit <- beta0 + beta1 * x
  prob <- 1 / (1 + exp(-logit))
  return(prob)
}

# Plot output
plot(x = xx, y = inf_prob(xx),
     xlab = "log10(Bacteria ingested)", ylab = "Probability of infection",
     type = "l")

# In-person bacterial growth

# Infected survival rate 

# Shedding rate

# Reservoir bacterial growth


# Construct IPM -----------------------------------------------------------


