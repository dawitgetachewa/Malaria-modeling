# Load the library

library(deSolve)
library(tidyverse)

# Set parameters for SIS model
parameters <- c(beta = 0.3,
                gamma = 0.1,
                delta= 0.3,
                N=500)

# Intial state
initial_state_values <- c(S = 490,
                          E= 9,
                          I = 1)
# Define time of one year by days
times <- seq(0, 365, by = 1)

# model
SEIS_MODEL <- function(time, state, parameters){
    with(as.list(c(state, parameters)), {
        N <- S + E + I
        lambda <- beta * I / N
        
        dS <- -lambda * S + gamma * I
        dE <-  lambda * S - delta * E
        dI <-  delta * E - gamma * I
        
        list(c(dS, dE, dI))
    })
}

# convert it to data drame 
out_put <- as.data.frame(ode(y = initial_state_values,
                             times = times,
                             func = SIS_MODEL,
                             parms = parameters))

head(out_put)

# Pivot everything

out_put_long <- out_put %>%
    pivot_longer(
        cols = -time,              # pivot everything except time
        names_to = "Compartment",  # new column for variable name
        values_to = "Count"        # new column for values
    )

head(out_put_long)

# Plot SIS MODEL 
ggplot(out_put_long, aes(x = time, y = Count, color = Compartment)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("S" = "steelblue", "I" = "firebrick", E ="yellow")) +
    theme_minimal(base_size = 14) +
    labs(title = "SIS Model Dynamics Over Time",
         subtitle = "Population split into Susceptible (S) and Infected (I)",
         x = "Time (days)",
         y = "Number of individuals",
         color = "Compartment") +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12),
        legend.position = "top"
    )

