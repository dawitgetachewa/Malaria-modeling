# Load the library

library(deSolve)
library(tidyverse)

parameters <- c(beta = 0.25,
                gamma = 1/20,
                alpha= 0.12,
                mu_m= 1/15)

# Intial state
initial_state_values <- c(Sm = 40000,
                          Im = 3000,
                          S= 5000,
                          I= 1000)

# Define time of one year by days
times <- seq(0, 100, by = 0.01)

# model
SIS_MODEL <- function(time, state, parameters){
    with(as.list(c(state, parameters)), {
        M <- Sm + Im
        H<- S+I
        
        lambda_M <- beta * Im / H
        lambda_H<- alpha*I/H
        
        dSm <- mu_m*M-lambda_M * Sm - mu_m*Sm
        dIm<-  lambda_M * Sm-mu_m*Im
        
        dS<- -lambda_H*S+ gamma*I
        dI<- lambda_H*S- gamma*I
        return(list(c(dSm,dIm, dS, dI)))
    })
}

# convert it to data drame 
library(deSolve)
library(tidyverse)

# ---- Parameters ----
K <- 3   # mosquito-to-human ratio
parameters <- c(
    beta  = 0.5,      # human -> mosquito infection probability
    gamma = 1/25,     # human recovery rate
    alpha = 9e-05,    # mosquito -> human transmission scaling
    mu_m  = 1/10,     # mosquito birth/death rate
    Nh    = 1000,     # humans
    Nm    = 1000 * K  # mosquitoes
)

# ---- Initial state ----
initial_state_values <- c(
    Sm = parameters["Nm"] * 0.98,  # susceptible mosquitoes
    Im = parameters["Nm"] * 0.02,  # infectious mosquitoes
    Sh = 990,                      # susceptible humans
    Ih = 10                        # infectious humans
)

# ---- Time ----
times <- seq(0, 365, by = 1)

# ---- Model ----
SIS_MODEL <- function(time, state, pars) {
    with(as.list(c(state, pars)), {
        Mh <- Sh + Ih   # human population
        Mv <- Sm + Im   # mosquito population
        
        # Forces of infection
        lambda_H <- alpha * Im / Mv   # mosquito -> human
        lambda_M <- beta  * Ih / Mh   # human -> mosquito
        
        # Humans: SIS
        dSh <- -lambda_H * Sh + gamma * Ih
        dIh <-  lambda_H * Sh - gamma * Ih
        
        # Mosquitoes: SI with turnover
        dSm <- mu_m * Nm - lambda_M * Sm - mu_m * Sm
        dIm <- lambda_M * Sm - mu_m * Im
        
        list(c(dSm, dIm, dSh, dIh))
    })
}

# ---- Solve ----
out_put <- as.data.frame(ode(
    y     = initial_state_values,
    times = times,
    func  = SIS_MODEL,
    parms = parameters
))

head(out_put)

head(out_put)

# Pivot everything

out <- out_put %>%
    pivot_longer(
        cols = -time,              # pivot everything except time
        names_to = "Compartment",  # new column for variable name
        values_to = "Count"        # new column for values
    )

head(out)

# Plot SIS MODEL 
ggplot(out, aes(x = time, y = Count, color = Compartment)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("S" = "steelblue", "I" = "firebrick", Sm ="yellow", Im="green")) +
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

library(dplyr)
library(ggplot2)

out %>%
    ggplot(aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~Compartment, scales = "free_y") +
    labs(x = "Days", y = "Count", title = "Compartment Dynamics") +
    theme_minimal()

# SEASSONALITY
parameters <- c(beta = 0.25,
                gamma = 1/20,
                alpha= 0.12,
                mu_m= 1/15,
                amp= 0.5,
                phi= 1)

# Intial state
initial_state_values <- c(Sm = 1970,
                          Im = 30,
                          S= 990,
                          I= 10)

# Define time of one year by days
times <- seq(0, 1000, by = 1)

# model
SIS_MODEL <- function(time, state, parameters){
    with(as.list(c(state, parameters)), {
        M <- Sm + Im
        H<- S+I
        
        lambda_M <- beta * Im / H
        lambda_H<- alpha*I/H
        
        seas<-1+amp*cos(2*pi*(time/365-phi))
                        
        dSm <- mu_m*M-seas*lambda_M* Sm - mu_m*Sm
        dIm<-  lambda_M * Sm-mu_m*Im
        
        dS<- -seas*lambda_H*S+ gamma*I
        dI<- lambda_H*S- gamma*I
        
        return(list(c(dSm,dIm, dS, dI)))
    })
}

# convert it to data drame 
out_put <- as.data.frame(ode(y = initial_state_values,
                             times = times,
                             func = SIS_MODEL,
                             parms = parameters))

head(out_put)

# Pivot everything

out <- out_put %>%
    pivot_longer(
        cols = -time,              # pivot everything except time
        names_to = "Compartment",  # new column for variable name
        values_to = "Count"        # new column for values
    )

  # plot 
out %>%
    ggplot(aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~Compartment, scales = "free_y") +
    labs(x = "Days", y = "Count", title = "Compartment Dynamics") +
    theme_minimal()

