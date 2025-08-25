# Load the library

library(deSolve)
library(tidyverse)
library(plotly)

parameters <- c(beta = 0.4,
                gamma = 0.1,
                alpha= 0.5)

# Intial state
initial_state_values <- c(S = 499,
                          I = 1
                          
                         )

# Define time of one year by days
times <- seq(0, 365, by = 1)

# model
SIS_MODEL <- function(time, state, parameters){
    with(as.list(c(state, parameters)), {
        N <- S + I
        
        
        lambda <- beta * I / N
        
        
        dS <- -lambda*S - gamma*S
        dI<-  lambda*S- gamma*I
        
        
        return(list(c(dS,dI)))
    })
}
# convert it to data drame 
out_put <- as.data.frame(ode(y = initial_state_values,
                             times = times,
                             func = SIS_MODEL,
                             parms = parameters))

head(out_put)

# Pivot everything
out<- out_put%>%
    group_by(time)%>%
    mutate(total_pop= sum(S,I))%>%
    pivot_longer(S:I)

head(out)

out <- out_put %>%
    pivot_longer(
        cols = -time,              # pivot everything except time
        names_to = "Compartment",  # new column for variable name
        values_to = "Count"        # new column for values
    )

head(out)
out[366,]

p1<-out %>%
    ggplot(aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~Compartment, scales = "free_y") +
    labs(x = "Days", y = "Count", 
         title = "Compartment Dynamics") +
    theme_minimal()

plotly::ggplotly(p1)
