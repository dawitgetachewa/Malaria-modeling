library(deSolve)

# ---------------------------
# Parameters
# ---------------------------
params <- list(
    Lam_h = 0.033,        # Humans × Days^-1
    psi_h = 1.1e-4,       # Days^-1
    psi_v = 0.13,         # Days^-1
    sigma_v = 0.50,       # Days^-1
    sigma_h = 19.0,       # Days^-1
    beta_hv = 0.022,      # prob mosquito->human
    beta_vh = 0.48,       # prob human->mosquito (Ih)
    beta_tilde_vh = 0.048,# prob human->mosquito (Rh carrier)
    nu_h = 0.10,          # Days^-1 (Eh -> Ih)
    nu_v = 0.091,         # Days^-1 (Ev -> Iv)
    gamma_h = 0.0035,     # Days^-1 (Ih -> Rh)
    delta_h = 9.0e-5,     # Days^-1 (disease-induced human death)
    rho_h = 5.5e-4,       # Days^-1 (Rh -> Sh)
    mu1h = 1.6e-5,        # Days^-1 (human density-independent death)
    mu2h = 3.0e-7,        # Humans^-1 Days^-1
    mu1v = 0.033,         # Days^-1 (mosquito density-independent death)
    mu2v = 2.0e-5         # Mosquitoes^-1 Days^-1
)

# ---------------------------
# Initial conditions
# ---------------------------
state <- c(
    Sh = 600,
    Eh = 20,
    Ih = 3,
    Rh = 0,
    Sv = 2400,
    Ev = 30,
    Iv = 5
)

# ---------------------------
# ODE system
# ---------------------------
model <- function(t, y, p) {
    with(as.list(c(y, p)), {
        
        # Populations
        Nh <- Sh + Eh + Ih + Rh
        Nv <- Sv + Ev + Iv
        
        # Density-dependent death rates
        mu_h <- mu1h + mu2h * Nh
        mu_v <- mu1v + mu2v * Nv
        
        # Force of infection
        lambda_h <- (sigma_v * sigma_h * Nv) / (sigma_v * Nv + sigma_h * Nh) * 
            beta_hv * (Iv / Nv)
        
        lambda_v <- (sigma_v * sigma_h * Nh) / (sigma_v * Nv + sigma_h * Nh) * 
            ((beta_vh * Ih + beta_tilde_vh * Rh) / Nh)
        
        # Human compartments
        dSh <- Lam_h + psi_h * Nh + rho_h * Rh - lambda_h * Sh - mu_h * Sh
        dEh <- lambda_h * Sh - nu_h * Eh - mu_h * Eh
        dIh <- nu_h * Eh - gamma_h * Ih - mu_h * Ih - delta_h * Ih
        dRh <- gamma_h * Ih - rho_h * Rh - mu_h * Rh
        
        # Vector compartments
        dSv <- psi_v * Nv - lambda_v * Sv - mu_v * Sv
        dEv <- lambda_v * Sv - nu_v * Ev - mu_v * Ev
        dIv <- nu_v * Ev - mu_v * Iv
        
        list(c(dSh, dEh, dIh, dRh, dSv, dEv, dIv))
    })
}

# ---------------------------
# Simulation
# ---------------------------
times <- seq(0, 10000, by = 2)  # 1 year, daily steps
out <- ode(y = state, times = times, func = model, parms = params)

out <- as.data.frame(out)

# ---------------------------
# Plot results
# ---------------------------
par(mfrow = c(1, 2))

# Humans
plot(out$time, out$Ih, type="l", col="blue", ylim=c(0, max(out$Sh, out$Eh, out$Ih, out$Rh)),
     xlab="Days", ylab="Humans", main="Human compartments")
lines(out$time, out$Eh, col="orange")
lines(out$time, out$Ih, col="red")
lines(out$time, out$Rh, col="green")
legend("right", legend=c("Sh","Eh","Ih","Rh"),
       col=c("blue","orange","red","green"), lty=1, cex=0.8)

# Mosquitoes
plot(out$time, out$Iv, type="l", col="blue", ylim=c(0, max(out$Sv, out$Ev, out$Iv)),
     xlab="Days", ylab="Mosquitoes", main="Mosquito compartments")
lines(out$time, out$Ev, col="orange")
lines(out$time, out$Iv, col="red")
legend("right", legend=c("Sv","Ev","Iv"),
       col=c("blue","orange","red"), lty=1, cex=0.8)

