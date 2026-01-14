getwd()
lv_data <- read.csv("./tmapper_blog/hudson-bay-lynx-hare.csv", skip = 2, header = TRUE)
lv_data

library(deSolve)

# Plot 1: Hare and Lynx time series
svg("./tmapper_blog/images/hare_lynx.svg", width = 10, height = 5)
par(mar = c(4.5, 4.5, 2, 2), mgp = c(2.5, 0.8, 0))

plot(lv_data$Year, lv_data$Hare,
     type = "o",
     col = "#2E86C1",
     pch = 16,
     cex = 1.2,
     xlab = "Year",
     ylab = "Number of pelts x 1000",
     main = "Hudson Bay Hare and Lynx Populations (1821-1914)",
     ylim = range(c(lv_data$Hare, lv_data$Lynx)),
     lwd = 2,
     cex.main = 1.3,
     cex.axis = 1.1,
     cex.lab = 1.1)

lines(lv_data$Year, lv_data$Lynx,
      type = "o",
      col = "#C0392B",
      pch = 16,
      cex = 1.2,
      lwd = 2)

legend("topleft",
       legend = c("Hare", "Lynx"),
       col = c("#2E86C1", "#C0392B"),
       lwd = 2,
       pch = 16,
       cex = 1.0,
       bty = "n")

grid(nx = NULL, ny = NULL, col = "gray85", lty = "solid", lwd = 0.5)

dev.off()


# LV system
lv_ode <- function(t, state, parameters) {
    # Ensure all inputs are numeric to avoid type errors
    x <- as.numeric(state["x"])
    y <- as.numeric(state["y"])
    alpha <- as.numeric(parameters["alpha"])
    beta  <- as.numeric(parameters["beta"])
    delta <- as.numeric(parameters["delta"])
    gamma <- as.numeric(parameters["gamma"])

    dx_dt <- alpha * x - beta * x * y
    dy_dt <- delta * x * y - gamma * y
    list(c(dx_dt, dy_dt))
}


# time sequence (1821-1914 is 93 years)
times <- seq(0, 93, by = 0.1)

# Initial condition [prey, predator] acc. to lv_data
state <- c(x = 30, y = 4)

# Negative log likelihood function 
nll <- function(pars) {
    parameters <- c(
        alpha = pars[1],
        beta = pars[2],
        delta = pars[3],
        gamma = pars[4]
    )
    
    # Solve ODE with given parameters
    lv_sol <- try(ode(y = state, times = times, func = lv_ode, 
                  parms = parameters, method = "rk4"), silent = TRUE)
    if (inherits(lv_sol, "try-error")) return(1e9)
    lv_sol_df <- as.data.frame(lv_sol)
    colnames(lv_sol_df) <- c("Time", "Hare", "Lynx")
    # Penalize if ODE produced non-finite values
    if (any(!is.finite(lv_sol_df$Hare)) || any(!is.finite(lv_sol_df$Lynx))) return(1e9)
    
    # interpolate model predictions to observed data time points
    observed_times <- lv_data$Year - lv_data$Year[1]
    hare_pred <- approx(lv_sol_df$Time, lv_sol_df$Hare, 
                        xout = observed_times)$y
    lynx_pred <- approx(lv_sol_df$Time, lv_sol_df$Lynx, 
                        xout = observed_times)$y
    
    # Remove NAs i.e. deal with edge cases
    valid_idx <- !is.na(hare_pred) & !is.na(lynx_pred)
    # Penalize if too few points match
    if (sum(valid_idx) < 0.8 * length(observed_times)) return(1e9)
    hare_obs <- lv_data$Hare[valid_idx]
    lynx_obs <- lv_data$Lynx[valid_idx]
    hare_pred <- hare_pred[valid_idx]
    lynx_pred <- lynx_pred[valid_idx]
    
    # standard deviations
    sigma_hare <- sd(lv_data$Hare, na.rm = TRUE) 
    sigma_lynx <- sd(lv_data$Lynx, na.rm = TRUE) 
    
    # Negative log-likelihood (Gaussian errors)
    nll_hare <- -sum(dnorm(hare_obs, mean = hare_pred, sd = sigma_hare, log = TRUE))
    nll_lynx <- -sum(dnorm(lynx_obs, mean = lynx_pred, sd = sigma_lynx, log = TRUE))
    
    return(nll_hare + nll_lynx)
}

# Fit parameters using optim with bounds for stability
set.seed(123)
initial_pars <- pmax(abs(rnorm(n = 4, mean = c(0.5, 0.03, 0.03, 0.8), sd = 0.2)), 1e-3)
fit <- optim(
    par = initial_pars,
    fn = nll,
    method = "L-BFGS-B",
    lower = c(1e-4, 1e-4, 1e-4, 1e-4),
    upper = c(2, 1, 1, 2),
    control = list(maxit = 10000, factr = 1e7)
)

fit$convergence # should return zero i.e. successful completion
fit$value # neg. log likelihood value

# Extract fitted parameters
fitted_pars <- fit$par
names(fitted_pars) <- c("alpha", "beta", "delta", "gamma")
print(fitted_pars)


# Solve ODE with fitted parameters
parameters <- c(
    alpha = fitted_pars[1],
    beta = fitted_pars[2],
    delta = fitted_pars[3],
    gamma = fitted_pars[4]
)

lv_solution <- ode(y = state, times = times, func = lv_ode, parms = parameters, method = "rk4")
lv_solution_df <- as.data.frame(lv_solution)
colnames(lv_solution_df) <- c("Time", "Hare", "Lynx")
lv_solution_df$Year <- lv_data$Year[1] + lv_solution_df$Time
lv_solution_df <- subset(lv_solution_df, is.finite(Hare) & is.finite(Lynx))

# Lotka-Volterra model fit
svg("./tmapper_blog/images/lv_solution2.svg", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.5, 0.8, 0))

# Hare subplot
plot(lv_solution_df$Year, lv_solution_df$Hare,
     type = "l",
     col = "#85C1E2",
     lwd = 2.5,
     xlab = "Year",
     ylab = "Number (x10^3)",
     main = "Hare Population",
     ylim = c(0, max(lv_data$Hare, lv_solution_df$Hare, na.rm = TRUE) * 1.1),
     cex.main = 1.3,
     cex.axis = 1.0,
     cex.lab = 1.0)
points(lv_data$Year, lv_data$Hare,
       col = "#2E86C1",
       pch = 16,
       cex = 1.2)
grid(nx = NULL, ny = NULL, col = "gray85", lty = "solid", lwd = 0.5)
legend("topleft",
       legend = c("Observed", "Fitted"),
       col = c("#2E86C1", "#85C1E2"),
       lwd = c(NA, 2.5),
       pch = c(16, NA),
       cex = 0.95,
       bty = "n")

# Lynx subplot
plot(lv_solution_df$Year, lv_solution_df$Lynx,
     type = "l",
     col = "#F5B7B1",
     lwd = 2.5,
     xlab = "Year",
     ylab = "Number (x10^3)",
     main = "Lynx Population",
     ylim = c(0, max(lv_data$Lynx, lv_solution_df$Lynx, na.rm = TRUE) * 1.1),
     cex.main = 1.3,
     cex.axis = 1.0,
     cex.lab = 1.0)
points(lv_data$Year, lv_data$Lynx,
       col = "#C0392B",
       pch = 16,
       cex = 1.2)
grid(nx = NULL, ny = NULL, col = "gray85", lty = "solid", lwd = 0.5)
legend("topleft",
       legend = c("Observed", "Fitted"),
       col = c("#C0392B", "#F5B7B1"),
       lwd = c(NA, 2.5),
       pch = c(16, NA),
       cex = 0.95,
       bty = "n")

mtext("Lotka-Volterra Model Fit: Observed vs MLE-Fitted Model (1821-1914)",
      outer = TRUE, line = -1.5, cex = 1.2, font = 2)

dev.off()



# Sample parameters (example)
tutorial_pars <- c(alpha = 1, beta = 1, delta = 1, gamma = 1)

# Create combined plot
svg("./tmapper_blog/images/lv_phase_plane.svg", width = 16, height = 8)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 2), mgp = c(2.5, 0.8, 0))


# Velocity field with trajectories
# Define range for the phase plane
hare_range <- seq(0.1, 5, length.out = 20)
lynx_range <- seq(0.1, 5, length.out = 20)

# Create meshgrid
hare_mesh <- rep(hare_range, length(lynx_range))
lynx_mesh <- rep(lynx_range, each = length(hare_range))

# calculate derivatives at each point 
dx <- tutorial_pars["alpha"] * hare_mesh - tutorial_pars["beta"] * hare_mesh * lynx_mesh
dy <- tutorial_pars["delta"] * hare_mesh * lynx_mesh - tutorial_pars["gamma"] * lynx_mesh

# normalize arrows for visualization
magnitude <- sqrt(dx^2 + dy^2)
dx_norm <- dx / (magnitude + 1e-6) * 0.15
dy_norm <- dy / (magnitude + 1e-6) * 0.15


plot(NULL,
     xlim = c(0, 5),
     ylim = c(0, 5),
     xlab = "Hare Population (x)",
     ylab = "Lynx Population (y)",
     main = "Velocity Field with Trajectories",
     cex.main = 1.2,
     cex.axis = 1.0,
     cex.lab = 1.0)

# add arrows
arrows(hare_mesh - dx_norm / 2, lynx_mesh - dy_norm / 2,
       hare_mesh + dx_norm / 2, lynx_mesh + dy_norm / 2,
       length = 0.08, col = "gray60", lwd = 1.5)

# multiple trajectories from different initial conditions
initial_conditions <- rbind(
  c(2, 2),
  c(0.5, 0.5),
  c(0.5, 1.5),
  c(1.5, 0.5),
  c(3, 3)
)

for (i in 1:nrow(initial_conditions)) {
  # Solve ODE from this initial condition
  ic <- initial_conditions[i, ]
  names(ic) <- c("x", "y")
  
  traj <- ode(y = ic, times = seq(0, 30, by = 0.1), 
              func = lv_ode, parms = tutorial_pars, method = "rk4")
  traj_df <- as.data.frame(traj)
  colnames(traj_df) <- c("Time", "Hare", "Lynx")
  
  # trajectory
  lines(traj_df$Hare, traj_df$Lynx, col = "#1F77B4", lwd = 1.8, alpha = 0.7)
  
  # add starting point
  points(ic[1], ic[2], col = "#D62728", pch = 16, cex = 1.2)
}

# add equilibrium point 
eq_point <- c(tutorial_pars["gamma"] / tutorial_pars["delta"],
              tutorial_pars["alpha"] / tutorial_pars["beta"])
points(eq_point[1], eq_point[2], col = "darkgreen", pch = 8, cex = 2.5, lwd = 2)

grid(nx = NULL, ny = NULL, col = "gray85", lty = "solid", lwd = 0.5)
legend("topright",
       legend = c("Trajectories", "Initial condition", "Equilibrium (1,1)"),
       col = c("#1F77B4", "#D62728", "darkgreen"),
       lwd = c(1.8, NA, NA),
       pch = c(NA, 16, 8),
       cex = 0.8,
       bty = "n")


# Nullclines

plot(NULL,
     xlim = c(0, 3),
     ylim = c(0, 3),
     xlab = "Hare Population (x)",
     ylab = "Lynx Population (y)",
     main = "Nullclines",
     cex.main = 1.2,
     cex.axis = 1.0,
     cex.lab = 1.0)

# x-nullcline
abline(h = 1, col = "#1F77B4", lwd = 2.5, lty = 1)  # y = 1

# y-nullcline
abline(v = 1, col = "#C0392B", lwd = 2.5, lty = 1)  # x = 1

# Add equilibrium point at (1,1)
points(1, 1, col = "darkred", pch = 8, cex = 2.5, lwd = 2.5)

# Add velocity field 
hare_range_small <- seq(0.1, 2.9, length.out = 15)
lynx_range_small <- seq(0.1, 2.9, length.out = 15)

hare_mesh <- rep(hare_range_small, length(lynx_range_small))
lynx_mesh <- rep(lynx_range_small, each = length(hare_range_small))

dx <- tutorial_pars["alpha"] * hare_mesh - tutorial_pars["beta"] * hare_mesh * lynx_mesh
dy <- tutorial_pars["delta"] * hare_mesh * lynx_mesh - tutorial_pars["gamma"] * lynx_mesh

magnitude <- sqrt(dx^2 + dy^2)
dx_norm <- dx / (magnitude + 1e-6) * 0.08
dy_norm <- dy / (magnitude + 1e-6) * 0.08

arrows(hare_mesh - dx_norm / 2, lynx_mesh - dy_norm / 2,
       hare_mesh + dx_norm / 2, lynx_mesh + dy_norm / 2,
       length = 0.05, col = "gray70", lwd = 1)

# Add sample trajectories in positive quadrant
initial_conditions_pos <- rbind(
  c(0.5, 0.5),
  c(1.5, 0.5),
  c(0.5, 1.5),
  c(1, 1.5)
)

for (i in 1:nrow(initial_conditions_pos)) {
  ic <- initial_conditions_pos[i, ]
  names(ic) <- c("x", "y")
  
  traj <- ode(y = ic, times = seq(0, 20, by = 0.1), 
              func = lv_ode, parms = tutorial_pars, method = "rk4")
  traj_df <- as.data.frame(traj)
  colnames(traj_df) <- c("Time", "Hare", "Lynx")
  
  # Plot full trajectory
  lines(traj_df$Hare, traj_df$Lynx, 
        col = "#FF7F0E", lwd = 1.5, alpha = 0.8)
}

grid(nx = NULL, ny = NULL, col = "gray85", lty = "solid", lwd = 0.5)
legend("topright",
       legend = c("x-nullcline (y=1)",
                  "y-nullcline (x=1)",
                  "Center point (1,1)",
                  "Sample trajectories"),
       col = c("#1F77B4", "#C0392B", "darkred", "#FF7F0E"),
       lwd = c(2.5, 2.5, NA, 1.5),
       pch = c(NA, NA, 8, NA),
       cex = 0.8,
       bty = "n")

mtext("Lotka-Volterra Phase Plane Analysis (α=β=γ=δ=1)",
      outer = TRUE, line = 0, cex = 1.3, font = 2)

dev.off()