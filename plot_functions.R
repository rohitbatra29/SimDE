# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# This file contains the plot functions needed for the plot outputs in Server

# Load packages

library(deSolve)
library(phaseR)
library(Sim.DiffProc)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(purrr)
library(ggpubr)
library(tidyr)
library(JuliaCall)
library(diffeqr)

# install.packages(c("deSolve", "dplyr", "ggplot2", "tidyr", "phaseR", "Sim.DiffProc", "ggh4x", "purrr", "ggpubr"))
# Julia's Path ------------------------------------------------------------

path <- "/path/to/your/julia-1.10.4/bin"

# Univariate Plots --------------------------------------------------------

univariate_plot <- function(plot_data){
  
  k <- ggplot(data = plot_data) + 
    geom_line(aes(x = time, y = Y)) +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
          axis.title = element_text(size = 25, face = "bold.italic", color = "#0044A4"),
          axis.text = element_text(color = "#0044A4", size = 20),
          panel.border = element_blank(),
          axis.line = element_line(colour = "#0044A4", linewidth = 1),
          axis.ticks.length = unit(.3, "cm"),
          axis.ticks = element_line(color = "#0044A4")
    ) +
    labs(x = "Time", y = "Y") +
    scale_x_continuous(limits = c(first(plot_data$time), last(plot_data$time)))
  
  return(k)
  
}


## 1. Linear ODE ---------------------------------------------------------

plot_ode <- function(input){
  
  beta <- input$beta 
  mean_y <- input$mean
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  start_value <- input$start_value
  
  # output calculation
  parameters <- c(a = beta,
                  c = mean_y)
  
  state <- c(Y = start_value)
  
  linear_ode <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      # differentials
      dY = a*(Y-c)
      # return list
      list(dY)
    })
  }
  
  times <- seq(period_start, period_end, by = step_size)
  
  out <- ode(y = state, times = times, func = linear_ode, parms = parameters)
  
  plot_y <- univariate_plot(data.frame(out))
  
  return(plot_y)

}

## 2. Polynomial ODE -----------------------------------------------------

plot_pde <- function(input){
  
  beta1 <- input$beta1
  beta2 <- input$beta2
  beta3 <- input$beta3
  beta4 <- input$beta4
  mean_y <- input$mean
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  start_value <- input$start_value
  
  parameters <- c(b1 = beta1,
                  b2 = beta2,
                  b3 = beta3,
                  b4 = beta4,
                  c = mean_y)
  
  state <- c(Y = start_value)
  
  poly_ode <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      # differentials
      dY = b1*Y + b2*(Y^2) + b3*(Y^3) + b4*(Y^4) - c
      # return list
      list(dY)
    })
  }
  
  times <- seq(period_start, period_end, by = step_size)
  
  out_poly <- ode(y = state, times = times, func = poly_ode, parms = parameters)
  
  plot_y <- univariate_plot(data.frame(out_poly))
  
  return(plot_y)
}

## 3. Linear SDE (OU) -----------------------------------------------------------------

plot_ou <- function(input){
  
  beta <- input$beta
  mean_y <- input$mean
  sigma <- input$sigma
  seed_y <- input$seed
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  start_value <- input$start_value
  
  set.seed(seed_y)
  
  drift_expr <- substitute(b * (x - m), list(b = beta, m = mean_y))
  diff_expr <- substitute(sig, list(sig = sigma))
  
  ou_mod <- snssde1d(drift = as.expression(drift_expr),
                     diffusion = as.expression(diff_expr), M = 1,
                     t0 = period_start, T = period_end, x0 = start_value,
                     N = floor((period_end - period_start) / step_size) + 1,
                     Dt = step_size,
                     type = "str",
                     method = "rk1")
  
  times <- seq(period_start, period_end, by = step_size)
  
  plot_data <- data.frame("time" = times, "Y" = ou_mod$X)
  
  plot_y <- univariate_plot(plot_data)
  
  return(plot_y)
  
}

## 4. Polynomial SDE Plot -----------------------------------------------------

plot_p_sde <- function(input){
  
  b1 <- input$beta1
  b2 <- input$beta2
  b3 <- input$beta3
  b4 <- input$beta4
  mean_y <- input$mean
  sigma <- input$sigma
  
  seed_y <- input$seed
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  start_value <- input$start_value
  
  set.seed(seed_y)
  
  drift_expr <- substitute(b1*x + b2*(x^2) + b3*(x^3) + b4*(x^4) - m, list(b1 = b1, b2 = b2, b3 = b3, b4 = b4, m = mean_y))
  diff_expr <- substitute(sig, list(sig = sigma))
  
  poly_sde_mod <- snssde1d(drift = as.expression(drift_expr),
                           diffusion = as.expression(diff_expr), M = 1,
                           t0 = period_start, T = period_end, x0 = start_value,
                           N = floor((period_end - period_start) / step_size) + 1,
                           Dt = step_size,
                           type = "str", # does it need to be stratonovich?
                           method = "rk1")
  
  times <- seq(period_start, period_end, by = step_size)
  plot_data <- data.frame("time" = times, "Y" = poly_sde_mod$X)
  
  plot_y <- univariate_plot(plot_data)
  
  return(plot_y)
}


## 5. Harmonic ODE Plot -------------------------------------------------------

plot_harmonic_ode <- function(input){
  
  omega <- input$omega
  mu <- input$mean
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  start_value <- input$start_value
  
  parameters <- c(w = omega,
                  m = mu)
  
  state <- c(Y = start_value,
             P = 0)
  
  # P is the additional variable to reduce the order of Y.
  
  harmonic_ode <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      # differentials
      dY = P
      dP = w*(Y - m)
      # return list
      list(c(dY, dP))
    })
  }
  
  times <- seq(period_start, period_end, by = step_size)
  
  out <- ode(y = state, times = times, func = harmonic_ode, parms = parameters)
  
  plot_y <- univariate_plot(as.data.frame(out))
  
  return(plot_y)
}


## 6. Damped Oscillator ODE ---------------------------------------------------

plot_damp_ode <- function(input){
  
  omega <- input$omega
  zeta <- input$zeta
  mu <- input$mean
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  start_value <- input$start_value
  
  parameters <- c(w = omega,
                  m = mu,
                  z = zeta)
  
  state <- c(Y = start_value,
             P = 0)
  
  harmonic_ode <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      # differentials
      dY = P
      dP = w*(Y - m) + z*P
      # return list
      list(c(dY, dP))
    })
  }
  
  times <- seq(period_start, period_end, by = step_size)
  
  out <- ode(y = state, times = times, func = harmonic_ode, parms = parameters)
  
  plot_y <- univariate_plot(as.data.frame(out))
  
  return(plot_y)
}


## 7. Harmonic (Undamped) Oscillator SDE --------------------------------------

plot_harmonic_sde <- function(input){
  
  omega_y <- input$omega
  mu_y <- input$mean
  sigma_y <- input$sigma
  seed_y <- input$seed
  start_value <- input$start_value
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  
  set.seed(seed_y)
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  f <- function(u,p,t) {
    omega_y <- p[1]
    mu_y <- p[2]
    
    du1 = u[2]
    du2 = omega_y*(u[1] - mu_y)
    return(c(du1, du2))
  }
  g <- function(u,p,t) {
    return(c(0, sigma_y))
  }
  
  u0 <- c(start_value, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_y)
  sol <- de$solve(prob, saveat = steps)
  
  plot_data <- data.frame("time" = sol$t, as.data.frame(t(sapply(sol$u, identity))))
  
  head(plot_data)
  
  colnames(plot_data) <- c("time", "Y", "P")
  
  plot_y <- univariate_plot(plot_data)
  
  return(plot_y)
}


## 8. Damped Oscillator SDE ---------------------------------------------------

plot_damp_sde <- function(input){
  
  omega_y <- input$omega
  zeta_y <- input$zeta
  mu_y <- input$mean
  sigma_y <- input$sigma
  seed_y <- input$seed
  start_value <- input$start_value
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  f <- function(u,p,t) {
    omega_y <- p[1]
    zeta_y <- p[2]
    mu_y <- p[3]
    
    du1 = u[2]
    du2 = omega_y*(u[1] - mu_y) + zeta_y*u[2]
    return(c(du1, du2))
  }
  g <- function(u,p,t) {
    return(c(0, sigma_y))
  }
  
  u0 <- c(start_value, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_y, zeta_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_y)
  sol <- de$solve(prob, saveat = steps)
  
  plot_data <- data.frame("time" = sol$t, as.data.frame(t(sapply(sol$u, identity))))
  
  head(plot_data)
  
  colnames(plot_data) <- c("time", "Y", "P")
  
  plot_y <- univariate_plot(plot_data)
  
  return(plot_y)
}

# Bivariate Plots ---------------------------------------------------------

bivariate_plot <- function(plot_data){
  k <- ggplot(data = plot_data) + 
    geom_line(aes(x = time, y = value, color = variable, linetype = variable)) +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
          axis.title = element_text(size = 25, face = "bold.italic", color = "#0044A4"),
          axis.text = element_text(color = "#0044A4", size = 20),
          panel.border = element_blank(),
          axis.line = element_line(colour = "#0044A4", linewidth = 1),
          axis.ticks.length = unit(.3, "cm"),
          axis.ticks = element_line(color = "#0044A4")
    ) +
    labs(x = "Time", y = "X and Y") +
    scale_x_continuous(limits = c(first(plot_data$time), last(plot_data$time))) + 
    scale_color_manual(values = c("black", "blue"), name = "Variables") +
    scale_linetype_manual(values = c("solid", "longdash"), name = "Variables")
  
  return(k)
}


## 1. Both are Linear ODEs -------------------------------------------------

plot_xy_ode2 <- function(input){

  
  parameters <- c(beta_x = input$beta_x,
                  gamma_x = input$gamma_x,
                  beta_y = input$beta_y,
                  gamma_y = input$gamma_y,
                  mu_x = input$mean_x,
                  mu_y = input$mean_y) 
  
  state <- c(X = input$start_x,
             Y = input$start_y)
  
  if(input$xy_coupling == "coup_diff"){
    linear_ode2 <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = beta_x * (X - mu_x) + gamma_x * (X - Y)
        dY = beta_y * (Y - mu_y) + gamma_y * (Y - X)
        # return list
        list(c(dX, dY))
      })
    }
  } else if(input$xy_coupling == "coup_mean"){
    linear_ode2 <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = beta_x * (X - mu_x) + gamma_x * (Y - mu_y)
        dY = beta_y * (Y - mu_y) + gamma_y * (X - mu_x)
        # return list
        list(c(dX, dY))
      })
    }
  }
  
  
  times <- seq(input$period[1], input$period[2], by = input$step_size) 
  
  # LSODA
  l_ode2 <- ode(y = state, times = times, func = linear_ode2, parms = parameters)
  
  plot_data <- as.data.frame(l_ode2) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(time),
           value = as.numeric(value),
           variable = as.factor(variable))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


## 2. Linear SDE and Linear ODE -----------------------------------------------
# For simulation, we need v1/x as a Linear SDE and v2/y as a linear ODE. This might be different than the user specified x and y .

plot_xy_l2_ode_sde <- function(input){
  
  set.seed(input$seed)
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size

  
  if (input$modeltype_x == "ou"){
    
    # User specified X as Linear SDE and Y as Linear ODE
    
    beta_x <- input$beta_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    sigma_x <- input$sigma_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "ou"){
    
    # User specified Y as Linear SDE and X as Linear ODE but the simulation codes it with flipped notation.
    beta_x <- input$beta_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    sigma_x <- input$sigma_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  # print(flipped)
  
  if(input$xy_coupling == "coup_diff"){
    
    drift_expr_x <- substitute(beta_x*(x - mu_x) + gamma_x * (x - y), list(beta_x = beta_x, mu_x = mu_x, 
                                                                           gamma_x = gamma_x)) 
    drift_expr_y <- substitute(beta_y*(y - mu_y) + gamma_y * (y - x), list(beta_y = beta_y, mu_y = mu_y,
                                                                           gamma_y = gamma_y))
    
  } else if(input$xy_coupling == "coup_mean"){
    
    drift_expr_x <- substitute(beta_x*(x - mu_x) + gamma_x * (y - mu_y), list(beta_x = beta_x, mu_x = mu_x, 
                                                                           gamma_x = gamma_x, mu_y = mu_y)) 
    drift_expr_y <- substitute(beta_y*(y - mu_y) + gamma_y * (x - mu_x), list(beta_y = beta_y, mu_y = mu_y,
                                                                           gamma_y = gamma_y, mu_x = mu_x))
    
  }
  
  diff_expr_x <- substitute(sig, list(sig = sigma_x))
  
  # For simulation, x is the Linear SDE and y is the Linear ODE
  l2_ode_sde <- snssde2d(drift = c(as.expression(drift_expr_x), as.expression(drift_expr_y)),
                         diffusion = c(as.expression(diff_expr_x), expression(0)),
                         t0 = period_start, T = period_end, x0 = c(start_x, start_y),
                         N = (period_end - period_start) / step_size,
                         Dt = step_size,
                         type = "str", # does it need to be stratonovich?
                         method = "rk1")
  
  times <- seq(period_start, period_end, by = step_size)
  sim_data <- data.frame("time" = times, "X" = l2_ode_sde$X, "Y" = l2_ode_sde$Y)
  
  if(flipped == TRUE){
    colnames(sim_data) <- c("time", "Y", "X")
  }

  plot_data <- as_tibble(sim_data) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable")
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


## 3. Both are Linear SDEs ----------------------------------------------------

plot_xy_l2_ode_sde <- function(input){
  
  set.seed(input$seed)
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  
  beta_x <- input$beta_x
  gamma_x <- input$gamma_x
  start_x <- input$start_x 
  mu_x <- input$mean_x
  sigma_x <- input$sigma_x
  
  beta_y <- input$beta_y
  gamma_y <- input$gamma_y
  start_y <- input$start_y
  mu_y <- input$mean_y
  sigma_y <- input$sigma_y
  
  if(input$xy_coupling == "coup_diff"){
    
    drift_expr_x <- substitute(beta_x*(x - mu_x) + gamma_x * (x - y), list(beta_x = beta_x, mu_x = mu_x, 
                                                                           gamma_x = gamma_x)) 
    drift_expr_y <- substitute(beta_y*(y - mu_y) + gamma_y * (y - x), list(beta_y = beta_y, mu_y = mu_y,
                                                                           gamma_y = gamma_y))
    
  } else if(input$xy_coupling == "coup_mean"){
    
    drift_expr_x <- substitute(beta_x*(x - mu_x) + gamma_x * (y - mu_y), list(beta_x = beta_x, mu_x = mu_x, 
                                                                           gamma_x = gamma_x, mu_y = mu_y)) 
    drift_expr_y <- substitute(beta_y*(y - mu_y) + gamma_y * (x - mu_x), list(beta_y = beta_y, mu_y = mu_y,
                                                                           gamma_y = gamma_y, mu_x = mu_x))
    
  }
  
  diff_expr_x <- substitute(sig_x, list(sig_x = sigma_x))
  diff_expr_y <- substitute(sig_y, list(sig_y = sigma_y))
  
  # Both X and Y are SDE
  
  l_sde2 <- snssde2d(drift = c(as.expression(drift_expr_x), as.expression(drift_expr_y)),
                     diffusion = c(as.expression(diff_expr_x), as.expression(diff_expr_y)),
                     t0 = period_start, T = period_end, x0 = c(start_x, start_y),
                     N = (period_end - period_start) / step_size,
                     Dt = step_size,
                     type = "str", # does it need to be stratonovich?
                     method = "rk1")
  
  times <- seq(period_start, period_end, by = step_size)
  sim_data <- data.frame("time" = times, "X" = l_sde2$X, "Y" = l_sde2$Y)
  plot_data <- as_tibble(sim_data) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable")
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 4. Linear ODE and Harmonic Oscillator ODE -------------------------------
# For simulation, v1/x is a Linear ODE and v2/y is the harmonic oscillator ODE. This might be different than the user specified x and y.

plot_xy_l_harmonic_odes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size 
  

  if (input$modeltype_x == "l_ode"){
    
    # User specified X as Linear ODE and Y as Harmonic Oscillator
    beta_x <- input$beta_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    
    omega_y <- input$omega_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "l_ode"){
    
    # User specified Y as Linear ODE and X as Harmonic Oscillator but the simulation codes it with flipped notation.
    beta_x <- input$beta_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    
    omega_y <- input$omega_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  parameters <- c(beta_x = beta_x,
                  gamma_x = gamma_x,
                  mu_x = mu_x,
                  omega_y = omega_y,
                  gamma_y = gamma_y,
                  mu_y = mu_y) 
  
  state <- c(X = start_x,
             Y = start_y,
             Z = 0)
  
  #print(flipped)
  
  if(input$xy_coupling == "coup_diff"){
    
    l_ode_harmonic_mod <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = beta_x * (X - mu_x) + gamma_x * (X - Y)
        dY = Z
        dZ = omega_y * (Y - mu_y) + gamma_y * (Y - X)
        # return list
        list(c(dX, dY, dZ))
      })
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    l_ode_harmonic_mod <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = beta_x * (X - mu_x) + gamma_x * (Y - mu_y)
        dY = Z
        dZ = omega_y * (Y - mu_y) + gamma_y * (X - mu_x)
        # return list
        list(c(dX, dY, dZ))
      })
    }
    
  }
  
  times <- seq(period_start, period_end, by = step_size) 
  
  # LSODA
  l_ode_harmonic <- ode(y = state, times = times, func = l_ode_harmonic_mod, parms = parameters)
  
  if(flipped == TRUE){
    colnames(l_ode_harmonic) = c("time", "Y", "X", "Z")
  }
  
  plot_data <- as_tibble(l_ode_harmonic) %>%
    select(-Z) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(time),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 5. Harmonic Oscillator ODE and Linear SDE -------------------------------
# For simulation, X/U1 is a Harmonic SDE and Y/U3 is the Linear SDE. This might be different than the user specified x and y.

plot_xy_l_sde_harmonic <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "harmonic_ode"){
    
    # User specified X as Harmonic ODE and Y as Linear SDE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    sigma_y <- input$sigma_y
    
    
    flipped = FALSE # X is being simulated as user X and Y is being simulated as user Y
    
  } else if (input$modeltype_y == "harmonic_ode"){
    
    # User specified Y as Harmonic ODE and X as Linear SDE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    sigma_y <- input$sigma_x
    
    flipped = TRUE # Y is being simulated as user X and X is being simulated as user Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      
      beta_y <- p[4]
      gamma_y <- p[5]
      mu_y <- p[6]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[1] - u[3])
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3))
    }
    
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      
      beta_y <- p[4]
      gamma_y <- p[5]
      mu_y <- p[6]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[3] - mu_y)
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[1] -  mu_x)
      return(c(du1, du2, du3))
    }
    
  }
  
  #Since only Y is stochastic, we feed the sigma_y to U3
  g <- function(u,p,t) {
    return(c(0, 0, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x,
         beta_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-P) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)

}


## 6. Both are Harmonic Oscillators ODEs -----------------------------------

plot_xy_harmonics_ode <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size 
  
  omega_x <- input$omega_x
  gamma_x <- input$gamma_x
  start_x <- input$start_x
  mu_x <- input$mean_x
  
  omega_y <- input$omega_y
  gamma_y <- input$gamma_y
  start_y <- input$start_y
  mu_y <- input$mean_y
  
  # output calculation
  parameters <- c(omega_x = omega_x,
                  gamma_x = gamma_x,
                  mu_x = mu_x,
                  omega_y = omega_y,
                  gamma_y = gamma_y,
                  mu_y = mu_y) 
  
  state <- c(X = start_x,
             Y = start_y,
             P = 0, # p is used to reduce the order of x
             Q = 0) # q is used to reduce the order of y
  
  if(input$xy_coupling == "coup_diff"){
    
    harmonics_ode_mod <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = P
        dP = omega_x * (X - mu_x) + gamma_x * (X - Y)
        dY = Q
        dQ = omega_y * (Y - mu_y) + gamma_y * (Y - X)
        # return list
        list(c(dX, dY, dP, dQ))
      })
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    harmonics_ode_mod <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = P
        dP = omega_x * (X - mu_x) + gamma_x * (Y - mu_y)
        dY = Q
        dQ = omega_y * (Y - mu_y) + gamma_y * (X - mu_x)
        # return list
        list(c(dX, dY, dP, dQ))
      })
    }
    
  }
  
  times <- seq(period_start, period_end, by = step_size) 
  
  harmonics_ode <- ode(y = state, times = times, func = harmonics_ode_mod, parms = parameters)
  
  
  plot_data <- as_tibble(harmonics_ode) %>%
    select(-P, -Q) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(time),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 7. Harmonic Oscillator SDE and Linear ODE -------------------------------
# For simulation, v1/x is a Harmonic SDE and v2/y is the Linear ODE. This might be different than the user specified x and y.

plot_xy_harmonic_sde_l_ode <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "harmonic_sde"){
    
    # User specified X as Harmonic SDE and Y as Linear ODE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    sigma_x <- input$sigma_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    flipped = FALSE # X is being simulated as user X and Y is being simulated as user Y
    
  } else if (input$modeltype_y == "harmonic_sde"){
    
    # User specified Y as Harmonic SDE and X as Linear ODE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    sigma_x <- input$sigma_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as user X and X is being simulated as user Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      
      beta_y <- p[4]
      gamma_y <- p[5]
      mu_y <- p[6]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[1] - u[3])
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      
      beta_y <- p[4]
      gamma_y <- p[5]
      mu_y <- p[6]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[3] - mu_y)
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[1] -  mu_x)
      return(c(du1, du2, du3))
    }
  
  }
  
  # Since only x is stochastic, we feed the sigma_x to U2
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0))
  }
  
  u0 <- c(start_x, 0, start_y)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x,
         beta_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-P) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 8. Harmonic Oscillator SDE and Linear SDE -------------------------------
# For simulation, X/U1 is a Harmonic SDE and Y/U3 is the Linear SDE. This might be different than the user specified x and y.

plot_xy_harmonic_l_sdes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "harmonic_sde"){
    
    # User specified X as Harmonic SDE and Y as Linear SDE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    sigma_x <- input$sigma_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    sigma_y <- input$sigma_y
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "harmonic_sde"){
    
    # User specified Y as Harmonic SDE and X as Linear SDE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    sigma_x <- input$sigma_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    sigma_y <- input$sigma_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      
      beta_y <- p[4]
      gamma_y <- p[5]
      mu_y <- p[6]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[1] - u[3])
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      
      beta_y <- p[4]
      gamma_y <- p[5]
      mu_y <- p[6]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[3] - mu_y)
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[1] -  mu_x)
      return(c(du1, du2, du3))
    }
    
  }
  
  # Since both x and y are stochastic, we feed the sigma_x to U2 and sigma_y to U3
  g <- function(u,p,t) {
    return(c(0, sigma_x, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x,
         beta_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-P) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 9. Damped Oscillator ODE and Linear ODE ---------------------------------
# For simulation, v1/x is a Damped Oscillator ODE and v2/y is the Linear ODE. This might be different than the user specified x and y.

plot_xy_damped_linear_odes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size 
  
  
  if (input$modeltype_x == "damp_ode"){
    
    # User specified X as Damped Oscillator ODE and Y as Linear ODE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "damp_ode"){
    
    # User specified Y as Damped Oscillator ODE and X as Linear ODE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  parameters <- c(omega_x = omega_x,
                  gamma_x = gamma_x,
                  mu_x = mu_x,
                  zeta_x = zeta_x,
                  beta_y = beta_y,
                  gamma_y = gamma_y,
                  mu_y = mu_y) 
  
  state <- c(X = start_x,
             Y = start_y,
             Z = 0)
  
  if(input$xy_coupling == "coup_diff"){
    
    damped_linear_odes_mod <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = Z
        dZ = omega_x * (X - mu_x) + gamma_x * (X - Y) + zeta_x * Z
        dY = beta_y * (Y - mu_y) + gamma_y * (Y - X)
        # return list
        list(c(dX, dY, dZ))
      })
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    damped_linear_odes_mod <- function(t, state, parameters){
      with(as.list(c(state, parameters)),{
        # differentials
        dX = Z
        dZ = omega_x * (X - mu_x) + gamma_x * (Y - mu_y) + zeta_x * Z
        dY = beta_y * (Y - mu_y) + gamma_y * (X - mu_x)
        # return list
        list(c(dX, dY, dZ))
      })
    }
    
  }
  
  times <- seq(period_start, period_end, by = step_size) 
  
  # LSODA
  damped_linear_odes <- ode(y = state, times = times, func = damped_linear_odes_mod, parms = parameters)
  
  if(flipped == TRUE){
    colnames(damped_linear_odes) = c("time", "Y", "X", "Z")
  }
  
  plot_data <- as_tibble(damped_linear_odes) %>%
    select(-Z) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(time),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 10. Damped Oscillator ODE and Linear SDE --------------------------------
# For simulation, v1/x is a Damped Oscillator ODE and v2/y is the Linear SDE. This might be different than the user specified x and y.

plot_xy_l_sde_damped <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_ode"){
    
    # User specified X as Damped Oscillator ODE and Y as Linear SDE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    sigma_y <- input$sigma_y
    
    
    flipped = FALSE # X is being simulated as user X and Y is being simulated as user Y
    
  } else if (input$modeltype_y == "damp_ode"){
    
    # User specified Y as Damped Oscillator ODE and X as Linear SDE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    sigma_y <- input$sigma_x
    
    flipped = TRUE # Y is being simulated as user X and X is being simulated as user Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      beta_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      beta_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[1] -  mu_x)
      return(c(du1, du2, du3))
    }
    
  }
  
  # Since only Y is stochastic, we feed only sigma_y to U3
  g <- function(u,p,t) {
    return(c(0, 0, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         beta_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-P) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 11. Damped and Harmonic Oscillators ODEs --------------------------------
# For simulation, v1/x is a Damped Oscillator ODE and v2/y is the Harmonic Oscillator ODE. This might be different than the user specified x and y.

plot_xy_damped_harmonic_odes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size 
  
  
  if (input$modeltype_x == "damp_ode"){
    
    # User specified X as Linear ODE and Y as Harmonic Oscillator
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    
    omega_y <- input$omega_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "damp_ode"){
    
    # User specified Y as Linear ODE and X as Harmonic Oscillator but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    
    omega_y <- input$omega_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  parameters <- c(omega_x = omega_x,
                  gamma_x = gamma_x,
                  mu_x = mu_x,
                  zeta_x = zeta_x,
                  omega_y = omega_y,
                  gamma_y = gamma_y,
                  mu_y = mu_y) 
  
  state <- c(X = start_x,
             Y = start_y,
             P = 0,
             Q = 0)
  
  # Here P is for the drift part of X and Q is for the drift part of Y since both X and Y need two equations.
  
  if(input$xy_coupling == "coup_diff"){
    
        damped_harmonic_odes_mod <- function(t, state, parameters){
          with(as.list(c(state, parameters)),{
            # differentials
            dX = P
            dP = omega_x*(X - mu_x) + gamma_x * (X - Y) + zeta_x*P
            dY = Q
            dQ = omega_y * (Y - mu_y) + gamma_y * (Y - X)
            # return list
            list(c(dX, dY, dP, dQ))
          })
        }
    
  } else if(input$xy_coupling == "coup_mean"){
    
        damped_harmonic_odes_mod <- function(t, state, parameters){
          with(as.list(c(state, parameters)),{
            # differentials
            dX = P
            dP = omega_x*(X - mu_x) + gamma_x * (Y - mu_y) + zeta_x*P
            dY = Q
            dQ = omega_y * (Y - mu_y) + gamma_y * (X - mu_x)
            # return list
            list(c(dX, dY, dP, dQ))
          })
        }
    
  }
  
  times <- seq(period_start, period_end, by = step_size) 
  
  # LSODA
  damped_harmonic_odes <- ode(y = state, times = times, func = damped_harmonic_odes_mod, parms = parameters)
  
  if(flipped == TRUE){
    colnames(damped_harmonic_odes) = c("time", "Y", "X", "P", "Q")
  }
  
  plot_data <- as_tibble(damped_harmonic_odes) %>%
    select(-P, -Q) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(time),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


## 12. Both are Damped Oscillator ODEs -------------------------------------

plot_xy_damped_odes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  
  omega_x <- input$omega_x
  gamma_x <- input$gamma_x
  start_x <- input$start_x
  mu_x <- input$mean_x
  zeta_x <- input$zeta_x
  
  omega_y <- input$omega_y
  gamma_y <- input$gamma_y
  start_y <- input$start_y
  mu_y <- input$mean_y
  zeta_y <- input$zeta_y
  
  parameters <- c(omega_x = omega_x,
                  gamma_x = gamma_x,
                  mu_x = mu_x,
                  zeta_x = zeta_x,
                  omega_y = omega_y,
                  gamma_y = gamma_y,
                  mu_y = mu_y,
                  zeta_y = zeta_y) 
  
  state <- c(X = start_x,
             Y = start_y,
             P = 0,
             Q = 0)
  
  # Here P is for the drift part of X and Q is for the drift part of Y since both X and Y need two equations.
  
  if(input$xy_coupling == "coup_diff"){
    
        damped_odes_mod <- function(t, state, parameters){
          with(as.list(c(state, parameters)),{
            # differentials
            dX = P
            dP = omega_x * (X - mu_x) + gamma_x * (X - Y) + zeta_x*P
            dY = Q
            dQ = omega_y * (Y - mu_y) + gamma_y * (Y - X) + zeta_y*Q
            # return list
            list(c(dX, dY, dP, dQ))
          })
        }
    
  } else if(input$xy_coupling == "coup_mean"){
    
        damped_odes_mod <- function(t, state, parameters){
          with(as.list(c(state, parameters)),{
            # differentials
            dX = P
            dP = omega_x * (X - mu_x) + gamma_x * (Y - mu_y) + zeta_x*P
            dY = Q
            dQ = omega_y * (Y - mu_y) + gamma_y * (X - mu_x) + zeta_y*Q
            # return list
            list(c(dX, dY, dP, dQ))
          })
        }
    
  }
  
  times <- seq(period_start, period_end, by = step_size) 
  
  damped_odes <- ode(y = state, times = times, func = damped_odes_mod, parms = parameters)
  
  plot_data <- as_tibble(damped_odes) %>%
    select(-P, -Q) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(time),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


# 13. Damped Oscillator SDE and Linear ODE --------------------------------
# For simulation, v1/x is a Damped Oscillator SDE and v2/y is the Linear ODE. This might be different than the user specified x and y.

plot_xy_damp_sde_l_ode <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_sde"){
    
    # User specified X as Damped Oscillator SDE and Y as Linear ODE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    sigma_x <- input$sigma_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    flipped = FALSE # X is being simulated as user X and Y is being simulated as user Y
    
  } else if (input$modeltype_y == "damp_sde"){
    
    # User specified Y as Damped Oscillator SDE and X as Linear ODE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    sigma_x <- input$sigma_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as user X and X is being simulated as user Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      beta_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      beta_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[1] -  mu_x)
      return(c(du1, du2, du3))
    }
    
  }
  
  # Since only X is stochastic, we feed only sigma_x to U2
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0))
  }
  
  u0 <- c(start_x, 0, start_y)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         beta_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-P) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


# 14. Damped Oscilator SDE and Linear SDE ---------------------------------
# For simulation, v1/x is a Damped Oscillator SDE and v2/y is the Linear SDE. This might be different than the user specified x and y.

plot_xy_damp_sde_l_sde <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_sde"){
    
    # User specified X as Damped Oscillator SDE and Y as Linear SDE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    sigma_x <- input$sigma_x
    
    beta_y <- input$beta_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    sigma_y <- input$sigma_y
    
    flipped = FALSE # X is being simulated as user X and Y is being simulated as user Y
    
  } else if (input$modeltype_y == "damp_sde"){
    
    # User specified Y as Damped Oscillator SDE and X as Linear ODE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    sigma_x <- input$sigma_y
    
    beta_y <- input$beta_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    sigma_y <- input$sigma_x
    
    flipped = TRUE # Y is being simulated as user X and X is being simulated as user Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      beta_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      beta_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = beta_y*(u[3] - mu_y) + gamma_y*(u[1] -  mu_x)
      return(c(du1, du2, du3))
    }
    
  }
  
  # Since both X and Y are stochastic, we feed sigma_x to U2 and sigma_y to U3
  g <- function(u,p,t) {
    return(c(0, sigma_x, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         beta_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-P) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


# 15. Harmonic Oscillator SDE and ODE -------------------------------------
# For simulation, X/U1 is a harmonic SDE and Y/U3 is a harmonic ode. This might be different than the user specified x and y.
plot_xy_harmonic_sde_ode <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "harmonic_sde"){
    
        # User specified X as Harmonic SDE and Y as Harmonic ODE
        omega_x <- input$omega_x
        gamma_x <- input$gamma_x
        start_x <- input$start_x
        mu_x <- input$mean_x
        sigma_x <- input$sigma_x
        
        omega_y <- input$omega_y
        gamma_y <- input$gamma_y
        start_y <- input$start_y
        mu_y <- input$mean_y
        
        
        flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "harmonic_sde"){
    
        # User specified X as Harmonic ODE and Y as Harmonic SDE but the simulation codes it with flipped notation.
        omega_x <- input$omega_y
        gamma_x <- input$gamma_y
        start_x <- input$start_y
        mu_x <- input$mean_y
        sigma_x <- input$sigma_y
        
        omega_y <- input$omega_x
        gamma_y <- input$gamma_x
        start_y <- input$start_x
        mu_y <- input$mean_x
        
        flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
          f <- function(u,p,t) {
            omega_x <- p[1]
            gamma_x <- p[2]
            mu_x <- p[3]
            
            omega_y <- p[4]
            gamma_y <- p[5]
            mu_y <- p[6]
            
            du1 = u[2]
            du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[1] - u[3])
            du3 = u[4]
            du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
            return(c(du1, du2, du3, du4))
          }
    
  } else if(input$xy_coupling == "coup_mean"){
    
        f <- function(u,p,t) {
          omega_x <- p[1]
          gamma_x <- p[2]
          mu_x <- p[3]
          
          omega_y <- p[4]
          gamma_y <- p[5]
          mu_y <- p[6]
          
          du1 = u[2]
          du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[3] - mu_y)
          du3 = u[4]
          du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[1] - mu_x)
          return(c(du1, du2, du3, du4))
        }
    
  }
  
  # Since only x/U1 is harmonic SDE, we feed the sigma_x to U2.
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0, 0))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x,
         omega_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X", "Q")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}

## 16. Both are Harmonic Oscillators SDE ----------------------------------------

plot_xy_harmonic_sdes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  omega_x <- input$omega_x
  gamma_x <- input$gamma_x
  start_x <- input$start_x
  mu_x <- input$mean_x
  sigma_x <- input$sigma_x
  
  omega_y <- input$omega_y
  gamma_y <- input$gamma_y
  start_y <- input$start_y
  mu_y <- input$mean_y
  sigma_y <- input$sigma_y
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
        f <- function(u,p,t) {
          omega_x <- p[1]
          gamma_x <- p[2]
          mu_x <- p[3]
          
          omega_y <- p[4]
          gamma_y <- p[5]
          mu_y <- p[6]
          
          du1 = u[2]
          du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[1] - u[3])
          du3 = u[4]
          du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
          return(c(du1, du2, du3, du4))
        }
    
  } else if(input$xy_coupling == "coup_mean"){
    
        f <- function(u,p,t) {
          omega_x <- p[1]
          gamma_x <- p[2]
          mu_x <- p[3]
          
          omega_y <- p[4]
          gamma_y <- p[5]
          mu_y <- p[6]
          
          du1 = u[2]
          du2 = omega_x*(u[1] - mu_x) + gamma_x*(u[3] - mu_y)
          du3 = u[4]
          du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[1] - mu_x)
          return(c(du1, du2, du3, du4))
        }
    
  }
  
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x,
         omega_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 17. Damped Oscillator ODE and Harmonic Oscillator SDE -----------------
# For simulation, X/U1 is a damped ODE and Y/U3 is a Harmonic SDE. This might be different than the user specified x and y.

plot_xy_damp_ode_harmonic_sde <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_ode"){
    
    # User specified X as Damped ODE and Y as Harmonic SDE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    
    omega_y <- input$omega_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    sigma_y <- input$sigma_y
    
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "damp_ode"){
    
    # User specified X as Harmonic SDE and Y as Damped ODE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    
    omega_y <- input$omega_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    sigma_y <- input$sigma_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3, du4))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[1] - mu_x)
      return(c(du1, du2, du3, du4))
    }
  }
  
  # Since only Y/U3 is harmonic SDE, we feed the sigma_x to U4.
  g <- function(u,p,t) {
    return(c(0, 0, 0, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         omega_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X", "Q")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


## 18. Damped Oscillator SDE and Harmonic Oscillator ODE -------------------
# For simulation, X/U1 is a Damped SDE and Y/U3 is a Harmonic ODE. This might be different than the user specified x and y.

plot_xy_damp_sde_harmonic_ode <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_sde"){
    
    # User specified X as Damped SDE and Y as Harmonic ODE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    sigma_x <- input$sigma_x
    
    omega_y <- input$omega_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "damp_sde"){
    
    # User specified X as Harmonic ODE and Y as Damped SDE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    sigma_x <- input$sigma_y
    
    omega_y <- input$omega_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3, du4))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[1] - mu_x)
      return(c(du1, du2, du3, du4))
    }
  }
  
  # Since only X/U1 is harmonic SDE, we feed the sigma_x to U2.
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0, 0))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         omega_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X", "Q")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}


## 19. Damped Oscillator SDE and Harmonic Oscillator SDE -----------------
# Here X/U1 is Damped SDE and Y/U3 is Harmonic SDE. This might be different than the user specified x and y.

plot_xy_damp_harmonic_sdes <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_sde"){
    
    # User specified X as Damped SDE and Y as Harmonic SDE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    sigma_x <- input$sigma_x
    
    omega_y <- input$omega_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    sigma_y <- input$sigma_y
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "damp_sde"){
    
    # User specified X as Harmonic SDE and Y as Damped SDE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    sigma_x <- input$sigma_y
    
    omega_y <- input$omega_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    sigma_y <- input$sigma_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
        f <- function(u,p,t) {
          omega_x <- p[1]
          gamma_x <- p[2]
          mu_x <- p[3]
          zeta_x <- p[4]
          
          omega_y <- p[5]
          gamma_y <- p[6]
          mu_y <- p[7]
          
          du1 = u[2]
          du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
          du3 = u[4]
          du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[3] - u[1])
          return(c(du1, du2, du3, du4))
        }
    
  } else if(input$xy_coupling == "coup_mean"){
    
        f <- function(u,p,t) {
          omega_x <- p[1]
          gamma_x <- p[2]
          mu_x <- p[3]
          zeta_x <- p[4]
          
          omega_y <- p[5]
          gamma_y <- p[6]
          mu_y <- p[7]
          
          du1 = u[2]
          du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
          du3 = u[4]
          du4 = omega_y*(u[3] - mu_y) + gamma_y*(u[1] - mu_x)
          return(c(du1, du2, du3, du4))
        }
  }
  
  # Both x and y are stochastic so sigma_x goes to U2 and sigma_y goes to U4. 
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         omega_y, gamma_y, mu_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X", "Q")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}



## 20. Damped Oscillator SDE and ODE ---------------------------------------


plot_xy_damp_sde_ode <- function(input){
  
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  if (input$modeltype_x == "damp_sde"){
    
    # User specified X as Damped SDE and Y as Damped ODE
    omega_x <- input$omega_x
    gamma_x <- input$gamma_x
    start_x <- input$start_x
    mu_x <- input$mean_x
    zeta_x <- input$zeta_x
    sigma_x <- input$sigma_x
    
    omega_y <- input$omega_y
    gamma_y <- input$gamma_y
    start_y <- input$start_y
    mu_y <- input$mean_y
    zeta_y <- input$zeta_y
    
    flipped = FALSE # X is being simulated as X and Y is being simulated as Y
    
  } else if (input$modeltype_y == "damp_sde"){
    
    # User specified X as Damped ODE and Y as Damped SDE but the simulation codes it with flipped notation.
    omega_x <- input$omega_y
    gamma_x <- input$gamma_y
    start_x <- input$start_y
    mu_x <- input$mean_y
    zeta_x <- input$zeta_y
    sigma_x <- input$sigma_y
    
    omega_y <- input$omega_x
    gamma_y <- input$gamma_x
    start_y <- input$start_x
    mu_y <- input$mean_x
    zeta_y <- input$zeta_x
    
    flipped = TRUE # Y is being simulated as X and X is being simulated as Y
  }
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      zeta_y <- p[8]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + zeta_y*u[4] + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3, du4))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      zeta_y <- p[8]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + zeta_y*u[4] + gamma_y*(u[1] - mu_x)
      return(c(du1, du2, du3, du4))
    }
  }
  
  # Since only the x/U1 is stochastic, sigma_x goes to U2.
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0, 0))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         omega_y, gamma_y, mu_y, zeta_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  if(flipped == TRUE){
    colnames(sim_data) = c("times", "Y", "P", "X", "Q")
  }
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
}


## 21. Both are Damped Oscillator SDEs -------------------------------------


plot_xy_damp_sdes <- function(input){
  period_start <- input$period[1]
  period_end <- input$period[2]
  step_size <- input$step_size
  seed_xy <- input$seed
  
  omega_x <- input$omega_x
  gamma_x <- input$gamma_x
  start_x <- input$start_x
  mu_x <- input$mean_x
  zeta_x <- input$zeta_x
  sigma_x <- input$sigma_x
  
  omega_y <- input$omega_y
  gamma_y <- input$gamma_y
  start_y <- input$start_y
  mu_y <- input$mean_y
  zeta_y <- input$zeta_y
  sigma_y <- input$sigma_y
  
  de <- diffeqr::diffeq_setup(JULIA_HOME = path)
  
  if(input$xy_coupling == "coup_diff"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      zeta_y <- p[8]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[1] - u[3])
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + zeta_y*u[4] + gamma_y*(u[3] - u[1])
      return(c(du1, du2, du3, du4))
    }
    
  } else if(input$xy_coupling == "coup_mean"){
    
    f <- function(u,p,t) {
      omega_x <- p[1]
      gamma_x <- p[2]
      mu_x <- p[3]
      zeta_x <- p[4]
      
      omega_y <- p[5]
      gamma_y <- p[6]
      mu_y <- p[7]
      zeta_y <- p[8]
      
      du1 = u[2]
      du2 = omega_x*(u[1] - mu_x) + zeta_x*u[2] + gamma_x*(u[3] - mu_y)
      du3 = u[4]
      du4 = omega_y*(u[3] - mu_y) + zeta_y*u[4] + gamma_y*(u[1] - mu_x)
      return(c(du1, du2, du3, du4))
    }
  }
  
  # Since both are stochastic, U2 will get sigma_x and U4 will get sigma_y
  g <- function(u,p,t) {
    return(c(0, sigma_x, 0, sigma_y))
  }
  
  u0 <- c(start_x, 0, start_y, 0)
  tspan <- c(period_start, period_end)
  p <- c(omega_x, gamma_x, mu_x, zeta_x,
         omega_y, gamma_y, mu_y, zeta_y)
  steps <- seq(period_start, period_end, by = step_size)
  prob <- de$SDEProblem(f,g,u0,tspan, p, seed = seed_xy)
  sol <- de$solve(prob, saveat = steps)
  
  sim_data <- data.frame("t" = sol$t, as.data.frame(t(sapply(sol$u, identity))) )
  colnames(sim_data) <- c("times", "X", "P", "Y", "Q")
  
  plot_data <- as_tibble(sim_data) %>%
    select(-c(P, Q)) %>%
    pivot_longer(cols = c("X", "Y"), values_to = "value", names_to = "variable") %>%
    mutate(time = as.numeric(times),
           value = as.numeric(value))
  
  plot_xy <- bivariate_plot(plot_data)
  
  return(plot_xy)
  
}








