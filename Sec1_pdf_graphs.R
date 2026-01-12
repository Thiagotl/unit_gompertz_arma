##############################################################
#
# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br)
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
# December/2025
#
##############################################################


dUGo <- function(x, mu=.5, sigma=1.2, tau=.5){
  log(tau)/(1 - mu^(-sigma)) * sigma * x^(-(1 + sigma)) *
    exp((log(tau)/(1 - mu^(-sigma))) * (1 - x^(-sigma)))
}

# Sigma  ---
plot_UGo_vsigma <- function(sigma_vals = c(0.5, 1, 3, 5, 7),
                            mu = 0.5, tau = 0.5,
                            legend_pos = "topright",
                            legend_inset = c(0.02, 0.02),
                            n = 2000, eps = 1e-6,
                            lwd = 2.4,
                            outfile = NULL, w = 9, h = 4.5){
  
  x  <- seq(eps, 1 - eps, length.out = n)
  fx <- sapply(sigma_vals, function(s) dUGo(x, mu = mu, sigma = s, tau = tau))
  
  lty <- c(1, 3, 4, 2, 5)[seq_along(sigma_vals)]
  leg <- lapply(sigma_vals, function(s) bquote(sigma == .(s)))
  
  if(!is.null(outfile)) pdf(outfile, width = w, height = h)
  
  matplot(x, fx, type = "l", col = 1, lty = lty, lwd = lwd,
          xlab = "y", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
  
  axis(1) 
  axis(2, las = 1)  
  
  
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  
  box()
  
  
  mtext("Probability Density Function", side = 2, line = 3)
  
  legend(legend_pos, legend = leg, lty = lty, col = 1, lwd = lwd,
         bty = "n", inset = legend_inset)
  
  if(!is.null(outfile)) dev.off()
  invisible(list(x = x, fx = fx))
}

#  mu  ---
plot_UGo_vmu <- function(mu_vals = c(0.1, 0.3, 0.5, 0.7, 0.9),
                         sigma = 1.2, tau = 0.5,
                         legend_pos = "topright",
                         legend_inset = c(0.02, 0.02),
                         n = 2000, eps = 1e-6,
                         lwd = 2.4,
                         outfile = NULL, w = 9, h = 4.5){
  
  x  <- seq(eps, 1 - eps, length.out = n)
  fx <- sapply(mu_vals, function(m) dUGo(x, mu = m, sigma = sigma, tau = tau))
  
  lty <- c(1, 3, 4, 2, 5)[seq_along(mu_vals)]
  leg <- lapply(mu_vals, function(m) bquote(mu == .(m)))
  
  if(!is.null(outfile)) pdf(outfile, width = w, height = h)
  
  matplot(x, fx, type = "l", col = 1, lty = lty, lwd = lwd,
          xlab = "y", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
  
  axis(1)
  axis(2, las = 1)
  
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  
  box()
  
  mtext("Probability Density Function", side = 2, line = 3)
  
  legend(legend_pos, legend = leg, lty = lty, col = 1, lwd = lwd,
         bty = "n", inset = legend_inset)
  
  if(!is.null(outfile)) dev.off()
  invisible(list(x = x, fx = fx))
}

#plot_UGo_vsigma()
#plot_UGo_vmu()
