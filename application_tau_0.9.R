library(forecast)
library(BTSR)
library(tidyverse)
library(e1071)
source("ugo_fit.R")
# source("functions.R")
library(readxl)
library(lubridate)
library(lmtest)

# rm(list = ls())
# gc()

# DATASETS ------

# INTEREST RATE
dados1 <- read_excel("dataset.xlsx", na = "-")

# dados1 <- read_excel("STP-20250509160743592.xlsx", na = "-")

dados <- na.omit(dados1[2]) / 100

# SERIE
Y <- ts(dados, frequency = 12)

m <- length(Y)
h1 <- 6
n <- m - h1

# TRAIN SERIE
y_train <- ts(Y[1:n], frequency = 12)

# TEST SERIE
y_test <- Y[(n + 1):m]

# tend_determ(Y)
# raiz_unit(Y)
# sazonalidade(Y)

#### REGRESSORS MATRIX ----


t <- 1:length(y_train)
t_hat <- (n + 1):(n + h1)

C <- cos(2 * pi * t / 12)
C_hat <- cos(2 * pi * t_hat / 12)

X <- as.matrix(C)
X_hat <- as.matrix(C_hat)

nX <- dim(X)[2]
X0 <- rbind(X, X_hat)

a01 <- auto.arima(y_train, seasonal = FALSE)
new1 <- Arima(y_test, model = a01) # one-step-ahead

a02 <- auto.arima(y_train, xreg = X, seasonal = FALSE)
new2 <- Arima(y_test, xreg = X_hat, model = a02) # one-step-ahead


#### BEST MODEL UGO ARMA----
source("best_ugoarma2.R")

pmax <- 3
qmax <- 3

y <- y_train

# FIT WITH REGRESSORS -----
best_ugoarma <- best_ugo_2(y_train,
  pmax = pmax, qmax = qmax,
  nbest = 8, X = X, X_hat = X_hat
)


if (best_ugoarma$q[1] == 0) {
  fit_ugoarma <- uGoarma.fit(y_train, ar = 1:best_ugoarma$p[1], ma = NA, X = X, X_hat = X_hat)
} else {
  fit_ugoarma <- uGoarma.fit(y_train, ar = 1:best_ugoarma$p[1], ma = 1:best_ugoarma$q[1], X = X, X_hat = X_hat)
}


# FIT WITHOUT REGRESSORS -------

best_ugoarma_sr <- best_ugo_2(y_train,
  pmax = pmax, qmax = qmax,
  nbest = 8
)

if (best_ugoarma_sr$q[1] == 0) {
  fit_ugoarma_sr <- uGoarma.fit(y_train, ar = 1:best_ugoarma_sr$p[1], ma = NA)
} else {
  fit_ugoarma_sr <- uGoarma.fit(y_train, ar = 1:best_ugoarma_sr$p[1], ma = 1:best_ugoarma_sr$q[1])
}

#### BETA E KW APPLICATION ----

quant <- .9 # quantil
# matrix deresiduals # matrix de resultados
order <- matrix(NA, nrow = 16, ncol = 6)
cont <- 1

colnames(order) <- c("p", "q", "barma.AIC", "karma.AIC", "barmax.AIC", "karmax.AIC")

for (i in 0:3) {
  for (j in 0:3) {
    barma <- summary(BARFIMA.fit(y_train,
      p = i, d = FALSE, q = j, info = TRUE,
      start = list(
        phi = rep(0, i),
        theta = rep(0, j)
      ),
      report = FALSE
    ))

    karma1 <- suppressWarnings(KARFIMA.fit(y_train,
      p = i, d = FALSE, q = j, info = TRUE,
      rho = quant,
      control = list(method = "Nelder-Mead", stopcr = 1e-2),
      report = FALSE
    ))
    karma <- summary(karma1)


    barmax <- summary(BARFIMA.fit(y_train,
      p = i, d = FALSE, q = j, info = TRUE,
      start = list(
        phi = rep(0, i),
        theta = rep(0, j),
        beta = coef(lm(y_train ~ X + 0))
      ),
      control = list(method = "Nelder-Mead", stopcr = 1e-2),
      xreg = X,
      report = FALSE
    ))

    karmax1 <- suppressWarnings(KARFIMA.fit(y_train,
      p = i, d = FALSE, q = j, info = TRUE,
      xreg = X, rho = quant,
      control = list(method = "Nelder-Mead", stopcr = 1e-2),
      report = FALSE
    ))
    karmax <- summary(karmax1)

    if (karma1$convergence == 1 || is.nan(karma$aic) == 1) karma$aic <- 0
    if (karmax1$convergence == 1 || is.nan(karmax$aic) == 1) karmax$aic <- 0

    order[cont, ] <- c(
      i, j, barma$aic, karma$aic,
      barmax$aic, karmax$aic
    )
    cont <- cont + 1
  }
}

order <- order[-1, ]
print(order)

orbarma <- order[which(order[, 3] == min(order[, 3])), c(1:3)]
orkarma <- order[which(order[, 4] == min(order[, 4])), c(1:2, 4)]
orbarmax <- order[which(order[, 5] == min(order[, 5])), c(1:2, 5)]
orkarmax <- order[which(order[, 6] == min(order[, 6])), c(1:2, 6)]

barma <- BARFIMA.fit(y_train,
  p = orbarma[1], d = F, q = orbarma[2],
  info = T, report = F
)
karma <- KARFIMA.fit(y_train,
  p = orkarma[1], d = F, q = orkarma[2], rho = quant,
  control = list(method = "Nelder-Mead", stopcr = 1e-2),
  info = T, report = F
)

barmax <- BARFIMA.fit(y_train,
  p = orbarmax[1], d = F, q = orbarmax[2],
  xreg = X, info = T, report = F
)
orbarmax[2] <- 0
orbarmax[1] <- 1

barmax <- BARFIMA.fit(y_train, 
                      p = 1, 
                      d = F, 
                      xreg = X, 
                      info = T, 
                      report = F)

karmax <- KARFIMA.fit(y_train,
  p = orkarmax[1], d = F, q = orkarmax[2], rho = quant,
  control = list(method = "Nelder-Mead", stopcr = 1e-2),
  xreg = X, info = T, report = F
)


results_insample <- rbind(
  forecast::accuracy(fit_ugoarma$fitted, y_train),
  forecast::accuracy(barmax$fitted.values, y_train),
  forecast::accuracy(karmax$fitted.values, y_train),
  forecast::accuracy(a02$fitted, y_train),
  forecast::accuracy(barma$fitted.values, y_train),
  forecast::accuracy(karma$fitted.values, y_train),
  forecast::accuracy(a01$fitted, y_train)
)[, c(3, 2, 5)]

rownames(results_insample) <- c(
  "UGOARMA", "BARMAX", "KARMAX", "a02",
  "BARMA", "KARMA", "a01"
)

round(results_insample, 4)


### out-of-sample

barma_out <- BARFIMA.extract(
  yt = Y,
  coefs = list(
    alpha = barma$coefficients[1],
    phi = barma$coefficients[2:(orbarma[1] + 1)],
    theta = if (orbarma[2] == 0) {
      NULL
    } else {
      barma$coefficients[(orbarma[1] + 2):(orbarma[1] + 1 + orbarma[2])]
    },
    nu = barma$coefficients[(orbarma[1] + 2 + orbarma[2])]
  )
)


karma_out <- KARFIMA.extract(
  yt = Y,
  coefs = list(
    alpha = karma$coefficients[1],
    phi = karma$coefficients[2:(orkarma[1] + 1)],
    theta = if (orkarma[2] == 0) {
      NULL
    } else {
      karma$coefficients[(orkarma[1] + 2):(orkarma[1] + 1 + orkarma[2])]
    },
    nu = karma$coefficients[(orkarma[1] + 2 + orkarma[2])]
  )
)


barmax_out <- BARFIMA.extract(
  yt = Y, xreg = X0,
  coefs = list(
    alpha = barmax$coefficients[1],
    beta = barmax$coefficients[2:(nX + 1)],
    phi = barmax$coefficients[(nX + 2):(orbarmax[1] + nX + 1)],
    theta = if (orbarmax[2] == 0) {
      NULL
    } else {
      barmax$coefficients[(orbarmax[1] + (nX + 2)):(orbarmax[1] + nX + 1 + orbarmax[2])]
    },
    nu = barmax$coefficients[(orbarmax[1] + (nX + 2) + orbarmax[2])]
  )
)


karmax_out <- KARFIMA.extract(
  yt = Y, xreg = X0, rho = quant,
  coefs = list(
    alpha = karmax$coefficients[1],
    beta = karmax$coefficients[2:(nX + 1)],
    phi = karmax$coefficients[(nX + 2):(orkarmax[1] + nX + 1)],
    theta = if (orkarmax[2] == 0) {
      NULL
    } else {
      karmax$coefficients[(orkarmax[1] + nX + 2):(orkarmax[1] + nX + 1 + orkarmax[2])]
    },
    nu = karmax$coefficients[(orkarmax[1] + nX + 2 + orkarmax[2])]
  )
)


# OUT OF SAMPLE FIT WITHOUT REGRESSORS ------

ugoarma_out1 <- KARFIMA.extract(
  yt = Y, rho = quant,
  coefs = list(
    alpha = fit_ugoarma_sr$coeff[1],
    phi = fit_ugoarma_sr$coeff[2:(best_ugoarma_sr$p[1] + 1)],
    theta = if (best_ugoarma_sr$q[1] == 0) {
      NULL
    } else {
      fit_ugoarma_sr$coeff[(best_ugoarma_sr$p[1] + 2):(best_ugoarma_sr$p[1] + 1 + best_ugoarma_sr$q[1])]
    },
    nu = fit_ugoarma_sr$coeff[(best_ugoarma_sr$p[1] + 2 + best_ugoarma_sr$q[1])]
  )
)


# OUT OF SAMPLE FIT WITH REGRESSORS  ------

ugoarma_out2 <- KARFIMA.extract(
  yt = Y, xreg = X0, rho = quant,
  coefs = list(
    alpha = fit_ugoarma$coeff[1],
    beta = fit_ugoarma$coeff[2:(nX + 1)],
    phi = fit_ugoarma$coeff[(nX + 2):(best_ugoarma$p[1] + 1 + nX)],
    theta = if (best_ugoarma$q[1] == 0) {
      NULL
    } else {
      fit_ugoarma$coeff[(best_ugoarma$p[1] + 2 + nX):(best_ugoarma$p[1] + 1 + best_ugoarma$q[1] + nX)]
    },
    nu = fit_ugoarma$coeff[(best_ugoarma$p[1] + 2 + best_ugoarma$q[1] + nX)]
  )
)


a <- n + 1:length(y_test)

results_outsample <- rbind(
  forecast::accuracy(ugoarma_out2$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(barmax_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(karmax_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(new2$fitted, y_test),
  forecast::accuracy(ugoarma_out1$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(barma_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(karma_out$mut[(n + 1)]:length(y_test), y_test),
  forecast::accuracy(new1$fitted, y_test)
)[, c(3, 2, 5)]
#


names_rows <- c(
  "UGOMAX", "BARMAX", "KARMAX", "ARIMAX",
  "UGO", "BARMA", "KARMA", "ARIMA"
)


row.names(results_outsample) <- names_rows

# xtable::xtable((results_outsample),digits=4)

print(round(results_outsample, 4))

round(results_outsample[, ], 4)

# df_percent_diff<-sweep(sweep(results_outsample, 2, results_outsample[1, ], FUN = "-"),
#                        2, results_outsample[1, ], FUN = "/")*100
#
#
# df_percent_diff<-df_percent_diff[1:4,]
#
# #
# # df<-data.frame(values=as.vector(df_percent_diff),
# #                model=rep(names_rows,3),
# #                measure=rep(c("MAE","MAPE","RMSE"),4)
# # )
#
#
# df <- data.frame(
#   values = as.vector(df_percent_diff),  # MAE1, MAE2, ..., RMSE1, RMSE2, ...
#   model = rep(rownames(df_percent_diff), times = ncol(df_percent_diff)),
#   measure = rep(colnames(df_percent_diff), each = nrow(df_percent_diff))
# )
#
# df_percent_diff
# # Remover o ARIMAX
# #df <- df %>% filter(model != "ARIMAX")
#
# ggplot(df,
#        aes(y = values, x = measure,
#            fill = factor(model))
# ) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.9) +
#   labs(fill = "", y = "Percentage differences", x = "") +
#   scale_fill_manual(values = c("#222222","#666666","#aaaaaa" ,"#ffffff")) +
#   ylim(0, 850) +
#   geom_text(aes(label = ifelse(values >= 0,
#                                sprintf("%.1f", values),
#                                sprintf("%.2f", values))),
#             position = position_dodge(width = 0.9),
#             fontface = "bold",
#             vjust = ifelse(df$values >= 0, -0.3, 1.2),  # Ajuste fino para posicionamento
#             hjust = 0.5,
#             angle = 0,
#             color = "black",
#             size = 2.5) +  # Aumentei um pouco o tamanho
#   theme(legend.position = "bottom",
#         strip.text = element_text(face = "bold", size = 8),
#         plot.title = element_text(face = "bold", size = 8),
#         legend.text = element_text(face = "bold", size = 8),
#         legend.key.size = unit(0.3, "cm"),
#         legend.spacing.x = unit(0.3, 'cm'),
#         legend.margin = margin(t = -13, unit = "pt"),
#         axis.title.y = element_text(face = "bold", color = "black", size = 8),
#         axis.title.x = element_text(face = "bold", color = "black", size = 8),
#         axis.text.x = element_text(face = "bold", color = "black", size = 8),
#         axis.text.y = element_text(face = "bold", color = "black", size = 8),
#         panel.background = element_rect(fill = "white", colour = "white"))

# which(fit_ugoarma$residuals < -3)

# RESIDUALS
# acf<-ggAcf(fit_ugoarma$residuals) +
#   ggtitle(NULL) +  # remove título
#   theme_bw() +     # fundo branco
#   theme(
#     panel.grid.major = element_blank(),  # remove grade maior
#     panel.grid.minor = element_blank(),  # remove grade menor
#     plot.title = element_blank(),        # garante título removido
#     panel.border = element_rect(color = "black", fill = NA),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14)
#   )
#
# pacf<-ggPacf(fit_ugoarma$residuals) +
#   ggtitle(NULL) +
#   ylab("Partial ACF") +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.title = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14)
#   )


# round(fit_ugoarma$model,4)
#
#
# round(karmax_out$coefs,4)
#
# coeftest()

# yt <- as.vector(t(ugoarma_out2$yt))
# mut <- ugoarma_out2$mut
#
#
# burn_in <- 6  # ajuste conforme o modelo
# mut_aligned <- c(rep(NA, burn_in), mut[(burn_in + 1):length(mut)])
#
# time <- seq(as.Date("2011-03-01"), by = "month", length.out = length(yt))
#
#
# df_obs <- data.frame(
#   Time = time,
#   Value = yt,
#   Type = "Observed data"
# )
#
# df_pred <- data.frame(
#   Time = time,
#   Value = mut_aligned,
#   Type = "Predicted Median"
# )
#
# df_plot <- bind_rows(df_obs, df_pred)
#
# ggplot(df_plot, aes(x = Time, y = Value, color = Type, linetype = Type)) +
#   geom_line(size = 0.6) +
#   scale_color_manual(values = c("Observed data" = "black", "Predicted Median" = "#D2042D")) + #D2042D , "#0047AB
#   scale_linetype_manual(values = c("Observed data" = "solid", "Predicted Median" = "dashed")) +
#   labs(x = "Time", y = "Rate of credit operations") +
#   scale_x_date(date_breaks = "2 year", date_labels = "%Y")+
#   theme_minimal(base_size = 12) +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 14, color = "black")
#         )
