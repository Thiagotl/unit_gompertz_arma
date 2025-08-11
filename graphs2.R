
graph_forcast <- function(df){
  
  h <- length(df)#y_test
  
  df_plot <- data.frame(
    h = 1:h,
    Actual = as.numeric(df),
    KARMA = tail(karmax_out$mut, h),
    UGOARMA = tail(ugoarma_out2$mut, h)
  )
  
  
  df_long <- df_plot |>
    pivot_longer(cols = -h, names_to = "Model", values_to = "Value") |>
    mutate(Model = recode(Model,
                          "Actual" = "Actual Values",
                          "UGOARMA" = "UGo-ARMA",
                          "KARMA" = "KARMA"))
  
}

# Horizonte de previsão
h <- length(y_test)

# Cria dataframe com os valores reais e predições dos modelos UGOARMA e KARMA
df_plot <- data.frame(
  h = 1:h,
  Actual = as.numeric(y_test),
  KARMA = tail(karmax_out$mut, h),
  UGOARMA = tail(ugoarma_out2$mut, h)
)

# Transforma para formato longo
df_long <- df_plot |>
  pivot_longer(cols = -h, names_to = "Model", values_to = "Value") |>
  mutate(Model = recode(Model,
                        "Actual" = "Actual Values",
                        "UGOARMA" = "UGo-ARMA",
                        "KARMA" = "KARMA"))
ggplot(df_long, aes(x = h, y = Value, linetype = Model, shape = Model, color = Model)) +
  geom_line(data = subset(df_long, Model != "Actual Values"), size = 0.9) +
  geom_point(data = subset(df_long, Model == "Actual Values"), size = 2) +
  labs(
    x = "h",
    y = "Rate of credit operations"
  ) +
  scale_color_manual(values = c(
    "Actual Values" = "black",
    "UGo-ARMA" = "#CC0000",
    "KARMA" = "#0018A8"
  )) +
  scale_linetype_manual(values = c(
    "UGo-ARMA" = "longdash",
    "KARMA" = "twodash",
    "Actual Values" = "blank"
  )) +
  scale_shape_manual(values = c(
    "Actual Values" = 1,
    "KARMA" = NA,
    "UGo-ARMA" = NA
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 17, color = "black"),
    axis.title = element_text(size = 17, color = "black"),
    legend.title = element_blank(),
    legend.position = c(0.15, 0.15), # <- canto inferior esquerdo
    legend.background = element_rect(fill = "white", color = NA), # fundo branco sem borda
    legend.key = element_blank()
  )


# Organiza os dados
yt <- as.vector(t(ugoarma_out2$yt))
mut <- ugoarma_out2$mut
mut_karmax <- karmax_out$mut

# Alinha os dados (se necessário)
burn_in <- 6  # ajuste conforme o modelo
mut_aligned <- c(rep(NA, burn_in), mut[(burn_in + 1):length(mut)])
mut_karmax_aligned <- c(rep(NA, burn_in), mut_karmax[(burn_in + 1):length(mut_karmax)])

# Sequência temporal
time <- seq(as.Date("2011-03-01"), by = "month", length.out = length(yt))

# Cria os dataframes
df_obs <- data.frame(
  Time = time,
  Value = yt,
  Type = "Observed data"
)

df_pred_ugoarma <- data.frame(
  Time = time,
  Value = mut_aligned,
  Type = "UGo-ARMA"
)

df_pred_karmax <- data.frame(
  Time = time,
  Value = mut_karmax_aligned,
  Type = "KARMA"
)

# Junta tudo
df_plot <- bind_rows(df_obs, df_pred_ugoarma, df_pred_karmax)

# Gera o gráfico
ggplot(df_plot, aes(x = Time, y = Value, color = Type, linetype = Type)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = c(
    "Observed data" = "black",
    "UGo-ARMA" = "#CC0000",
    "KARMA" = "#0018A8"
  )) +
  scale_linetype_manual(values = c(
    "Observed data" = "solid",
    "UGo-ARMA" = "longdash",  # Tracejado longo
    "KARMA" = "twodash"        # Tracejado duplo
  )) +
  labs(x = "Time", y = "Rate of credit operations") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.15), # <- canto inferior direito
    legend.background = element_rect(fill = "white", color = NA), # fundo branco sem borda
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_text(size = 17, color = "black")
  )

