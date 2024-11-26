library(ggplot2)
library(tidyr)
library(dplyr)

# PARA AS MÉDIAS --------------------------------------------------------------
data <- data.frame(
  n = c(70, 150, 300, 500, 1000),
  alpha = c(0.9993, 0.9975, 1.0009, 0.9987, 0.9998),
  phi = c(0.2034, 0.2025, 0.2004, 0.2005, 0.2005),
  theta = c(0.3945, 0.3964, 0.3993, 0.3990, 0.3995),
  sigma = c(6.3520, 6.1723, 6.0792, 6.0518, 6.0243)
)


data_long<- data |> 
  pivot_longer(cols = -n, names_to = "parameter", values_to = "value")


ggplot(data_long, aes(x=n, y=value, color = parameter, group = parameter)) +
  geom_line(size = 1) +
  geom_point(size = 3)+
  labs(
    title = "Média Estimada dos Parâmetros",
    x = "Tamanho Amostral (n)",
    y = "Média Estimada",
    color = "Parâmetro"
  ) +
  scale_x_continuous(breaks = c(70, 150, 300, 500, 1000)) +  
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

  
# PARA O VIÉS RELATIVO --------------------------------------------------------

rb_data <- data.frame(
  n = c(70, 150, 300, 500, 1000),
  alpha = c(0.0704, 0.2483, -0.0941, 0.1280, 0.0197),
  phi = c(-1.6905, -1.2445, -0.1955, -0.2648, -0.2553),
  theta = c(1.3703, 0.8981, 0.1632, 0.2509, 0.1333),
  sigma = c(-5.8671, -2.8712, -1.3193, -0.8636, -0.4042)
)


rb_data_long <- rb_data %>%
  pivot_longer(cols = -n, names_to = "parameter", values_to = "value")

ggplot(rb_data_long, aes(x = n, y = value, color = parameter, group = parameter)) +
  geom_line(size = 1) +  
  geom_point(size = 3) +  
  labs(
    title = "Viés Relativo (RB%) dos Parâmetros",
    x = "Tamanho Amostral (n)",
    y = "Viés Relativo (%)",
    color = "Parâmetro"
  ) +
  scale_x_continuous(breaks = c(70, 150, 300, 500, 1000)) +  
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# PARA ERRO QUADRÁTICO MÉDIO --------------------------------------------------
mse_data <- data.frame(
  n = c(70, 150, 300, 500, 1000),
  alpha = c(0.0672, 0.0278, 0.0133, 0.0078, 0.0040),
  phi = c(0.0275, 0.0111, 0.0052, 0.0031, 0.0016),
  theta = c(0.0258, 0.0098, 0.0044, 0.0025, 0.0013),
  sigma = c(1.2299, 0.4941, 0.2282, 0.1294, 0.0642)
)

mse_data_long <- mse_data %>%
  pivot_longer(cols = -n, names_to = "parameter", values_to = "value")


ggplot(mse_data_long, aes(x = as.factor(n), y = value, fill = parameter)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(
    title = "Erro Quadrático Médio (MSE) dos Parâmetros",
    x = "Tamanho Amostral (n)",
    y = "MSE",
    fill = "Parâmetro"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# PARA A COBERTURA DO INTERVALO DE CONFIANÇA ----------------------------------

cr_data <- data.frame(
  n = c(70, 150, 300, 500, 1000),
  alpha = c(0.9225, 0.9399, 0.9471, 0.9484, 0.9456),
  phi = c(0.9204, 0.9427, 0.9448, 0.9483, 0.9469),
  theta = c(0.9052, 0.9335, 0.9422, 0.9496, 0.9489),
  sigma = c(0.9347, 0.9411, 0.9468, 0.9501, 0.9488)
)


cr_data_long <- cr_data %>%
  pivot_longer(cols = -n, names_to = "parameter", values_to = "value")


ggplot(cr_data_long, aes(x = n, y = value, color = parameter, group = parameter)) +
  geom_line(size = 1) +  
  geom_point(size = 3) +  
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +  
  labs(
    title = "Cobertura do Intervalo de Confiança (CR%) dos Parâmetros",
    x = "Tamanho Amostral (n)",
    y = "CR (%)",
    color = "Parâmetro"
  ) +
  scale_x_continuous(breaks = c(70, 150, 300, 500, 1000)) +  
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.1))) +  
    
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

 # PARA Faceted Plots  --------------------------------------------------------

faceted_data <- data.frame(
  n = rep(c(70, 150, 300, 500, 1000), times = 3),  # Tamanhos amostrais replicados para cada modelo
  model = rep(c("UGo-ARMA(1,1)", "UGo-ARMA(1,0)", "UGo-ARMA(0,1)"), each = 5),  # Modelos
  alpha = c(0.9993, 0.9975, 1.0009, 0.9987, 0.9998, 1.0126, 1.0054, 1.0023, 1.0016, 1.0007, 1.0065, 1.0028, 1.0010, 1.0009, 1.0001),
  phi = c(0.2034, 0.2025, 0.2004, 0.2005, 0.2005, 0.1920, 0.1973, 0.1985, 0.1987, 0.1994, NA, NA, NA, NA, NA),  # Apenas onde aplicável
  theta = c(0.3945, 0.3964, 0.3993, 0.3990, 0.3995, NA, NA, NA, NA, NA, 0.3994, 0.3993, 0.4001, 0.4001, 0.4003),
  sigma = c(6.3520, 6.1723, 6.0792, 6.0518, 6.0243, 6.2997, 6.1372, 6.0695, 6.0451, 6.0200, 6.2341, 6.0978, 6.0555, 6.0277, 6.0178)
)


faceted_data_long <- faceted_data %>%
  pivot_longer(cols = c(alpha, phi, theta, sigma), names_to = "parameter", values_to = "value")




ggplot(faceted_data_long, aes(x = n, y = value, color = parameter, group = parameter)) +
  geom_line(size = 1) +  
  geom_point(size = 3) +  
  labs(
    title = "Métricas por Modelo (Faceted Plots)",
    x = "Tamanho Amostral (n)",
    y = "Valor",
    color = "Parâmetro"
  ) +
  scale_x_continuous(breaks = c(70, 150, 300, 500, 1000)) +  
  facet_wrap(~model, ncol = 1) +  
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Heatmap para Análise Comparativa  -------------------------------------------

heatmap_data <- data.frame(
  n = rep(c(70, 150, 300, 500, 1000), times = 3),
  model = rep(c("UGo-ARMA(1,1)", "UGo-ARMA(1,0)", "UGo-ARMA(0,1)"), each = 5),
  alpha = c(0.0704, 0.2483, -0.0941, 0.1280, 0.0197, -1.2553, -0.5437, -0.2255, -0.1572, -0.0727, -0.6533, -0.2841, -0.1033, -0.0863, -0.0104),
  phi = c(-1.6905, -1.2445, -0.1955, -0.2648, -0.2553, 3.9799, 1.3503, 0.7620, 0.6277, 0.3155, NA, NA, NA, NA, NA),
  theta = c(1.3703, 0.8981, 0.1632, 0.2509, 0.1333, NA, NA, NA, NA, NA, 0.1501, 0.1652, -0.0332, -0.0356, -0.0654),
  sigma = c(-5.8671, -2.8712, -1.3193, -0.8636, -0.4042, -4.9946, -2.2869, -1.1577, -0.7511, -0.3333, -3.9024, -1.6308, -0.9243, -0.4614, -0.2971)
)


heatmap_data_long <- heatmap_data %>%
  pivot_longer(cols = c(alpha, phi, theta, sigma), names_to = "parameter", values_to = "value")



ggplot(heatmap_data_long, aes(x = factor(n), y = parameter, fill = value)) +
  geom_tile(color = "white") +  
  facet_wrap(~model, ncol = 1) +  
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "RB (%)"
  ) +
  labs(
    title = "Heatmap: RB (%) por Modelo, Tamanho Amostral e Parâmetro",
    x = "Tamanho Amostral (n)",
    y = "Parâmetro"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    strip.text = element_text(size = 12),  
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "right"
  )

# Média Estimada e o Viés Relativo (RB%) --------------------------------------

combined_data <- data.frame(
  n = rep(c(70, 150, 300, 500, 1000), times = 2),
  metric = rep(c("mean", "rb"), each = 5),
  alpha = c(0.9993, 0.9975, 1.0009, 0.9987, 0.9998, 
            0.0704, 0.2483, -0.0941, 0.1280, 0.0197),
  phi = c(0.2034, 0.2025, 0.2004, 0.2005, 0.2005, 
          -1.6905, -1.2445, -0.1955, -0.2648, -0.2553),
  theta = c(0.3945, 0.3964, 0.3993, 0.3990, 0.3995, 
            1.3703, 0.8981, 0.1632, 0.2509, 0.1333),
  sigma = c(6.3520, 6.1723, 6.0792, 6.0518, 6.0243, 
            -5.8671, -2.8712, -1.3193, -0.8636, -0.4042)
)


combined_data_long <- combined_data %>%
  pivot_longer(cols = -c(n, metric), names_to = "parameter", values_to = "value")





# Gráfico combinado de Média e RB (%)
ggplot() +
  # Média Estimada
  geom_line(data = filter(combined_data_long, metric == "mean"),
            aes(x = n, y = value, color = parameter, group = parameter), size = 1) +
  geom_point(data = filter(combined_data_long, metric == "mean"),
             aes(x = n, y = value, color = parameter, group = parameter), size = 3) +
  # Viés Relativo (RB%)
  geom_line(data = filter(combined_data_long, metric == "rb"),
            aes(x = n, y = value * 10, color = parameter, group = parameter), 
            linetype = "dashed", size = 1) +
  geom_point(data = filter(combined_data_long, metric == "rb"),
             aes(x = n, y = value * 10, color = parameter, group = parameter), size = 3) +
  # Configuração dos eixos
  scale_y_continuous(
    name = "Média Estimada",
    sec.axis = sec_axis(~ . / 10, name = "Viés Relativo (RB%)")
  ) +
  labs(
    title = "Gráfico Combinado: Média Estimada e Viés Relativo",
    x = "Tamanho Amostral (n)",
    color = "Parâmetro"
  ) +
  scale_x_continuous(breaks = c(70, 150, 300, 500, 1000)) +  # Ajuste do eixo X
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )











