# Pacotes necessários
library(ggplot2)

# Função para gerar a densidade
dUGo <- function(x, mu = 0.5, sigma, tau = 0.5) {
  fx1 <- log(tau) / (1 - mu^(-sigma)) * sigma * x^(-(1 + sigma)) * 
    exp((log(tau) / (1 - mu^(-sigma))) * (1 - x^(-sigma)))
  return(fx1)
}

# Função para gerar o gráfico em preto e branco com ajustes estéticos
plot_density_bw <- function(mu = 0.5, tau = 0.5, sigma_values = c(0.5, 1, 3, 5, 7)) {
  # Criando sequência de valores de x
  x <- seq(0.01, 1, length.out = 1000) # Evitar x = 0 para evitar erros
  
  # Gerando densidades para diferentes sigma
  data <- data.frame()
  for (sigma in sigma_values) {
    y <- dUGo(x, mu = mu, sigma = sigma, tau = tau)
    data <- rbind(data, data.frame(x = x, y = y, sigma = paste("σ =", sigma)))
  }
  
  # Gerando o gráfico
  ggplot(data, aes(x = x, y = y, linetype = sigma)) +
    geom_line(size = 1.2, color = "black") +
    labs(
      title = "Probability Density Function",
      x = "y",
      y = "Density",
      linetype = "σ"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      panel.grid = element_blank(), # Remove as linhas do fundo
      panel.background = element_blank(), # Fundo branco
      axis.line = element_line(color = "black"), # Linhas dos eixos em preto
      legend.position = c(0.8, 0.8), # Coloca a legenda dentro do gráfico (80% x e y)
      legend.background = element_rect(fill = "white", color = "black"), # Fundo da legenda branco com borda preta
      legend.key.width = unit(2, "line") # Aumenta o tamanho dos ícones da legenda
    )
}

# Exemplo de uso
plot_density_bw()
