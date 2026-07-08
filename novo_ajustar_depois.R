load("rdata/MC_UGoARMA_ARMA11_MC_10k_success_parallel.RData")


library(tidyverse)

tamanhos <- names(MC_out)

# 2. Criar um data.frame único com todos os dados
df_list <- list()

for (n in tamanhos) {
  # Extrair os estimadores para este tamanho de amostra
  estim_temp <- MC_out[[n]]$estim
  
  # Converter para data.frame se for matriz
  if (is.matrix(estim_temp)) {
    estim_temp <- as.data.frame(estim_temp)
  }
  
  # Adicionar coluna identificando o tamanho da amostra
  estim_temp$tamanho <- n
  
  # Adicionar coluna com o ID da replicação (opcional, mas útil)
  estim_temp$replicacao <- 1:nrow(estim_temp)
  
  # Guardar na lista
  df_list[[n]] <- estim_temp
}

# 3. Combinar tudo em um único data.frame
dados_completos <- bind_rows(df_list)

dados_completos <- dados_completos %>%
  rename(
    alpha = V1,
    phi1 = V2,
    theta1 = V3,
    sigma = V4
  )
dados_longos <- dados_completos |> 
  pivot_longer(
    cols = -c(tamanho, replicacao),
    names_to = "parametro",
    values_to = "estimativa"
  )


niveis_ordenados <- sort(unique(as.numeric(dados_longos$tamanho)))

# Converter tamanho para fator com níveis ordenados
dados_longos <- dados_longos %>%
  mutate(
    tamanho = factor(tamanho, 
                     levels = as.character(niveis_ordenados),
                     labels = paste0("n = ", niveis_ordenados))
  )

# Verificar se funcionou
unique(dados_longos$tamanho)



parametro_escolhido <- "sigma"
cor_escolhida <- "#4C78A8"

# dados_longos <- dados_longos %>%
#   mutate(
#     tamanho = factor(
#       tamanho,
#       levels = as.character(sort(as.numeric(as.character(unique(tamanho))))),
#       labels = paste0(
#         "n = ",
#         sort(as.numeric(as.character(unique(tamanho))))
#       )
#     )
#   )




p4 <- ggplot(
  dados_longos %>% filter(parametro == parametro_escolhido),
  aes(x = estimativa)
) +
  geom_histogram(
    aes(y = stat(density)),
    bins = 30,
    fill = cor_escolhida,
    color = "black",
    alpha = 0.65
  ) +
  geom_density(
    color = cor_escolhida,
    linewidth = 1.2
  ) +
  facet_wrap(~ tamanho, scales = "free", ncol = 1) +
  labs(
    x = "Estimated value",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray85", color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )



plot_grid(p1, p2, p3, p4, ncol = 4)
