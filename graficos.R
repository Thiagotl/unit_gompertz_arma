# Gráficos



rm(list = ls())
gc()
library(forecast)
library(ggplot2)
library(ggfortify)
library(readxl)



### TAXA DESEMPREGO ITÁLIA

#dados1 <- read_excel("STP-20250509110327181.xlsx",na = "-")


#dados1 <- read_excel("STP-20250509171131887.xlsx", na = "-")


dados<-na.omit(dados1[4])/100

Y<-ts(dados, frequency = 12)

plot(decompose(Y))
title(main = "") 


Acf(Y)
pacf(Y)


decomp <- decompose(Y)


dec<-autoplot(decomp) +
  theme_bw() +  # fundo branco
  theme(
    panel.grid.major = element_blank(),  # remove grade principal
    panel.grid.minor = element_blank(),  # remove grade secundária
    strip.background = element_blank(),  # remove fundo da faixa dos títulos
    strip.text = element_text(size = 14) # aumenta texto dos títulos
  )

saz<-monthplot(Y, ylab = "")


acf<-ggAcf(Y) +
  ggtitle(NULL) +  # remove título
  theme_bw() +     # fundo branco
  theme(
    panel.grid.major = element_blank(),  # remove grade maior
    panel.grid.minor = element_blank(),  # remove grade menor
    plot.title = element_blank(),        # garante título removido
    panel.border = element_rect(color = "black", fill = NA)
  )

pacf<-ggPacf(Y) +
  ggtitle(NULL) +
  ylab("Partial ACF") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )



# # Garantir que os fatores estejam ordenados
# df_plot$model <- factor(df_plot$model, levels = unique(df_plot$model))
# df_plot$measure <- factor(df_plot$measure, levels = c("MAE", "MAPE", "RMSE"))
# 
# # Gráfico final com fundo branco e sem linhas de grade
# ggplot(df_plot, aes(x = model, y = values, fill = model)) +
#   geom_bar(stat = "identity", width = 0.7, color = "black", position = "dodge") +
#   geom_text(
#     aes(label = round(values, 2)),
#     fontface = "bold",
#     vjust = ifelse(df_plot$values >= 0, -0.3, 1.8),
#     color = "black",
#     size = 2
#   ) +
#   facet_wrap(~measure, scales = "free_y", nrow = 1) +
#   scale_fill_manual(values = gray.colors(length(unique(df_plot$model)))) +
#   coord_cartesian(ylim = c(min(df_plot$values) - 10, max(df_plot$values) + 30)) +
#   labs(y = "Percentage differences", x = "", fill = "") +
#   theme_minimal(base_size = 10) +
#   theme(
#     legend.position = "bottom",
#     axis.text.x = element_text(face = "bold", size = 8),
#     axis.text.y = element_text(face = "bold", size = 8),
#     axis.title.y = element_text(face = "bold", size = 8),
#     strip.text = element_text(face = "bold", size = 9),
#     panel.background = element_rect(fill = "white", colour = "white"),
#     plot.background = element_rect(fill = "white", colour = "white"),
#     panel.grid = element_blank()
#   )





