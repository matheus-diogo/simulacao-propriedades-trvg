# Simulação de Monte Carlo

## Salvar o tempo inicial de execução do script
tempo_inicial <- Sys.time()

## Importação dos pacotes necessários
library(tibble)
library(dplyr)
library(magrittr)

## Definição do Método Numérico para Estimação de alpha ####################################################################################

### Determinar a estimativa do parâmetro alpha pelo método de Brent aplicado a função f_alpha
estimar_alpha <- function(amostra) {
  ### Definir uma função em parâmetro alpha, dada uma amostra, para aplicar o processo iterativo
  f_alpha <- function(alpha, amostra) {
    log(alpha) - digamma(alpha) - log(mean(amostra)) + mean(log(amostra))
  }
  ### Aplicar o método de Brent
  intervalo <- c(1/10, 10)
  while (TRUE) {
    estimativa_alpha <- try(
        uniroot(f = f_alpha, interval = intervalo, check.conv = TRUE, amostra = amostra)$root,
        silent = TRUE
    )
    ifelse(
        is.numeric(estimativa_alpha),
        break,
        intervalo <- intervalo * c(1/2, 2)
    )
  }

  return(estimativa_alpha)
}

#### A estimativa do parâmetro beta é expressa como a média amostral dividida pela estimativa do parâmetro alpha

## Determinação da função para o TRVG dada as hipóteses H0: alpha == alpha0 e H1: alpha != alpha0 ##########################################
trvg <- function(amostra, alpha0, significancia) {
  #### Estimar alpha sob espaço paramétrico completo
  est_alpha_ir <- estimar_alpha(amostra = amostra)

  #### Estimar beta sob espaço paramétrico completo
  est_beta_ir <- mean(amostra) / est_alpha_ir

  #### Determinar alpha sob espaço paramétrico restrito
  est_alpha_r <- alpha0

  #### Determinar beta sob espaço paramétrico restrito
  est_beta_r <- mean(amostra) / est_alpha_r

  #### Calcular tamanho amostral e média amostral
  n = length(amostra)
  m = mean(amostra)

  #### Calcular o valor da estatística W
  w = 2*n * log((gamma(est_alpha_r)*(m/est_alpha_r)**est_alpha_r) / (gamma(est_alpha_ir)*(m/est_alpha_ir)**est_alpha_ir)) +
      2*n*est_alpha_r - 2*n*est_alpha_ir + 2*(est_alpha_ir-est_alpha_r) * sum(log(amostra))
  
  #### Calcular o p-valor assintótico
  p_valor <- pchisq(q = w, df = 1, lower.tail = FALSE)

  #### Decidir se deve rejeita H0
  rejeita <- ifelse(p_valor < significancia, TRUE, FALSE)

  #### Calcular o valor crítico da região de rejeição do teste
  c <- qchisq(significancia, df = 1, lower.tail = FALSE)

  return(list(
    estimativa_alpha_irrestrita = est_alpha_ir,
    estimativa_alpha_restrita = est_alpha_r,
    estimativa_beta_irrestrita = est_beta_ir,
    estimativa_beta_restrita = est_beta_r,
    estatistica_w = w,
    p_valor_assintotico = p_valor,
    rejeita_H0 = rejeita,
    valor_critico = c
  ))
}

## Critérios da Simulação ##################################################################################################################

### Definir os tamanhos amostrais em uso
tamanhos_amostrais <- c(20, 30, 50, 100, 200)

### Fixar a quantidade de réplicas para cada tamanho amostral e avaliação
R <- 10000

### Fixar a seed para reprodutibilidade
set.seed(2026)

## Avaliar a distribuição assintótica do teste #############################################################################################

### Criar uma tabela com os valores da estatística W, sob H0, para cada tamanho amostral
amostras_valores_estatistica_w <- tibble(
  tamanho_amostral = integer(), valor_estatistica_w = numeric()
)

### Determinar o nível de significância do teste
significancia <- 0.05

### Calcular os valores da estatística W, sob H0, para cada tamanho amostral
for (n in tamanhos_amostrais) {
  for (r in 1:R) {
    ### Gerar amostra aleatória da distribuição gama com forma == 1 e escala == 2
    amostra <- rgamma(n = n, shape = 1, scale = 2)

    ### Registrar o valor da estatística W dada a amostra gerada
    amostras_valores_estatistica_w <- amostras_valores_estatistica_w %>%
      add_row(
        tamanho_amostral = n,
        valor_estatistica_w = trvg(amostra = amostra, alpha0 = 1, significancia = significancia)$estatistica_w
      )
  }
}

#### Organizar os gráficos para cada tamanho amostral em 2 por 3
layout(matrix(c(
  1, 1, 2, 2, 3, 3,
  0, 4, 4, 5, 5, 0
), 2, 6, byrow = TRUE))

### Exibir os histogramas dos valores da estatística W, sob H0, para cada tamanho amostral
for (n in tamanhos_amostrais) {
  #### Filtrar a tabela de valores da estatística W, sob H0, por tamanho amostral
  tabela_filtrada <- amostras_valores_estatistica_w %>% filter(tamanho_amostral == n)

  #### Gerar o respectivo histograma da tabela filtrada
  hist(
    x = tabela_filtrada$valor_estatistica_w,
    main = paste0('Com ', as.character(n), ' observações'),
    xlab = 'Valor da Estatística W',
    ylab = '',
    freq = FALSE,
  )

  #### Adicionar uma curva da distribuição assintótica teórica da estatística W
  curve(dchisq(x, df = 1), lty = 2, col = 'red', add = TRUE)
}

### Exibir os gráficos Q-Q da estatística W, sob H0, para cada tamanho amostral
for (n in tamanhos_amostrais) {
  #### Filtrar a tabela de valores da estatística W, sob H0, por tamanho amostral
  tabela_filtrada <- amostras_valores_estatistica_w %>% filter(tamanho_amostral == n)

  #### Gerar o respectivo gráfico Q-Q da tabela filtrada
  qqplot(
    x = tabela_filtrada$valor_estatistica_w,
    y = qchisq(ppoints(rchisq(500, df = 1)), df = 1),
    main = paste0('Com ', as.character(n), ' observações'),
    xlab = 'Quantil Empírico',
    ylab = 'Quantil Teórico',
    xlim = c(0, 20),
    ylim = c(0, 20)
  ); abline(0, 1, col = 'red', lty = 2)
}

## Avaliar do tamanho do teste #############################################################################################################

### Criar uma tabela com o tamanho amostral e a decisão de rejeitar H0
amostras_trvg_sob_h0 <- tibble(
    tamanho_amostral = integer(), rejeita = logical()
)

### Determinar o nível de significância exato do teste
significancia_exata <- 0.05

### Calcular o nível de significância estimado do teste dada amostras de uma distribuição Exponencial com escala == 2
for (n in tamanhos_amostrais) {
  for (r in 1:R) {
    #### Gerar a amostra aleatória de determinada distribuição
    amostra <- rexp(n = n, rate = 1/2)

    #### Registrar o resultado do TRVG
    amostras_trvg_sob_h0 <- amostras_trvg_sob_h0 %>% 
      add_row(
        tamanho_amostral = n,
        rejeita = trvg(amostra = amostra, alpha0 = 1, significancia = significancia_exata)$rejeita_H0
    )
  }
}

#### Observe que a distribuição Exponencial é um caso da distribuição Gama com forma == 1

### Construir a tabela de níveis de significância estimado para cada tamanho amostral
amostras_trvg_sob_h0 %>% group_by(tamanho_amostral) %>% summarise(significancia_estimada = mean(rejeita)) -> tabela1

### Salvar o resultado em um .CSV
tabela1 %>% write.csv('tabela1.csv')

#### Restaurar a visualização de um único gráfico
dev.off()

### Exibir a evolução dos níveis de significância estimados com o aumento do tamanho amostral
plot(
  x = tabela1$tamanho_amostral,
  y = tabela1$significancia_estimada,
  main = 'Evolução dos Níveis de Significância Estimados por Tamanho Amostral',
  xlab = 'Tamanho Amostral',
  ylab = 'Nível de Significância Estimado',
  type = 'b'
)
abline(h = significancia_exata, lty = 2, col = 'red')
legend(x = 'topright', legend = 'Nível de Significância Exato', lty = 2, col = 'red')

## Avaliar do poder do teste ###############################################################################################################

### Criar uma tabela com o tamanho amostral e a decisão de rejeitar H0
amostras_trvg_sob_h1 <- tibble(
    tamanho_amostral = integer(), alpha_verdadeiro = numeric(), rejeita = logical()
)

### Determinar o nível de significância do teste
significancia <- 0.05

### Definir diferentes valores de alpha
alphas <- seq(from = 0.5, to = 1.5, by = 0.1)

### Calcular o nível de significância estimado do teste dada amostras de uma distribuição Gama com alpha variados e beta fixo
for (alpha in alphas) {
  for (n in tamanhos_amostrais) {
    for (r in 1:R) {
      #### Gerar a amostra aleatória de determinada distribuição
      amostra <- rgamma(n = n, shape = alpha, scale = 2)
      
      #### Registrar o resultado do TRVG
      amostras_trvg_sob_h1 <- amostras_trvg_sob_h1 %>%
        add_row(
          alpha_verdadeiro = alpha,
          tamanho_amostral = n,
          rejeita = trvg(amostra = amostra, alpha0 = 1, significancia = significancia)$rejeita_H0
        )
      }
  }
}

### Construir a tabela de probabilidades de rejeição estimadas para cada tamanho amostral e valor verdadeiro alpha
amostras_trvg_sob_h1 %>% group_by(tamanho_amostral, alpha_verdadeiro) %>% summarise(prob_rejeicao_estimada = mean(rejeita)) -> tabela2

### Salvar o resultado em um .CSV
tabela2 %>% write.csv('Tabela2.csv')

### Exibir a evolução das estimativas na mudança do valor verdadeiro de alpha com o aumento do tamanho amostral

#### Organizar os gráficos para cada tamanho amostral em 2 por 3
layout(matrix(c(
  1, 1, 2, 2, 3, 3,
  0, 4, 4, 5, 5, 0
), 2, 6, byrow = TRUE))

for (n in tamanhos_amostrais) {
  #### Filtrar a tabela para tamanho amostral n
  tabela2_filtrada <- tabela2 %>% filter(tamanho_amostral == n)
  
  #### Gerar as curvas da função poder dado tamanho amostral n
  plot(
    x = tabela2_filtrada$alpha_verdadeiro,
    y = tabela2_filtrada$prob_rejeicao_estimada,
    main = paste0('Com ', as.character(n), ' observações'),
    xlab = 'Valor Verdadeiro de Alfa',
    ylab = 'Probabilidade de Rejeição Estimada',
    type = 'b',
    col = 'red'
  )
}

## Salvar o tempo total de execução do script
tempo_total <- Sys.time() - tempo_inicial