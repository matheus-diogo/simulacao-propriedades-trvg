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
replicas <- 10000

## Avaliar do tamanho do teste #############################################################################################################

### Criar uma tabela com o tamanho amostral e a decisão de rejeitar H0
amostras_trvg_sob_h0 <- tibble(
    tamanho_amostral = integer(), rejeita = logical()
)

### Determinar o nível de significância exato do teste
significancia_exata <- 0.05

### Calcular o nível de significância estimado do teste dada amostras de uma distribuição Exponencial com escala == 2
for (n in tamanhos_amostrais) {
  for (r in 1:replicas) {
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

#### Observe que a distribuição Exponencial é um caso da distribuição Gama com alpha = 1

### Construir a tabela de níveis de significância estimado para cada tamanho amostral
amostras_trvg_sob_h0 %>% group_by(tamanho_amostral) %>% summarise(significancia_estimada = mean(rejeita)) -> t_amostrais_sig_estimada

## Avaliar do poder do teste ###############################################################################################################

### Criar uma tabela com o tamanho amostral e a decisão de rejeitar H0
amostras_trvg_sob_h1 <- tibble(
    tamanho_amostral = integer(), rejeita = logical()
)

### Determinar o nível de significância do teste
significancia <- 0.05

### Calcular o nível de significância estimado do teste dada amostras de uma distribuição Gama com alpha variados e beta fixo
for (n in tamanhos_amostrais) {
  for (r in replicas) {
    #### Gerar a amostra aleatória de determinada distribuição
    amostra <- rgamma(n = n, shape = alpha, scale = 2)
    
    #### Construir um gráfico da curva da função poder com o alpha variando de 0 até 2
    #### ...
    }
}

## Salvar o tempo total de execução do script
tempo_total <- Sys.time() - tempo_inicial