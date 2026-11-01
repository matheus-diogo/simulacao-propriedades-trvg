# Simulação de Monte Carlo

## Salvar o tempo inicial de execução do script
tempo_inicial <- Sys.time()

## Importação dos pacotes necessários
library(tibble)
library(dplyr)
library(magrittr)

## Definição do Método Numérico para Estimação de alpha e beta #############################################################################

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

  #### Calcular o valor da estatística W
  w = 2*length(amostra)*
    log((gamma(est_alpha_r)*(mean(amostra)/est_alpha_r)**est_alpha_r)/(gamma(est_alpha_ir)*(mean(amostra)/est_alpha_ir)**est_alpha_ir))+
    2*length(amostra)*est_alpha_r-2*length(amostra)*est_alpha_ir+2*(est_alpha_ir-est_alpha_r)*sum(log(amostra))
  
  #### Calcular o p-valor assintótico
  p_valor <- pchisq(q = w, df = 1, lower.tail = FALSE)

  #### Decidir se deve rejeita H0
  rejeita <- ifelse(p_valor < significancia, TRUE, FALSE)

  return(list(
    estimativa_alpha_irrestritra = est_alpha_ir,
    estimativa_alpha_restritra = est_alpha_r,
    estimativa_beta_irrestritra = est_beta_ir,
    estimativa_beta_restritra = est_beta_r,
    estatistica_w = w,
    p_valor_assintotico = p_valor,
    rejeita_H0 = rejeita
  ))
}

## Avaliar do tamanho do teste #############################################################################################################
## ...

## Avaliar do poder do teste ###############################################################################################################
## ...

## Salvar o tempo total de execução do script
tempo_total <- Sys.time() - tempo_inicial