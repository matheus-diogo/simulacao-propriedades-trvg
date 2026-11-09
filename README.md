# Simulação de Monte Carlo para o Estudo das Propriedades do Teste de Razão de Verossimilhaça
Trabalho acadêmico de elaboração da simulação de Monte Carlo para avaliar empiricamente o teste de Razão de Verossimilhaça (TRVG) em hipóteses acerca do parâmetro de forma do Modelo Probabilístico Gama, desenvolvido a partir da linguagem R.

## Objetivos

- Definir as hipóteses sob o parâmetro de forma da distribuição Gama.
- Obter a estatística do TRVG e a distribuição assintótica da estatística W.
- Determinar a região crítica e o p-valor dados os critérios estabelecidos.
- Implementar computacionalmente a função do TRVG que retorna diversas informações do teste.
- Avaliar o tamanho e o poder do teste, além da distribuição assintótica de W, por simulação.

## Tecnologias

|   Pacote  |    Versão |
|:---------:|:---------:|
|    base   |   4.5.0   |
|  datasets |   4.5.0   |
|   dplyr   |   1.2.1   |
|  graphics |   4.5.0   |
| grDevices |   4.5.0   |
|  magrittr |   2.0.5   |
|  methods  |   4.5.0   |
|   stats   |   4.5.0   |
|   tibble  |   3.3.1   |
|   utils   |   4.5.0   |

Ademais, a apresentação foi estruturada pelo LaTeX.

## Resultados

Em síntese, os resultados encontrados na simulação de Monte Carlo foram os esperados do que tem-se no estudo teórico do TRVG. Do que foi observado, é importante ressaltar que: o tamanho empı́rico do teste demonstrou convergência assintótica para o tamanho exato do teste determinado por 5%, mesmo com os ruı́dos de Monte Carlo; o poder do teste aumentou conforme foram considerados diferentes valores verdadeiros do parâmetro de forma alfa que afastaram-se de alfa zero igual à 1, definido pela hipótese nula. Além de que o mesmo ocorre com n aumentado para qualquer valor verdadeiro alternativo de alpha, ou seja, sob hipótese alternativa; a distribuição empı́rica da estatı́stica W apresentou melhores aderências à distribuição Qui-Quadrado com grau de liberdade 1 conforme o aumento do tamanho amostral n, desde a região central da distribuição até a sua cauda à direita, de acordo com o que é visualizável nos histogramas e nos gráficos Q-Q.

Além disso, o script R teve tempo de execução de aproximadamente 33 minutos.

## Licença

Este projeto está sob a licença MIT.
