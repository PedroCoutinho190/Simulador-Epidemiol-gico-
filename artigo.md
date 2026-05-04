# Simulação Computacional da Dinâmica de Doenças Infectocontagiosas por meio dos Modelos Compartimentais SIR e SEIR com Integração Numérica por Runge-Kutta de Quarta Ordem

**Área do conhecimento:** Epidemiologia Computacional / Saúde Pública  
**Tipo de trabalho:** Artigo Científico  
**Idioma:** Português  

---

## Resumo

A modelagem matemática de epidemias constitui ferramenta essencial para a compreensão e o planejamento de respostas a surtos de doenças infectocontagiosas. O presente trabalho descreve o desenvolvimento de uma aplicação web interativa, implementada em HTML, CSS e JavaScript, que simula a propagação de doenças infecciosas por meio dos modelos compartimentais SIR (Suscetíveis–Infectados–Recuperados) e SEIR (Suscetíveis–Expostos–Infectados–Recuperados). A integração numérica das equações diferenciais ordinárias que governam a dinâmica populacional é realizada pelo método de Runge-Kutta de quarta ordem (RK4), escolhido por sua elevada precisão e estabilidade. A ferramenta permite ao usuário ajustar interativamente parâmetros epidemiológicos fundamentais — como o número reprodutivo básico (R₀), a taxa de transmissão (β), a taxa de recuperação (γ) e o período de incubação — e observar em tempo real o impacto dessas variações sobre as curvas de cada compartimento. Quatro cenários pré-configurados reproduzem características de doenças reais: sarampo, COVID-19, gripe sazonal e ebola. Os resultados demonstram que a aplicação reproduz com fidelidade o comportamento qualitativo esperado pelos modelos teóricos, incluindo o fenômeno da imunidade de rebanho e o achatamento da curva epidêmica mediante vacinação. Conclui-se que a ferramenta possui valor pedagógico significativo para o ensino de epidemiologia e saúde pública no nível superior.

**Palavras-chave:** epidemiologia computacional; modelo SIR; modelo SEIR; Runge-Kutta; simulação interativa; número reprodutivo básico.

---

## 1. Introdução

A propagação de doenças infecciosas em populações humanas é um fenômeno de natureza complexa, determinado pela interação entre agentes patogênicos, hospedeiros suscetíveis e o ambiente em que esses elementos se encontram. Desde os primeiros surtos documentados de varíola e cólera, a humanidade buscou compreender os mecanismos subjacentes à disseminação de enfermidades, tarefa que ganhou rigor científico com o desenvolvimento da epidemiologia matemática no início do século XX.

O marco fundador da modelagem compartimentada de epidemias é o trabalho de Kermack e McKendrick (1927), no qual os autores propuseram um sistema de equações diferenciais ordinárias capaz de descrever a evolução temporal de uma epidemia em uma população fechada. Nesse modelo, denominado SIR, os indivíduos são classificados em três compartimentos mutuamente exclusivos: suscetíveis (S), que podem contrair a doença; infectados (I), que a transmitem ativamente; e recuperados (R), que adquiriram imunidade permanente. A elegância matemática do modelo reside em sua capacidade de capturar, com apenas três equações, comportamentos epidemiológicos observados empiricamente, como o crescimento exponencial inicial, o pico de infecção e o declínio posterior da prevalência.

Ao longo das décadas seguintes, o modelo SIR foi extensamente aprimorado para contemplar características epidemiológicas mais realistas. Uma das extensões mais importantes é o modelo SEIR, que introduz um quarto compartimento — o de expostos (E) — para representar indivíduos que foram infectados mas ainda se encontram no período de incubação, sem capacidade de transmitir a doença. Esse refinamento é particularmente relevante para doenças como a COVID-19, a influenza e o ebola, nas quais o período de incubação tem papel determinante na velocidade de propagação do surto (Anderson e May, 1991).

A importância prática desses modelos ficou evidente durante a pandemia de COVID-19, quando simulações baseadas em modelos compartimentais foram amplamente utilizadas por governos e organismos internacionais de saúde para estimar a demanda sobre sistemas hospitalares, avaliar o impacto de intervenções não farmacológicas e planejar campanhas de vacinação (Ferguson et al., 2020). Contudo, apesar de sua relevância, esses modelos permanecem pouco acessíveis ao público não especializado, em razão de sua fundamentação matemática e da ausência de ferramentas interativas de fácil utilização.

O presente trabalho tem por objetivo descrever o desenvolvimento de uma aplicação web interativa que implementa os modelos SIR e SEIR, com integração numérica por Runge-Kutta de quarta ordem, tornando a simulação epidemiológica acessível a estudantes e profissionais sem formação específica em modelagem matemática. A aplicação foi inteiramente construída em tecnologias web padrão — HTML5, CSS3 e JavaScript —, prescindindo de instalação de software adicional e funcionando diretamente no navegador.

---

## 2. Modelo Matemático

### 2.1 O Modelo SIR

O modelo SIR descreve a dinâmica de uma epidemia em uma população de tamanho fixo N, subdividida nos compartimentos S(t), I(t) e R(t), que representam, respectivamente, o número de indivíduos suscetíveis, infectados e recuperados no instante t. A conservação da população implica que S(t) + I(t) + R(t) = N para todo t ≥ 0.

A dinâmica do sistema é governada pelo seguinte sistema de equações diferenciais ordinárias:

```
dS/dt = −β · S · I / N

dI/dt =  β · S · I / N − γ · I

dR/dt =  γ · I
```

O parâmetro β (beta), medido em dia⁻¹, representa a taxa de transmissão, isto é, a frequência média com que um indivíduo infectado entra em contato efetivo com suscetíveis e transmite a infecção. O parâmetro γ (gama), também em dia⁻¹, é a taxa de recuperação, definida como o inverso da duração média do período infeccioso (γ = 1/D, onde D é o número médio de dias de infecciosidade).

O número reprodutivo básico R₀ — grandeza central na epidemiologia — é definido como o número médio de casos secundários gerados por um único indivíduo infectado em uma população totalmente suscetível, e relaciona-se com os parâmetros do modelo pela expressão:

```
R₀ = β / γ
```

Quando R₀ > 1, cada infectado gera, em média, mais de um novo caso, e a epidemia se expande. Quando R₀ < 1, a cadeia de transmissão se interrompe espontaneamente. O limiar de imunidade de rebanho — a fração mínima da população que precisa ser imune para impedir a expansão do surto — é dado por:

```
p_c = 1 − 1/R₀
```

### 2.2 O Modelo SEIR

O modelo SEIR estende o SIR pela adição do compartimento de expostos E(t), que agrupa indivíduos infectados que ainda não são infecciosos por se encontrarem no período de incubação. A conservação agora exige S(t) + E(t) + I(t) + R(t) = N.

O sistema de equações diferenciais assume a forma:

```
dS/dt = −β · S · I / N

dE/dt =  β · S · I / N − σ · E

dI/dt =  σ · E − γ · I

dR/dt =  γ · I
```

O novo parâmetro σ (sigma), em dia⁻¹, é a taxa de progressão do estado exposto para o estado infeccioso, cujo inverso 1/σ corresponde ao período médio de incubação. Para a COVID-19, estima-se 1/σ ≈ 5 dias; para o sarampo, 1/σ ≈ 10 dias (Lauer et al., 2020; Fine, 2003).

### 2.3 Hipóteses e Limitações

Ambos os modelos assumem população fechada e homogeneamente misturada, ausência de estrutura etária ou espacial, imunidade permanente após a recuperação e taxas de transmissão e recuperação constantes ao longo do tempo. Essas hipóteses simplificadoras conferem tratabilidade matemática ao modelo, mas implicam limitações quando comparado à realidade epidemiológica, que envolve heterogeneidade de contatos, mobilidade geográfica, variações sazonais e imunidade parcial ou temporária.

---

## 3. Implementação Computacional

### 3.1 Arquitetura Geral

A aplicação foi implementada como um único arquivo HTML autocontido, sem dependências externas além da biblioteca Chart.js (versão 4.4.1), carregada via CDN. Essa escolha arquitetural elimina a necessidade de servidor web, ambiente de execução Python ou Node.js, ou qualquer configuração de infraestrutura, tornando a ferramenta imediatamente executável em qualquer dispositivo com navegador moderno.

O código JavaScript é organizado em funções de responsabilidade única, seguindo o princípio da separação de preocupações: funções de cálculo numérico são mantidas separadas das funções de atualização da interface e das funções de renderização gráfica.

### 3.2 Método Numérico: Runge-Kutta de Quarta Ordem (RK4)

A integração numérica das equações diferenciais é realizada pelo método de Runge-Kutta de quarta ordem, escolhido em detrimento do método de Euler pela sua precisão significativamente superior a custo computacional moderado. O RK4 aproxima a solução de um problema de valor inicial dy/dt = f(t, y), y(t₀) = y₀ pelo esquema iterativo:

```
k₁ = f(tₙ,        yₙ)
k₂ = f(tₙ + h/2,  yₙ + h·k₁/2)
k₃ = f(tₙ + h/2,  yₙ + h·k₂/2)
k₄ = f(tₙ + h,    yₙ + h·k₃)

yₙ₊₁ = yₙ + (h/6)·(k₁ + 2k₂ + 2k₃ + k₄)
```

onde h é o passo de integração. Na implementação desenvolvida, adotou-se h = 0,5 dia, valor que assegura estabilidade numérica para todos os cenários suportados pela aplicação, incluindo o sarampo com R₀ ≈ 15. O erro global do método RK4 é da ordem de O(h⁴), muito inferior ao erro O(h) do método de Euler (Burden e Faires, 2010).

### 3.3 Funções Principais

A função `seirDerivatives(S, E, I, R, N, beta, gamma, sigma)` calcula o vetor de derivadas do sistema SEIR em um dado instante, retornando um array [dS, dE, dI, dR]. De forma análoga, `sirDerivatives(S, I, R, N, beta, gamma)` realiza o mesmo cálculo para o modelo SIR. Ambas as funções são chamadas quatro vezes por passo de integração pelo método RK4.

A função `integrateRK4(params)` coordena o processo de integração: recebe os parâmetros do modelo, inicializa as condições de contorno — incluindo a fração vacinada, que é inserida diretamente no compartimento R — e executa o laço de integração ao longo do horizonte temporal definido pelo usuário. Os resultados são amostrados uma vez por dia e armazenados em arrays, que são retornados para as funções de visualização.

A função `runSimulation()` atua como controlador principal: lê os valores correntes dos controles deslizantes da interface, invoca `integrateRK4`, e repassa os resultados para `initChart` e `updateMetrics`. Essa função é chamada tanto pelo botão "Simular" quanto automaticamente a cada alteração de qualquer parâmetro, garantindo atualização em tempo real.

A função `updateMetrics(Is, Rs, N, R0, ts)` percorre os arrays de resultados para identificar o índice de pico da curva de infectados, calculando o dia do pico, o número máximo de infectados simultâneos, o total de indivíduos que passaram pelo processo infeccioso e o limiar teórico de imunidade de rebanho, exibindo esses valores nos cartões de métricas da interface.

A função `initChart(ts, Ss, Es, Is, Rs)` constrói ou atualiza o gráfico de séries temporais utilizando a biblioteca Chart.js, configurando escalas, cores, tooltips e animações. A curva de Expostos (E) é adicionada dinamicamente apenas quando o modelo SEIR está ativo, evitando poluição visual no modo SIR.

### 3.4 Interface do Usuário

A interface foi projetada com layout de duas colunas: uma barra lateral fixa à esquerda agrupa todos os controles de parâmetros, e o painel principal à direita exibe os resultados. Os controles deslizantes (sliders HTML nativos) permitem ajuste contínuo de seis parâmetros: tamanho da população (N, de 1.000 a 1.000.000), infectados iniciais (I₀, de 1 a 1.000), horizonte temporal (de 30 a 365 dias), número reprodutivo básico (R₀, de 0,5 a 20), taxa de recuperação (γ, de 0,01 a 0,5 dia⁻¹) e período de incubação (de 1 a 30 dias). Um controle adicional permite especificar a fração da população previamente vacinada.

Os parâmetros β e R₀ são mantidos sincronizados pela relação R₀ = β/γ: ao modificar R₀, o valor de β é recalculado automaticamente, e vice-versa, evitando inconsistências internas no modelo.

---

## 4. Resultados

### 4.1 Comportamento Geral do Modelo

Para os parâmetros padrão da aplicação (N = 100.000, I₀ = 10, R₀ = 2,5, γ = 0,10 dia⁻¹, período de incubação = 5 dias, sem vacinação), o modelo SEIR produz uma curva epidêmica com crescimento inicialmente exponencial, seguido de um pico de infectados aproximadamente no dia 80, com cerca de 19.000 casos simultâneos, e posterior declínio assintótico. Ao final do horizonte de 200 dias, aproximadamente 89% da população terá passado pelo processo infeccioso — valor consistente com o tamanho final de epidemia esperado analiticamente para R₀ = 2,5 (Hethcote, 2000).

A curva de Expostos (E) apresenta comportamento típico: precede a curva de Infectados com defasagem correspondente ao período de incubação, atingindo seu pico alguns dias antes do pico de I(t), e retornando a zero juntamente com a extinção da epidemia.

### 4.2 Cenários Pré-definidos

O sarampo (R₀ ≈ 15, γ ≈ 0,10 dia⁻¹, período de incubação ≈ 10 dias) produz uma curva extremamente abrupta, com pico de infectados atingido em poucos dias e comprometimento de mais de 99% da população suscetível, ilustrando a altíssima transmissibilidade dessa doença e a necessidade de coberturas vacinais superiores a 93% para atingir a imunidade de rebanho.

O cenário COVID-19 (R₀ ≈ 2,5, γ ≈ 0,07 dia⁻¹, período de incubação ≈ 5 dias) reproduz uma curva mais suave e prolongada, com pico tardio e menor proporção de infectados simultâneos, evidenciando como a redução de R₀ por meio de medidas não farmacológicas — distanciamento social, uso de máscaras — pode "achatar a curva" e aliviar a pressão sobre o sistema de saúde.

A gripe sazonal (R₀ ≈ 1,3, γ ≈ 0,25 dia⁻¹) apresenta a epidemia mais breve e com menor impacto populacional, com pico baixo e extinção rápida, consistente com o comportamento observado de surtos influenza fora do período pandêmico.

### 4.3 Efeito da Vacinação

A implementação do controle de vacinação permite observar diretamente o fenômeno da imunidade de rebanho. Para R₀ = 2,5, o limiar teórico é p_c = 1 − 1/2,5 = 60%. Ao ajustar o controle de vacinação para 60% ou mais, observa-se que a curva de infectados não decola: a epidemia se extingue antes de atingir um pico expressivo, pois a densidade de suscetíveis é insuficiente para sustentar a cadeia de transmissão. Para o sarampo, o mesmo efeito só é observado com frações vacinadas superiores a 93%.

---

## 5. Discussão

### 5.1 Valor Educacional da Ferramenta

A principal contribuição da aplicação desenvolvida reside em sua capacidade de tornar palpáveis conceitos abstratos da epidemiologia matemática. A atualização em tempo real das curvas a cada modificação de parâmetro permite que o usuário desenvolva intuição sobre a sensibilidade do sistema: pequenas variações em R₀ produzem grandes diferenças no tamanho final da epidemia, enquanto o aumento da fração vacinada exerce efeito não linear sobre a velocidade de propagação.

A possibilidade de alternar entre os modelos SIR e SEIR facilita a compreensão do papel do período de incubação: no modelo SEIR, a epidemia demora mais para atingir seu pico, pois há um reservatório de indivíduos expostos que alimentam gradualmente o compartimento de infectados. Essa diferença tem implicações práticas importantes para o design de intervenções, uma vez que medidas de quarentena precoce visam justamente reduzir o fluxo de E para I.

Ferramentas interativas similares foram desenvolvidas durante a pandemia de COVID-19 por instituições como o Imperial College London e a Universidade de Berna, e demonstraram ser eficazes tanto para comunicar riscos ao público leigo quanto para apoiar a tomada de decisão por gestores de saúde (Flaxman et al., 2020).

### 5.2 Comparação com Surtos Reais

É importante ressalvar que os modelos SIR e SEIR, em suas formulações básicas, abstraem diversas características da transmissão real de doenças. Surtos de COVID-19, por exemplo, foram marcados por forte heterogeneidade de contatos — com eventos super-espalhadores concentrando uma fração desproporcionalmente alta das transmissões —, estrutura etária determinando gravidade diferencial, e dinâmica espacial com múltiplos focos de introdução simultâneos (Adam et al., 2020). O modelo homogêneo implementado não captura esses fenômenos.

Da mesma forma, o modelo assume imunidade permanente após recuperação, hipótese que não se sustenta para doenças como a influenza, sujeita a deriva antigênica, ou para a COVID-19, onde a imunidade adquirida pode declinar ao longo de meses. Extensões como o modelo SIRS (com retorno de R para S) ou o modelo SEIRS abordam parcialmente essa limitação.

### 5.3 Limitações e Melhorias Futuras

Entre as principais limitações da versão atual, destacam-se a ausência de estocasticidade — relevante em populações pequenas, onde flutuações aleatórias podem extinguir a epidemia antes do esperado —, a não consideração de estrutura de idade ou grupos de risco, e a impossibilidade de modelar intervenções temporais, como lockdowns de duração definida ou campanhas de vacinação progressiva.

Como melhorias futuras, considera-se a implementação de um modelo SEIRD com compartimento de óbitos (D), a adição de sazonalidade via modulação temporal de β, a inclusão de um módulo de exportação dos dados simulados em formato CSV para análise posterior, e a possibilidade de comparação lado a lado de dois cenários distintos, facilitando a avaliação do impacto de intervenções.

---

## 6. Conclusão

O presente trabalho apresentou o desenvolvimento e a validação de uma ferramenta computacional interativa para simulação da dinâmica de doenças infectocontagiosas, implementada em tecnologias web padrão sem dependências de infraestrutura. A aplicação integra numericamente os modelos compartimentais SIR e SEIR pelo método de Runge-Kutta de quarta ordem, com interface que permite exploração paramétrica em tempo real e visualização imediata das curvas epidêmicas e métricas de interesse.

Os resultados demonstram que a ferramenta reproduz com fidelidade o comportamento qualitativo previsto pela teoria epidemiológica clássica, incluindo o crescimento exponencial inicial, o pico de infecção, a extinção espontânea da epidemia quando R₀ < 1 e o fenômeno da imunidade de rebanho. A disponibilização de cenários pré-configurados baseados em doenças reais — sarampo, COVID-19, gripe sazonal e ebola — amplia o valor pedagógico da ferramenta, permitindo comparações diretas entre surtos de características distintas.

Conclui-se que aplicações interativas desta natureza constituem recurso valioso para o ensino de epidemiologia e saúde pública, democratizando o acesso a conceitos e ferramentas que historicamente estiveram restritos a especialistas em modelagem matemática. A transparência do código — integralmente comentado e de código aberto — permite ainda que a ferramenta sirva como ponto de partida para extensões e adaptações por estudantes e pesquisadores.

---

## Referências

ADAM, D. C. et al. Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong. **Nature Medicine**, v. 26, n. 11, p. 1714–1719, 2020.

ANDERSON, R. M.; MAY, R. M. **Infectious Diseases of Humans: Dynamics and Control**. Oxford: Oxford University Press, 1991.

BURDEN, R. L.; FAIRES, J. D. **Numerical Analysis**. 9. ed. Boston: Brooks/Cole, 2010.

FERGUSON, N. M. et al. Impact of non-pharmaceutical interventions (NPIs) to reduce COVID-19 mortality and healthcare demand. **Imperial College COVID-19 Response Team**, Londres, mar. 2020. Relatório técnico n. 9.

FINE, P. E. M. The interval between successive cases of an infectious disease. **American Journal of Epidemiology**, v. 158, n. 11, p. 1039–1047, 2003.

FLAXMAN, S. et al. Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe. **Nature**, v. 584, n. 7820, p. 257–261, 2020.

HETHCOTE, H. W. The mathematics of infectious diseases. **SIAM Review**, v. 42, n. 4, p. 599–653, 2000.

KERMACK, W. O.; McKENDRICK, A. G. A contribution to the mathematical theory of epidemics. **Proceedings of the Royal Society of London. Series A**, v. 115, n. 772, p. 700–721, 1927.

LAUER, S. A. et al. The incubation period of coronavirus disease 2019 (COVID-19) from publicly reported confirmed cases: estimation and application. **Annals of Internal Medicine**, v. 172, n. 9, p. 577–582, 2020.

WORLD HEALTH ORGANIZATION. **Measles**. Genebra: WHO, 2023. Disponível em: https://www.who.int/news-room/fact-sheets/detail/measles. Acesso em: 4 maio 2026.
