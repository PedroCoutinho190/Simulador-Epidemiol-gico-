# Simulador Epidemiológico SIR/SEIR

Aplicação web interativa para simulação da propagação de doenças infectocontagiosas, desenvolvida como atividade acadêmica da disciplina de Computação.

---

## 🔗 Demo ao vivo

"https://github.com/PedroCoutinho190/Simulador-Epidemiol-gico-.git"

---

## 📁 Arquivos do repositório

| Arquivo | Descrição |
|---|---|
| `index.html` | Simulador interativo completo (HTML + CSS + JS em um único arquivo) |
| `artigo.md` | Artigo científico explicando o modelo matemático e a implementação |
| `README.md` | Este arquivo |

---

## 🤖 Ferramentas de IA utilizadas

| Entrega | IA utilizada | Modelo |
|---|---|---|
| Código (`index.html`) | Claude (Anthropic) | Claude Sonnet 4.6 |
| Artigo (`artigo.md`) | Claude (Anthropic) | Claude Sonnet 4.6 |

---

## 💬 Prompt utilizado para gerar o código

```
You are an expert epidemiological modeler and software developer. Create a
complete, self-contained, interactive web application (single HTML file with
embedded CSS and JavaScript) that simulates the spread of an infectious disease
using the SIR (Susceptible–Infected–Recovered) compartmental model, extended
to SEIR when appropriate.

The application must include:

1. MATHEMATICAL MODEL
   - Implement the classic SIR/SEIR differential equations
   - Parameters: population size (N), basic reproduction number (R0),
     infection rate (beta), recovery rate (gamma), incubation period (if SEIR),
     initial infected count
   - Numerical integration using Euler or Runge-Kutta 4th order method
   - Time horizon: up to 365 days

2. USER INTERFACE
   - Sliders or input fields for all parameters, with real-time updates
   - A time-series chart showing S(t), I(t), R(t) curves simultaneously
   - Key metrics displayed: peak infection day, peak infected count,
     total recovered (epidemic size), herd immunity threshold
   - A "Run simulation" button and a "Reset to defaults" button

3. VISUALIZATION
   - Use Chart.js (loaded from CDN) for the graph
   - Color code each curve clearly (Susceptible = blue, Infected = red,
     Recovered = green)
   - Responsive design that works on desktop and mobile
   - Display current parameter values next to each slider

4. CODE QUALITY
   - Add inline comments explaining each section of the code
   - Organize code into clearly named functions (initChart, runSimulation,
     updateMetrics, etc.)
   - Use descriptive variable names matching epidemiological terminology

5. PRESETS
   - Include at least 3 preset scenarios: Measles (high R0 ~15),
     COVID-19 (R0 ~2.5), Seasonal Flu (R0 ~1.3)
   - Each preset button loads the corresponding parameters automatically

Make the application visually polished, scientifically accurate, and
immediately usable by someone with no prior coding experience.
Output only the complete HTML file contents, ready to save and open in a browser.
Faça em português.
```

---

## 💬 Prompt utilizado para gerar o artigo

```
Você é um especialista em epidemiologia computacional e redação acadêmica.
Escreva um artigo científico em português explicando o código de simulação
epidêmica apresentado abaixo. O artigo deve seguir a estrutura acadêmica
padrão e ser adequado para submissão a um periódico de epidemiologia
computacional ou informática em saúde pública.

Estrutura:
- Título (informativo e acadêmico)
- Resumo (150–200 palavras): contexto, objetivo, métodos, resultados, conclusão
- 1. Introdução: contexto da modelagem epidêmica, histórico dos modelos SIR
  (Kermack & McKendrick, 1927), motivação
- 2. Modelo Matemático: descrição formal das equações diferenciais, definição
  dos parâmetros com unidades e faixas típicas, hipóteses e limitações
- 3. Implementação: descrição do método numérico utilizado (Runge-Kutta RK4),
  arquitetura do software, principais funções e seus papéis
- 4. Resultados: descrição do que a simulação produz, interpretação das curvas,
  explicação das métricas (R₀, limiar de imunidade de rebanho, pico de infecção)
- 5. Discussão: valor educacional da ferramenta, comparação com surtos reais
  (COVID-19, sarampo), limitações e melhorias futuras
- 6. Conclusão
- Referências: pelo menos 8 referências acadêmicas reais e citáveis no formato
  ABNT, incluindo o artigo original de Kermack-McKendrick

Diretrizes de estilo:
- Português formal acadêmico, terceira pessoa, voz passiva quando apropriado
- Sem listas com marcadores dentro do corpo do artigo — use prosa contínua
- Defina todos os símbolos matemáticos na primeira ocorrência
- Extensão mínima: 1.500 palavras

Produza o artigo em formato Markdown, pronto para colar em um editor de texto.
```

---

## 📐 Sobre o modelo matemático

O simulador implementa os modelos compartimentais **SIR** e **SEIR**, clássicos da epidemiologia matemática, propostos originalmente por Kermack e McKendrick (1927). A integração numérica é realizada pelo método de **Runge-Kutta de 4ª ordem (RK4)**, com passo de 0,5 dia.

**Parâmetros principais:**

| Símbolo | Nome | Descrição |
|---|---|---|
| N | População | Tamanho total da população simulada |
| R₀ | Número reprodutivo básico | Média de casos gerados por um infectado |
| β | Taxa de transmissão | Frequência de contatos infecciosos (dia⁻¹) |
| γ | Taxa de recuperação | Inverso da duração da infecciosidade (dia⁻¹) |
| σ | Taxa de progressão | Inverso do período de incubação — apenas SEIR (dia⁻¹) |

**Cenários disponíveis:** Sarampo (R₀ ≈ 15), COVID-19 (R₀ ≈ 2,5), Gripe sazonal (R₀ ≈ 1,3), Ebola (R₀ ≈ 1,8).

---

## 🚀 Como executar localmente

1. Baixe o arquivo `index.html`
2. Abra-o diretamente no navegador (Chrome, Firefox, Edge)
3. Não é necessária nenhuma instalação adicional

---

*Atividade acadêmica — fins exclusivamente educacionais.*
