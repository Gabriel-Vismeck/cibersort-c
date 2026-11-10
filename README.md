# CIBERSORT em C

Reimplementação em C99 do método **CIBERSORT** (Newman et al., 2015, *Nature Methods*) para estimar a
composição celular de amostras de tecido a partir de dados de expressão gênica ("deconvolução").
Sem dependências externas: só a biblioteca padrão do C e a `libm`. Roda em Linux, macOS e Windows.

> **Nota importante.** Este é um projeto independente, escrito a partir da descrição publicada do
> algoritmo. Não é a distribuição oficial e não tem vínculo com os autores do método. As proporções
> estimadas são determinísticas, mas o valor-p depende do gerador de números aleatórios e por isso
> pode diferir levemente do obtido com o script oficial em R. Se for usar em um trabalho a ser
> publicado, valide uma amostra de resultados contra a ferramenta oficial.

---

## Sumário

1. [O que o CIBERSORT faz](#1-o-que-o-cibersort-faz)
2. [Instalação e compilação](#2-instalação-e-compilação)
3. [Preparando os dados](#3-preparando-os-dados)
4. [Executando](#4-executando)
5. [Lendo o resultado](#5-lendo-o-resultado)
6. [Referência de opções](#6-referência-de-opções)
7. [Como o programa funciona por dentro](#7-como-o-programa-funciona-por-dentro)
8. [Solução de problemas](#8-solução-de-problemas)
9. [Estrutura do projeto](#9-estrutura-do-projeto)
10. [Limitações](#10-limitações)
11. [Referências e citação](#11-referências-e-citação)

---

## 1. O que o CIBERSORT faz

Uma amostra de tecido (tumor, sangue, biópsia) é uma **mistura** de vários tipos celulares. Ao medir a
expressão gênica do tecido inteiro (bulk RNA-seq ou microarray), obtém-se a soma das contribuições de
todas essas células. O CIBERSORT faz o caminho inverso: estima **quanto de cada tipo celular** existe
na mistura.

Para isso ele precisa de uma **matriz de assinatura**: uma tabela com a expressão característica de
cada tipo celular de referência. A mais usada é a **LM22**, com 547 genes que distinguem 22 tipos de
células imunes. A ideia central é

```
expressão da mistura  ≈  assinatura  ×  proporções celulares
```

e as proporções são estimadas por regressão com vetores de suporte (ν-SVR), repetida com três
parâmetros ν, ficando com o melhor ajuste. Em seguida, um **teste de permutação** avalia se o ajuste
é melhor do que o esperado ao acaso (o valor-p).

O que você fornece e o que recebe:

| Entrada | Saída |
|---|---|
| Matriz de **expressão** (genes × amostras) | Proporção de cada tipo celular em cada amostra |
| Matriz de **assinatura** (genes × tipos celulares), p. ex. LM22 | Valor-p, correlação e RMSE do ajuste, por amostra |

A matriz LM22 é distribuída pelos autores do método (portal CIBERSORTx, `cibersortx.stanford.edu`).
Obtenha-a lá, respeitando os termos de uso.

---

## 2. Instalação e compilação

Você precisa de um compilador C (gcc ou clang) e do `make`. O programa funciona sem OpenMP; ele só é
necessário se você quiser usar vários núcleos do processador (recomendado quando há muitas amostras
ou permutações).

Primeiro, baixe o projeto:

```sh
git clone https://github.com/Gabriel-Vismeck/cibersort_c.git
cd cibersort
```

### Linux

```sh
# Debian/Ubuntu
sudo apt update && sudo apt install build-essential

# Fedora
sudo dnf install gcc make

make           # compila em build/cibersort
make OPENMP=1  # alternativa: usa todos os núcleos
```

### macOS

```sh
xcode-select --install  # instala clang e make

make                    # compila em build/cibersort
```

O compilador padrão do macOS não inclui OpenMP. Para usar vários núcleos, instale o GCC pelo Homebrew
e aponte o `make` para ele (troque `NN` pela versão instalada, por exemplo `14`):

```sh
brew install gcc
make limpar
make CC=gcc-NN OPENMP=1
```

### Windows


**MSYS2 (gera um `.exe` nativo)**

1. Instale o [MSYS2](https://www.msys2.org/) e abra o terminal **MSYS2 UCRT64**.
2. Instale as ferramentas:
   ```sh
   pacman -S --needed mingw-w64-ucrt-x86_64-gcc make
   ```
3. Dentro da pasta do projeto, compile:
   ```sh
   make  # ou: make OPENMP=1
   ```
   O executável é gerado em `build/cibersort.exe`.

Se preferir não usar o `make`, o comando de compilação direta funciona em qualquer sistema com gcc:

```sh
mkdir -p build
gcc -O2 -std=c99 -Iinclude -o build/cibersort src/*.c -lm          # sem OpenMP
gcc -O2 -std=c99 -fopenmp -Iinclude -o build/cibersort src/*.c -lm # com OpenMP
```

### Conferindo a instalação

```sh
./build/cibersort --ajuda
```

Se a lista de opções aparecer, está tudo certo. Para apagar os arquivos compilados, use `make limpar`.

---

## 3. Preparando os dados

Os dois arquivos de entrada são **tabelas de texto** com genes nas linhas. O separador (tab, `;` ou
`,`) é detectado automaticamente a partir do cabeçalho.

**Assinatura** (genes × tipos celulares):

```
Gene symbol    B cells naive    B cells memory    Plasma cells    ...
ABCB4          555.7            110.2             6.1             ...
ACAP1          9.4              1025.6            12.3            ...
```

**Mistura / expressão** (genes × amostras):

```
Gene symbol    Amostra_01    Amostra_02    Amostra_03
ABCB4          12.5          30.1          8.7
ACAP1          450.2         120.9         600.4
```

Pontos de atenção, que cobrem a maioria dos problemas na prática:

- **Identificadores de gene iguais nas duas tabelas.** O cruzamento é feito pelo nome exato, com
  diferença entre maiúsculas e minúsculas. Use símbolos HGNC (como `ABCB4`) nas duas.
- **Escala linear, não log.** O programa verifica o maior valor da mistura: se for menor que 50, assume
  que os dados estão em log2 e aplica 2^x automaticamente (e avisa no terminal). Valores como TPM,
  FPKM ou contagens normalizadas em escala linear são o formato esperado.
- **Ponto como separador decimal** (`12.5`, não `12,5`).
- **Sem valores ausentes** (`NA`, células vazias). Remova ou substitua esses genes antes.
- **Genes duplicados na mistura:** apenas a primeira ocorrência é usada.
- **Genes duplicados na assinatura** não são permitidos (o programa para com mensagem de erro).
- Genes que existem só em uma das tabelas são ignorados. O programa informa quantos genes ficaram em
  comum; confira se o número é razoável (para a LM22, espera-se a grande maioria dos 547).
- O cabeçalho pode ter a célula do canto (`Gene symbol`) ou não (formato do `write.table()` do R).

---

## 4. Executando

Comando mínimo:

```sh
./build/cibersort --assinatura LM22.txt --expressao expressao.txt --saida resultado.tsv
```

### Receitas comuns

**RNA-seq (TPM, FPKM)** - o padrão, sem normalização por quantis:

```sh
./build/cibersort --assinatura LM22.txt --expressao expressao.txt --saida resultado.tsv --perm 100
```

**Microarray** - ative a normalização por quantis:

```sh
./build/cibersort --assinatura LM22.txt --expressao expressao.txt --saida resultado.tsv --perm 100 --quantis
```

**Valor-p mais estável**, com mais permutações e usando todos os núcleos (binário compilado com OpenMP):

```sh
./build/cibersort --assinatura LM22.txt --expressao expressao.txt --saida resultado.tsv --perm 1000
```

**Limitando o número de núcleos:**

```sh
./build/cibersort --assinatura LM22.txt --expressao expressao.txt --saida resultado.tsv --threads 4
```

**Modo absoluto** (escore de infiltração comparável entre amostras):

```sh
./build/cibersort --assinatura LM22.txt --expressao expressao.txt --saida resultado.tsv --absoluto
```

**Reprodutibilidade:** o sorteio das permutações usa uma semente fixa (padrão `1`), então rodar duas
vezes com os mesmos dados e opções dá exatamente o mesmo arquivo. Use `--semente N` para variar o
sorteio.

### Quantis ou sem:

| Tipo de dado | Recomendação |
|---|---|
| RNA-seq (TPM/FPKM/normalizado) | padrão, **sem** quantis |
| Microarray | `--quantis` |

A normalização por quantis força todas as amostras a terem a mesma distribuição de valores. Isso
corrige diferenças técnicas entre chips, mas em RNA-seq pode apagar diferenças biológicas reais, razão
pela qual a recomendação do método é desligá-la nesse caso. A normalização é aplicada só à mistura,
nunca à assinatura.

---

## 5. Lendo o resultado

O arquivo de saída é um TSV com uma linha por amostra:

```
Mistura   B cells naive  ...  Neutrophils   Valor_p   Correlacao   RMSE
Amostra_01    0.0312     ...      0.0000     0.0100      0.8421   0.5390
```

| Coluna | Significado |
|---|---|
| Um por tipo celular | Proporção estimada. No modo padrão, **a soma por amostra é 1**. |
| `Valor_p` | Confiabilidade da deconvolução (teste de permutação). Quanto **menor**, melhor. |
| `Correlacao` | Correlação de Pearson entre a expressão reconstruída e a observada. Quanto **maior**, melhor. |
| `RMSE` | Erro quadrático médio da reconstrução. Quanto **menor**, melhor. |
| `Escore_absoluto(...)` | Só no modo `--absoluto`: soma das proporções ajustadas, usada como índice geral de infiltração. |

**Como interpretar o valor-p.** Ele responde: "a correlação obtida é maior do que a que se obtém ao
embaralhar os dados?" Na prática, é comum filtrar amostras com `Valor_p` acima de 0,05 antes de
analisar as proporções, pois nelas a deconvolução não é confiável. O valor-p tem resolução limitada
pelo número de permutações: com `--perm 100`, os valores variam em passos de cerca de 0,01. Use mais
permutações (por exemplo 500 a 1000) quando precisar de valores mais finos.

**Dois cálculos de valor-p:**

- `--valor-p original` (padrão): replica o cálculo do script oficial. Pode resultar em exatamente 0.
- `--valor-p classico`: `(1 + nº de permutações com correlação ≥ observada) / (1 + total)`. Nunca dá 0,
  sendo mais conservador.

**Sem permutações** (`--perm 0`) a coluna `Valor_p` recebe o valor `9999`, convenção do original para
"não calculado".

**Modo absoluto.** Em `sig.score` (padrão), as proporções são reescaladas pela mediana de expressão da
amostra em relação à mediana global, e o escore é a soma. Em `no.sumto1`, as proporções não são
divididas pela soma. Escolha com `--metodo-absoluto`.

---

## 6. Referência de opções

| Opção | Descrição | Padrão |
|---|---|---|
| `--assinatura ARQ` | Tabela da assinatura (genes × tipos celulares), escala linear. **Obrigatória.** | — |
| `--expressao ARQ` | Tabela de expressão (genes × amostras). **Obrigatória.** (`--mistura` é aceito como sinônimo.) | — |
| `--saida ARQ` | Arquivo de resultado (`.tsv`). **Obrigatória.** | — |
| `--perm N` | Número de permutações. `0` desliga o valor-p. | `100` |
| `--quantis` | Normaliza a mistura por quantis (microarray). | desligado |
| `--sem-quantis` | Desliga a normalização por quantis (mantida por compatibilidade). | — |
| `--absoluto` | Ativa o modo absoluto. | desligado |
| `--metodo-absoluto M` | `sig.score` ou `no.sumto1`. | `sig.score` |
| `--valor-p T` | `original` ou `classico`. | `original` |
| `--semente N` | Semente do sorteio das permutações. | `1` |
| `--tolerancia X` | Critério de parada do ν-SVR. | `0.001` |
| `--threads N` | Núcleos usados (só em binário compilado com OpenMP). | todos |
| `--ajuda` | Mostra a ajuda. | — |

---

## 7. Como o programa funciona por dentro

```
 Mistura (genes × amostras)           Assinatura (genes × tipos celulares)
          │                                        │
          ▼                                        ▼
 1. genes em comum  ◄──────────────────────────────┤
          │                                        │
          ▼                                        │
 2. pré-processamento                              │
    (desfaz log2; quantis, se pedido)              │
          │                                        ▼
          ├──────────────────────►  3. padronização da assinatura
          │                                        │
          ▼                                        ▼
 4. distribuição nula  ◄──── embaralha a mistura e repete o ajuste
          │
          ▼
 5. deconvolução de cada amostra (ν-SVR, 3 valores de ν, melhor RMSE)
          │
          ▼
 6. proporções + valor-p + correlação + RMSE  →  resultado.tsv
```

**Etapa 5, em detalhe** (módulo `central`). Para cada ν em {0,25; 0,5; 0,75}: ajusta-se a ν-SVR, os
pesos negativos são zerados, os pesos são divididos pela soma (viram proporções) e a amostra é
reconstruída para medir RMSE e correlação. Fica o ν de menor RMSE.

**A ν-SVR vem do aprendizado de máquina, mas aqui funciona como estimador, não
como modelo preditivo: a cada amostra o ajuste é refeito do zero, e os pesos obtidos são o próprio
resultado (as proporções). Nada é aprendido de uma amostra para outra. A assinatura é fixa.

**Sobre as permutações.** Para medir o que o acaso produz, os valores da mistura são embaralhados, o
mesmo algoritmo é rodado e a correlação resultante é guardada. O conjunto dessas correlações forma a
distribuição nula. O valor-p compara a correlação real com ela.

**Paralelismo.** As permutações e as amostras são independentes entre si, então são distribuídas entre
os núcleos (OpenMP). O resultado numérico é idêntico com 1 ou com vários núcleos. Só o tempo muda.

---

## 8. Solução de problemas

| Sintoma | Causa provável e solução |
|---|---|
| `nenhum gene em comum entre a assinatura e a mistura` | Identificadores diferentes (Ensembl × símbolo, maiúsculas/minúsculas). Converta para o mesmo padrão. |
| `valor não numérico '...' no gene '...'` | Vírgula decimal, `NA` ou célula vazia. Use ponto decimal e remova valores ausentes. |
| `linha N: X campos (esperava Y)` | Linha com número de colunas diferente do cabeçalho. Verifique separadores soltos e quebras de linha. |
| `a assinatura tem o gene '...' repetido` | Remova ou colapse as linhas duplicadas na assinatura. |
| Aviso `assumindo log2` | Máximo da mistura < 50. Se seus dados **já** estão em escala linear com valores todos baixos, o programa os interpretará de forma errada; revise a escala dos dados. |
| `ν-SVR não convergiu` | Raro. Aumente a tolerância (`--tolerancia 0.01`) e verifique se há genes com valores extremos. |
| `--threads ignorado` | O binário foi compilado sem OpenMP. Recompile com `make OPENMP=1`. |
| Muitas amostras com `Valor_p` alto | Verifique o tipo de dado (use `--quantis` em microarray), a escala e se a assinatura é compatível com o tecido. |

---

## 9. Estrutura do projeto

```
cibersort/
├── Makefile
├── include/        interfaces públicas (cada .h documenta o seu módulo)
└── src/            implementações
```

| Módulo | Responsabilidade |
|---|---|
| `main` | Ponto de entrada: argumentos → tabelas → pipeline. |
| `cli` | Interpreta a linha de comando e define os valores padrão. |
| `opcoes` | Estrutura com as opções de execução. |
| `pipeline` | Orquestra as etapas (genes em comum, pré-processamento, nula, saída). |
| `central` | Algoritmo central para uma amostra: 3 valores de ν, fica o de menor RMSE. |
| `nusvr` | ν-SVR linear (solver SMO no estilo da LIBSVM). |
| `quantis` | Normalização por quantis (Bolstad et al., 2003). |
| `valor_p` | Valor-p a partir da distribuição nula (`original` ou `classico`). |
| `aleatorio` | Gerador xoshiro256\*\* e sorteio das permutações. |
| `estatistica` | Média, desvio, mediana, correlação, RMSE, padronização. |
| `tabela` | Leitura de tabelas TSV/CSV. |
| `util` | Erro fatal e alocação segura. |

Dependências entre módulos (sem ciclos): `main → cli, pipeline, tabela`;
`pipeline → central, quantis, valor_p, aleatorio, nusvr, estatistica, tabela`;
`central → nusvr, estatistica`; todos usam `util`.

**Por onde começar a ler o código:**

1. `include/pipeline.h` e a função `executar_cibersort()` no fim de `src/pipeline.c`: mostram as 6
   etapas em ordem.
2. `include/central.h` e `include/nusvr.h`: o método em si (a teoria da ν-SVR está documentada em
   `nusvr.h`).
3. Os demais módulos são peças pequenas e independentes.

---

## 10. Limitações

- Leitura de tabelas: não trata separador dentro de aspas.
- O valor-p pode diferir levemente do obtido no script oficial em R (gerador de números aleatórios
  diferente). As proporções, correlação e RMSE não dependem do sorteio.
- Não inclui a matriz LM22 nem outras assinaturas. Elas devem ser obtidas das fontes originais.

---

## 11. Referências e citação

Algoritmos usados na implementação:

- Fan R-E, Chen P-H, Lin C-J. Working set selection using second order information for training
  support vector machines. *JMLR* 6:1889–1918 (2005).
- Chang C-C, Lin C-J. Training ν-support vector regression: theory and algorithms. *Neural
  Computation* 14:1959–1977 (2002).
- Bolstad BM, et al. A comparison of normalization methods for high density oligonucleotide array
  data based on variance and bias. *Bioinformatics* 19:185–193 (2003).
