/*
 * pipeline.c - Pipeline completo do CIBERSORT. Veja pipeline.h.
 */
#include "pipeline.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "aleatorio.h"
#include "central.h"
#include "estatistica.h"
#include "nusvr.h"
#include "quantis.h"
#include "util.h"
#include "valor_p.h"

/* ---------------------------------------------------------------------------
 * Estruturas internas
 * ------------------------------------------------------------------------- */

/* Genes presentes na assinatura E na mistura (na ordem alfabética da assinatura). */
typedef struct {
    int n;                /* quantidade de genes em comum */
    int *pos_mistura;     /* n: linha de cada gene na tabela da mistura */
    int *pos_assinatura;  /* n: linha de cada gene na tabela da assinatura */
} GenesComuns;

/* Problema pronto para a deconvolução (só genes em comum). */
typedef struct {
    int l;           /* genes em comum */
    int c;           /* tipos celulares */
    int n_amostras;
    double *X;       /* l x c, assinatura padronizada (linha por linha) */
    double *K;       /* l x l, núcleo linear de X */
    double *Y;       /* n_amostras x l, uma amostra por vez: Y[amostra * l + gene] */
} Problema;

typedef struct {
    const char *nome;
    int indice;
} NomeIndice;

/* ---------------------------------------------------------------------------
 * Etapa 1: genes em comum
 * ------------------------------------------------------------------------- */

/* Ordena por nome; desempata pelo índice (o primeiro duplicado vem antes). */
static int comparar_nomes(const void *a, const void *b)
{
    const NomeIndice *x = (const NomeIndice *)a, *y = (const NomeIndice *)b;
    int c = strcmp(x->nome, y->nome);
    if (c != 0) return c;
    return (x->indice > y->indice) - (x->indice < y->indice);
}

/* Vetor (nome, índice) ordenado por nome. Liberar com free(). */
static NomeIndice *ordenar_nomes(char **nomes, int n)
{
    NomeIndice *v = (NomeIndice *)alocar((size_t)n * sizeof(NomeIndice));
    int i;
    for (i = 0; i < n; i++) {
        v[i].nome = nomes[i];
        v[i].indice = i;
    }
    qsort(v, (size_t)n, sizeof(NomeIndice), comparar_nomes);
    return v;
}

/*
 * Cruza os nomes de genes das duas tabelas. Genes duplicados na mistura: só o
 * primeiro casa com a assinatura. A assinatura não pode ter genes repetidos.
 */
static GenesComuns encontrar_genes_comuns(const Tabela *assinatura, const Tabela *mistura)
{
    GenesComuns g;
    int n_ass = assinatura->n_linhas, n_mis = mistura->n_linhas, i;
    NomeIndice *ord_ass = ordenar_nomes(assinatura->nomes_linhas, n_ass);
    NomeIndice *ord_mis = ordenar_nomes(mistura->nomes_linhas, n_mis);

    for (i = 1; i < n_ass; i++)
        if (strcmp(ord_ass[i].nome, ord_ass[i - 1].nome) == 0)
            morrer("a assinatura tem o gene '%s' repetido.", ord_ass[i].nome);

    g.n = 0;
    g.pos_mistura = (int *)alocar((size_t)n_ass * sizeof(int));
    g.pos_assinatura = (int *)alocar((size_t)n_ass * sizeof(int));
    for (i = 0; i < n_ass; i++) {
        int lo = 0, hi = n_mis; /* busca binária: primeiro nome igual na mistura */
        while (lo < hi) {
            int meio = (lo + hi) / 2;
            if (strcmp(ord_mis[meio].nome, ord_ass[i].nome) < 0) lo = meio + 1;
            else hi = meio;
        }
        if (lo < n_mis && strcmp(ord_mis[lo].nome, ord_ass[i].nome) == 0) {
            g.pos_mistura[g.n] = ord_mis[lo].indice;
            g.pos_assinatura[g.n] = ord_ass[i].indice;
            g.n++;
        }
    }
    free(ord_ass);
    free(ord_mis);
    if (g.n == 0) morrer("nenhum gene em comum entre a assinatura e a mistura.");
    return g;
}

/* ---------------------------------------------------------------------------
 * Etapa 2: pré-processamento da mistura
 * ------------------------------------------------------------------------- */

/*
 * Reorganiza a mistura para "uma amostra por vez" (colunas[amostra * n_genes + gene]),
 * desfaz o log2 se o máximo for < 50 e aplica a normalização por quantis (opcional).
 * Devolve a matriz (liberar com free()) e a mediana global em *mediana_global.
 */
static double *preparar_mistura(const Tabela *mistura, const Opcoes *op, double *mediana_global)
{
    int n_mis = mistura->n_linhas, n_amostras = mistura->n_colunas, i, j;
    size_t total = (size_t)n_mis * (size_t)n_amostras, k;
    double *colunas = (double *)alocar(total * sizeof(double));
    double maximo = -INFINITY;

    for (i = 0; i < n_mis; i++)
        for (j = 0; j < n_amostras; j++) {
            double v = mistura->dados[(size_t)i * n_amostras + j];
            colunas[(size_t)j * n_mis + i] = v;
            if (v > maximo) maximo = v;
        }
    if (maximo < 50) { /* parece log2: desfaz */
        fprintf(stderr, "Máximo da mistura < 50: assumindo log2 e aplicando 2^x.\n");
        for (k = 0; k < total; k++) colunas[k] = pow(2.0, colunas[k]);
    }
    if (op->normalizar_por_quantis) normalizar_quantis(colunas, n_mis, n_amostras);

    *mediana_global = mediana(colunas, total); /* usada só no modo absoluto */
    if (*mediana_global < 1.0) *mediana_global = 1.0;
    return colunas;
}

/* ---------------------------------------------------------------------------
 * Etapa 3: seleção de genes e padronização da assinatura
 * ------------------------------------------------------------------------- */

/*
 * Monta X (padronizada com média e desvio de toda a matriz), K e Y.
 * 'colunas' é a mistura já pré-processada, com n_mis genes por amostra.
 */
static Problema montar_problema(const Tabela *assinatura, const Tabela *mistura, const GenesComuns *genes, const double *colunas)
{
    Problema p;
    int i, j;
    size_t k, total_x;
    double m, dp;

    p.l = genes->n;
    p.c = assinatura->n_colunas;
    p.n_amostras = mistura->n_colunas;

    total_x = (size_t)p.l * (size_t)p.c;
    p.X = (double *)alocar(total_x * sizeof(double));
    for (i = 0; i < p.l; i++)
        for (j = 0; j < p.c; j++)
            p.X[(size_t)i * p.c + j] = assinatura->dados[(size_t)genes->pos_assinatura[i] * p.c + j];
    m = media(p.X, total_x);
    dp = desvio_padrao(p.X, total_x);
    if (dp == 0.0) morrer("a assinatura não tem variação.");
    for (k = 0; k < total_x; k++) p.X[k] = (p.X[k] - m) / dp;
    p.K = matriz_nucleo_linear(p.X, p.l, p.c);

    /* Y: só os genes em comum, na mesma ordem de X. */
    p.Y = (double *)alocar((size_t)p.l * (size_t)p.n_amostras * sizeof(double));
    for (j = 0; j < p.n_amostras; j++)
        for (i = 0; i < p.l; i++)
            p.Y[(size_t)j * p.l + i] = colunas[(size_t)j * mistura->n_linhas + genes->pos_mistura[i]];
    return p;
}

static void liberar_problema(Problema *p)
{
    free(p->X);
    free(p->K);
    free(p->Y);
}

/* ---------------------------------------------------------------------------
 * Etapa 4: distribuição nula por permutação
 * ------------------------------------------------------------------------- */

/*
 * Embaralha os valores da mistura, roda o algoritmo central e guarda a
 * correlação de cada permutação. Devolve a nula ordenada, sem NaN
 * (liberar com free()); *n_validos recebe o tamanho. Sem permutações: NULL.
 */
static double *calcular_nula(const Problema *pb, const Opcoes *op, int *n_validos)
{
    int n_perm = op->n_perm, p, validos = 0;
    double *nula;
    int *sorteios;

    *n_validos = 0;
    if (n_perm <= 0) return NULL;

    nula = (double *)alocar((size_t)n_perm * sizeof(double));
    sorteios = sortear_indices(n_perm, pb->l, pb->l * pb->n_amostras, op->semente);
    fprintf(stderr, "Rodando %d permutações...\n", n_perm);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (p = 0; p < n_perm; p++) {
        double *yr = (double *)alocar((size_t)pb->l * sizeof(double));
        double *prop = (double *)alocar((size_t)pb->c * sizeof(double));
        int q;
        for (q = 0; q < pb->l; q++) yr[q] = pb->Y[sorteios[(size_t)p * pb->l + q]];
        if (!padronizar(yr, (size_t)pb->l)) {
            nula[p] = NAN;
        } else {
            ResultadoCentral r = algoritmo_central(pb->X, pb->K, pb->l, pb->c, yr, op->metodo, op->tolerancia, prop);
            nula[p] = r.correlacao;
        }
        free(yr);
        free(prop);
    }
    free(sorteios);

    for (p = 0; p < n_perm; p++) /* tira NaN */
        if (!isnan(nula[p])) nula[validos++] = nula[p];
    qsort(nula, (size_t)validos, sizeof(double), comparar_doubles);
    *n_validos = validos;
    return nula;
}

/* ---------------------------------------------------------------------------
 * Etapa 5: deconvolução de cada amostra
 * ------------------------------------------------------------------------- */

/*
 * Devolve tabela n_amostras x n_colunas_saida (liberar com free()).
 * Colunas: c proporções | valor-p | correlação | RMSE | [escore absoluto]
 */
static double *deconvoluir_amostras(const Problema *pb, const Opcoes *op, const double *nula, int n_nula,
                                    double mediana_global, int n_colunas_saida)
{
    double *tabela_saida = (double *)alocar((size_t)pb->n_amostras * (size_t)n_colunas_saida * sizeof(double));
    int j;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (j = 0; j < pb->n_amostras; j++) {
        double *y = (double *)alocar((size_t)pb->l * sizeof(double));
        double *linha = tabela_saida + (size_t)j * n_colunas_saida;
        int q;
        memcpy(y, pb->Y + (size_t)j * pb->l, (size_t)pb->l * sizeof(double));
        if (!padronizar(y, (size_t)pb->l)) {
            for (q = 0; q < n_colunas_saida; q++) linha[q] = NAN;
        } else {
            ResultadoCentral r = algoritmo_central(pb->X, pb->K, pb->l, pb->c, y, op->metodo,
                                                   op->tolerancia, linha);
            if (op->absoluto && op->metodo == METODO_RELATIVO) {
                double med = mediana(pb->Y + (size_t)j * pb->l, (size_t)pb->l);
                for (q = 0; q < pb->c; q++) linha[q] = linha[q] * med / mediana_global;
            }
            linha[pb->c] = valor_p(nula, n_nula, r.correlacao, op->tipo_p); /* 9999 se sem permutações */
            linha[pb->c + 1] = r.correlacao;
            linha[pb->c + 2] = r.rmse;
            if (op->absoluto) {
                double soma = 0.0;
                for (q = 0; q < pb->c; q++) soma += linha[q];
                linha[pb->c + 3] = soma;
            }
        }
        free(y);
    }
    return tabela_saida;
}

/* ---------------------------------------------------------------------------
 * Etapa 6: gravação do resultado
 * ------------------------------------------------------------------------- */

static void gravar_resultado(const char *arquivo_saida, const Tabela *assinatura, const Tabela *mistura,
                             const Opcoes *op, const double *tabela_saida, int n_colunas_saida)
{
    FILE *saida = fopen(arquivo_saida, "w");
    int j, k;
    if (!saida) morrer("falha ao criar o arquivo '%s'.", arquivo_saida);

    fprintf(saida, "Mistura");
    for (k = 0; k < assinatura->n_colunas; k++) fprintf(saida, "\t%s", assinatura->nomes_colunas[k]);
    fprintf(saida, "\tValor_p\tCorrelacao\tRMSE");
    if (op->absoluto)
        fprintf(saida, "\tEscore_absoluto(%s)", op->metodo == METODO_RELATIVO ? "sig.score" : "no.sumto1");
    fprintf(saida, "\n");

    for (j = 0; j < mistura->n_colunas; j++) {
        fprintf(saida, "%s", mistura->nomes_colunas[j]);
        for (k = 0; k < n_colunas_saida; k++)
            fprintf(saida, "\t%.12g", tabela_saida[(size_t)j * n_colunas_saida + k]);
        fprintf(saida, "\n");
    }
    fclose(saida);
    fprintf(stderr, "Resultado gravado em %s\n", arquivo_saida);
}

/* ---------------------------------------------------------------------------
 * Função pública: orquestra todas as etapas
 * ------------------------------------------------------------------------- */

void executar_cibersort(const Tabela *assinatura, const Tabela *mistura, const Opcoes *op,
                        const char *arquivo_saida)
{
    GenesComuns genes;
    Problema problema;
    double *colunas, *nula, *tabela_saida, mediana_global;
    int n_nula, n_colunas_saida = assinatura->n_colunas + 3 + (op->absoluto ? 1 : 0);

    genes = encontrar_genes_comuns(assinatura, mistura);                              /* 1 */
    fprintf(stderr, "Genes na assinatura: %d | genes em comum: %d | amostras: %d\n", assinatura->n_linhas,
            genes.n, mistura->n_colunas);

    colunas = preparar_mistura(mistura, op, &mediana_global);                         /* 2 */
    problema = montar_problema(assinatura, mistura, &genes, colunas);                 /* 3 */
    nula = calcular_nula(&problema, op, &n_nula);                                     /* 4 */
    tabela_saida = deconvoluir_amostras(&problema, op, nula, n_nula, mediana_global,  /* 5 */
                                        n_colunas_saida);
    gravar_resultado(arquivo_saida, assinatura, mistura, op, tabela_saida, n_colunas_saida); /* 6 */

    free(tabela_saida);
    free(nula);
    free(colunas);
    liberar_problema(&problema);
    free(genes.pos_mistura);
    free(genes.pos_assinatura);
}
