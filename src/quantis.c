/*
 * quantis.c - Normalização por quantis. Veja quantis.h.
 */
#include "quantis.h"

#include <stdlib.h>

#include "util.h"

typedef struct {
    double valor;
    int indice;
} ParValorIndice;

/* Ordena por valor, desempata pelo índice para o resultado ser determinístico. */
static int comparar_pares(const void *a, const void *b)
{
    const ParValorIndice *x = (const ParValorIndice *)a, *y = (const ParValorIndice *)b;
    if (x->valor < y->valor) return -1;
    if (x->valor > y->valor) return 1;
    return (x->indice > y->indice) - (x->indice < y->indice);
}

/* Copia a amostra j para 'pares' e ordena. */
static void ordenar_amostra(const double *m, int n_genes, int j, ParValorIndice *pares)
{
    int g;
    for (g = 0; g < n_genes; g++) {
        pares[g].valor = m[(size_t)j * n_genes + g];
        pares[g].indice = g;
    }
    qsort(pares, (size_t)n_genes, sizeof(ParValorIndice), comparar_pares);
}

void normalizar_quantis(double *m, int n_genes, int n_amostras)
{
    long double *soma_alvo = (long double *)alocar((size_t)n_genes * sizeof(long double));
    double *alvo = (double *)alocar((size_t)n_genes * sizeof(double));
    ParValorIndice *pares = (ParValorIndice *)alocar((size_t)n_genes * sizeof(ParValorIndice));
    int j, g;

    /* Passo 1-2: o alvo de cada posição é a média entre as amostras ordenadas. */
    for (g = 0; g < n_genes; g++) soma_alvo[g] = 0.0L;
    for (j = 0; j < n_amostras; j++) {
        ordenar_amostra(m, n_genes, j, pares);
        for (g = 0; g < n_genes; g++) soma_alvo[g] += pares[g].valor;
    }
    for (g = 0; g < n_genes; g++) alvo[g] = (double)(soma_alvo[g] / (long double)n_amostras);

    /* Passo 3: cada gene recebe o alvo do seu posto (empates: média dos alvos). */
    for (j = 0; j < n_amostras; j++) {
        int inicio = 0;
        ordenar_amostra(m, n_genes, j, pares);
        while (inicio < n_genes) {
            int fim = inicio, k;
            long double soma = 0.0L;
            double valor;
            while (fim + 1 < n_genes && pares[fim + 1].valor == pares[inicio].valor) fim++;
            for (k = inicio; k <= fim; k++) soma += alvo[k];
            valor = (double)(soma / (long double)(fim - inicio + 1));
            for (k = inicio; k <= fim; k++) m[(size_t)j * n_genes + pares[k].indice] = valor;
            inicio = fim + 1;
        }
    }
    free(soma_alvo);
    free(alvo);
    free(pares);
}
