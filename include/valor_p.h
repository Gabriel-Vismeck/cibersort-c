/*
 * valor_p.h - Valor-p a partir da distribuição nula (permutações).
 *
 * P_ORIGINAL: o cálculo do CIBERSORT, 1 - (posição do valor da nula mais próximo de r) / (tamanho da nula). Pode dar 0.
 * P_CLASSICO: (1 + nº de valores nulos >= r) / (1 + tamanho da nula).
 */
#ifndef CIBERSORT_VALOR_P_H
#define CIBERSORT_VALOR_P_H

#define P_ORIGINAL 0
#define P_CLASSICO 1

/*
 * 'nula' deve estar ordenada (crescente) e ter 'n' valores.
 * Sem permutações (n == 0) devolve 9999, como o CIBERSORT original.
 */
double valor_p(const double *nula, int n, double r, int tipo);

#endif /* CIBERSORT_VALOR_P_H */
