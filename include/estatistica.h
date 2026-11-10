/*
 * estatistica.h - Estatística básica.
 *
 * Convenções iguais às do R: desvio padrão amostral (n-1) e correlação de
 * Pearson. As somas usam long double para perder o mínimo de precisão.
 */
#ifndef CIBERSORT_ESTATISTICA_H
#define CIBERSORT_ESTATISTICA_H

#include <stddef.h>

double media(const double *v, size_t n);

/* Desvio padrão amostral. Encerra com erro se n < 2. */
double desvio_padrao(const double *v, size_t n);

/* Comparador para qsort() em ordem crescente de double. */
int comparar_doubles(const void *a, const void *b);

/* Mediana (trabalha numa cópia; o vetor original não é alterado). */
double mediana(const double *v, size_t n);

/* Padroniza no lugar: (v - média) / desvio. Devolve 0 se o desvio for zero, 1 caso contrário. */
int padronizar(double *v, size_t n);

/* Correlação de Pearson entre a e b. Devolve NAN se algum dos vetores for constante. */
double correlacao(const double *a, const double *b, size_t n);

/* Raiz do erro quadrático médio entre previsto e observado. */
double rmse(const double *previsto, const double *observado, size_t n);

#endif /* CIBERSORT_ESTATISTICA_H */
