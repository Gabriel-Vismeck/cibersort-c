/*
 * estatistica.c - Estatística básica. Veja estatistica.h.
 */
#include "estatistica.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

double media(const double *v, size_t n)
{
    long double soma = 0.0L;
    size_t i;
    for (i = 0; i < n; i++) soma += v[i];
    return (double)(soma / (long double)n);
}

double desvio_padrao(const double *v, size_t n)
{
    long double m = media(v, n), soma = 0.0L;
    size_t i;
    if (n < 2) morrer("o desvio padrão amostral exige pelo menos 2 valores.");
    for (i = 0; i < n; i++) soma += ((long double)v[i] - m) * ((long double)v[i] - m);
    return (double)sqrtl(soma / (long double)(n - 1));
}

int comparar_doubles(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

double mediana(const double *v, size_t n)
{
    double *c = (double *)alocar(n * sizeof(double));
    double resultado;
    memcpy(c, v, n * sizeof(double));
    qsort(c, n, sizeof(double), comparar_doubles);
    resultado = (n % 2 == 1) ? c[n / 2] : (c[n / 2 - 1] + c[n / 2]) / 2.0;
    free(c);
    return resultado;
}

int padronizar(double *v, size_t n)
{
    double m = media(v, n), dp = desvio_padrao(v, n);
    size_t i;
    if (dp == 0.0) return 0;
    for (i = 0; i < n; i++) v[i] = (v[i] - m) / dp;
    return 1;
}

double correlacao(const double *a, const double *b, size_t n)
{
    long double ma = media(a, n), mb = media(b, n), sab = 0.0L, saa = 0.0L, sbb = 0.0L;
    size_t i;
    for (i = 0; i < n; i++) {
        long double da = (long double)a[i] - ma, db = (long double)b[i] - mb;
        sab += da * db;
        saa += da * da;
        sbb += db * db;
    }
    if (saa == 0.0L || sbb == 0.0L) return NAN;
    return (double)(sab / sqrtl(saa * sbb));
}

double rmse(const double *previsto, const double *observado, size_t n)
{
    long double soma = 0.0L;
    size_t i;
    for (i = 0; i < n; i++) {
        long double d = (long double)previsto[i] - (long double)observado[i];
        soma += d * d;
    }
    return (double)sqrtl(soma / (long double)n);
}
