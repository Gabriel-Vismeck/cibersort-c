/*
 * central.c - Algoritmo central do CIBERSORT (3 valores de ν, melhor RMSE).
 */
#include "central.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "estatistica.h"
#include "nusvr.h"
#include "util.h"

ResultadoCentral algoritmo_central(const double *X, const double *K, int l, int c, const double *y, int metodo, double tolerancia, double *prop)
{
    static const double NUS[3] = {0.25, 0.5, 0.75};
    double *coef = (double *)alocar((size_t)l * sizeof(double));
    double *w = (double *)alocar((size_t)c * sizeof(double));
    double *melhor_w = (double *)alocar((size_t)c * sizeof(double));
    double *reconstruido = (double *)alocar((size_t)l * sizeof(double));
    ResultadoCentral melhor;
    int tem_melhor = 0, t, i, k;

    melhor.rmse = NAN;
    melhor.correlacao = NAN;
    melhor.nu = NAN;

    for (t = 0; t < 3; t++) {
        InfoNuSVR info;
        double soma = 0.0, erro, r;
        treinar_nusvr(K, l, y, NUS[t], 1.0, tolerancia, 100000000L, coef, &info);
        if (!info.convergiu) fprintf(stderr, "Aviso: ν-SVR não convergiu (ν=%g).\n", NUS[t]);

        /* Pesos w = Σ_i coef_i x_i; negativos viram 0. */
        for (k = 0; k < c; k++) w[k] = 0.0;
        for (i = 0; i < l; i++)
            if (coef[i] != 0.0)
                for (k = 0; k < c; k++) w[k] += coef[i] * X[(size_t)i * c + k];
        for (k = 0; k < c; k++) {
            if (w[k] < 0.0) w[k] = 0.0;
            soma += w[k];
        }
        if (soma == 0.0) continue; /* solução inútil */

        /* Reconstrução com as proporções normalizadas e qualidade do ajuste. */
        for (i = 0; i < l; i++) {
            double s = 0.0;
            for (k = 0; k < c; k++) s += X[(size_t)i * c + k] * (w[k] / soma);
            reconstruido[i] = s;
        }
        erro = rmse(reconstruido, y, (size_t)l);
        r = correlacao(reconstruido, y, (size_t)l);
        if (!tem_melhor || erro < melhor.rmse) { /* which.min: fica com o primeiro menor */
            tem_melhor = 1;
            melhor.rmse = erro;
            melhor.correlacao = r;
            melhor.nu = NUS[t];
            memcpy(melhor_w, w, (size_t)c * sizeof(double));
        }
    }

    if (!tem_melhor) {
        for (k = 0; k < c; k++) prop[k] = NAN;
    } else {
        double soma = 0.0;
        for (k = 0; k < c; k++) soma += melhor_w[k];
        for (k = 0; k < c; k++) prop[k] = (metodo == METODO_SEM_SOMA1) ? melhor_w[k] : melhor_w[k] / soma;
    }
    free(coef);
    free(w);
    free(melhor_w);
    free(reconstruido);
    return melhor;
}
