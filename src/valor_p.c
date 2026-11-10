/*
 * valor_p.c - Valor-p a partir da distribuição nula. Veja valor_p.h.
 */
#include "valor_p.h"

#include <math.h>

double valor_p(const double *nula, int n, double r, int tipo)
{
    int i;
    if (n == 0) return 9999.0; /* sem permutações: o CIBERSORT usa 9999 */
    if (tipo == P_ORIGINAL) {
        int melhor = 0;
        double menor = fabs(nula[0] - r);
        for (i = 1; i < n; i++) {
            double d = fabs(nula[i] - r);
            if (d < menor) { /* which.min: primeiro de menor diferença */
                menor = d;
                melhor = i;
            }
        }
        return 1.0 - (double)(melhor + 1) / (double)n;
    } else {
        int maiores = 0;
        for (i = 0; i < n; i++)
            if (nula[i] >= r) maiores++;
        return (1.0 + maiores) / (1.0 + n);
    }
}
