/*
 * central.h - Algoritmo central do CIBERSORT para uma amostra.
 *
 * Para cada ν em {0,25 ; 0,5 ; 0,75}:
 *   1. treina a ν-SVR e obtém os pesos w (um por tipo celular);
 *   2. pesos negativos viram 0;
 *   3. divide pela soma (proporções somam 1);
 *   4. reconstrói a amostra k = X w e mede RMSE e correlação com y.
 * Fica com o ν de menor RMSE.
 */
#ifndef CIBERSORT_CENTRAL_H
#define CIBERSORT_CENTRAL_H

#define METODO_RELATIVO 0  /* proporções somam 1 (padrão) */
#define METODO_SEM_SOMA1 1 /* modo absoluto "no.sumto1": não divide pela soma */

typedef struct {
    double rmse;       /* erro da melhor reconstrução */
    double correlacao; /* correlação de Pearson entre reconstrução e amostra */
    double nu;         /* ν vencedor */
} ResultadoCentral;

/*
 * X: l x c (linha por linha), padronizada;  K : núcleo de X;  y : l valores.
 * metodo: METODO_RELATIVO ou METODO_SEM_SOMA1.
 * prop: saída, c proporções. Se nenhum ν der solução útil, prop = NAN.
 */
ResultadoCentral algoritmo_central(const double *X, const double *K, int l, int c, const double *y, int metodo, double tolerancia, double *prop);

#endif /* CIBERSORT_CENTRAL_H */
