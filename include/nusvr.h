/*
 * nusvr.h - ν-SVR linear (o coração do método).
 *
 * Reescrita do algoritmo da LIBSVM (Chang & Lin) para a ν-SVR, o mesmo que o
 * pacote e1071 do R usa por baixo dos panos no CIBERSORT original.
 *
 * Problema: dada a assinatura X (l genes x c tipos celulares) e a expressão y
 * de uma amostra, achar os pesos w (c valores) com  y ≈ X w.
 *
 * Forma dual. 2l variáveis alfa (as l primeiras são α_i,
 * as l últimas são α*_i):
 *      minimizar   (1/2) αᵀ Q α + pᵀ α
 *      sujeito a   0 <= α_k <= C
 *                  soma dos α_i (k < l)  =  soma dos α*_i (k >= l)  =  C ν l / 2
 *   com Q_km = s_k s_m K(k mod l, m mod l),  s = +1 para k < l e -1 para k >= l,
 *   K = X Xᵀ (núcleo linear),  p_k = -y_k (k < l)  e  +y_(k-l) (k >= l).
 * No fim, coef_i = α_i - α*_i  e  w = Σ_i coef_i x_i.
 *
 * Método SMO: a cada passo escolhe dois α do mesmo grupo, resolve só esses
 * dois por fórmula fechada (a soma do grupo não muda) e atualiza o gradiente.
 * Para quando a maior violação das condições ótimas (KKT) é < tolerância.
 */
#ifndef CIBERSORT_NUSVR_H
#define CIBERSORT_NUSVR_H

typedef struct {
    double rho;        /* -intercepto (o CIBERSORT não usa) */
    double epsilon;    /* largura do tubo encontrada pelo ajuste */
    long iteracoes;    /* passos do SMO executados */
    int convergiu;     /* 1 se parou pelo critério de tolerância */
} InfoNuSVR;

/*
 * Matriz de núcleo linear: K[i*l + j] = produto escalar das linhas i e j de X
 * (X: l x c, linha por linha). Devolve l x l; liberar com free().
 */
double *matriz_nucleo_linear(const double *X, int l, int c);

/*
 * Treina a ν-SVR.
 * Entradas: K (l x l), y (l valores), nu, C, tolerância, limite de iterações.
 * Saídas:   coef[0..l-1] (coef_i = α_i - α*_i) e as informações em 'info'.
 */
void treinar_nusvr(const double *K, int l, const double *y, double nu, double C, double tolerancia, long max_iteracoes, double *coef, InfoNuSVR *info);

#endif /* CIBERSORT_NUSVR_H */
