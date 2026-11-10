/*
 * nusvr.c - ν-SVR linear (solver SMO estilo LIBSVM). Veja nusvr.h para a teoria.
 */
#include "nusvr.h"

#include <math.h>
#include <stdlib.h>

#include "util.h"

#define TAU 1e-12            /* no lugar de coeficientes quadráticos <= 0 */
#define NO_LIMITE_INFERIOR 0 /* α = 0 */
#define NO_LIMITE_SUPERIOR 1 /* α = C */
#define LIVRE 2              /* 0 < α < C */

double *matriz_nucleo_linear(const double *X, int l, int c)
{
    double *K = (double *)alocar((size_t)l * (size_t)l * sizeof(double));
    int i, j, k;
    for (i = 0; i < l; i++)
        for (j = i; j < l; j++) {
            double v = 0.0;
            for (k = 0; k < c; k++) v += X[(size_t)i * c + k] * X[(size_t)j * c + k];
            K[(size_t)i * l + j] = v;
            K[(size_t)j * l + i] = v;
        }
    return K;
}

static unsigned char situacao_do_alfa(double alfa, double C)
{
    if (alfa >= C) return NO_LIMITE_SUPERIOR;
    if (alfa <= 0.0) return NO_LIMITE_INFERIOR;
    return LIVRE;
}

void treinar_nusvr(const double *K, int l, const double *y, double nu, double C, double tolerancia, long max_iteracoes, double *coef, InfoNuSVR *info)
{
    int n = 2 * l, k, m;
    double *alfa = (double *)alocar((size_t)n * sizeof(double));
    double *G = (double *)alocar((size_t)n * sizeof(double)); /* gradiente Q α + p */
    unsigned char *estado = (unsigned char *)alocar((size_t)n);
    double *diagonal = (double *)alocar((size_t)l * sizeof(double)); /* K_ii */
    double *Kdif = (double *)alocar((size_t)l * sizeof(double));
    double soma;
    long iteracoes = 0;
    int convergiu = 0;

    for (k = 0; k < l; k++) diagonal[k] = K[(size_t)k * l + k];

    /* Ponto de partida viável (como na LIBSVM): distribui C*ν*l/2 em cada grupo. */
    soma = C * nu * (double)l / 2.0;
    for (k = 0; k < l; k++) {
        double a = soma < C ? soma : C;
        alfa[k] = a;
        alfa[k + l] = a;
        soma -= a;
    }

    /* Gradiente inicial: G = p + Q α, com (Qα)_k = s_k Σ_m K[k][m](α_m - α*_m). */
    for (k = 0; k < l; k++) {
        double s = 0.0;
        for (m = 0; m < l; m++) s += K[(size_t)k * l + m] * (alfa[m] - alfa[m + l]);
        Kdif[k] = s;
    }
    for (k = 0; k < l; k++) {
        G[k] = -y[k] + Kdif[k];
        G[k + l] = y[k] - Kdif[k];
    }
    for (k = 0; k < n; k++) estado[k] = situacao_do_alfa(alfa[k], C);

    while (iteracoes < max_iteracoes) {
        double gmax_p = -INFINITY, gmax_p2 = -INFINITY, gmax_n = -INFINITY, gmax_n2 = -INFINITY;
        int ip = -1, in_ = -1, j_escolhido = -1, i, j;
        double menor_queda = INFINITY;
        const double *linha_ip = NULL, *linha_in = NULL, *linha_i, *linha_j;
        double coef_quad, delta, soma_ij, antigo_i, antigo_j, a_i, a_j;

        /* ===== 1) Escolha de i: a variável que mais viola as condições ótimas
         *          em cada grupo ("+" = 0..l-1 ; "-" = l..2l-1). ===== */
        for (k = 0; k < l; k++)
            if (estado[k] != NO_LIMITE_SUPERIOR && -G[k] >= gmax_p) {
                gmax_p = -G[k];
                ip = k;
            }
        for (k = l; k < n; k++)
            if (estado[k] != NO_LIMITE_INFERIOR && G[k] >= gmax_n) {
                gmax_n = G[k];
                in_ = k;
            }
        if (ip != -1) linha_ip = K + (size_t)ip * l;
        if (in_ != -1) linha_in = K + (size_t)(in_ - l) * l;

        /* ===== 1b) Escolha de j: a que mais reduz o objetivo (2a ordem). ===== */
        if (ip != -1) {
            for (j = 0; j < l; j++) {
                if (estado[j] != NO_LIMITE_INFERIOR) {
                    double dif_grad = gmax_p + G[j];
                    if (G[j] >= gmax_p2) gmax_p2 = G[j];
                    if (dif_grad > 0) {
                        double queda;
                        coef_quad = diagonal[ip] + diagonal[j] - 2.0 * linha_ip[j];
                        queda = (coef_quad > 0) ? -(dif_grad * dif_grad) / coef_quad
                                                : -(dif_grad * dif_grad) / TAU;
                        if (queda <= menor_queda) {
                            j_escolhido = j;
                            menor_queda = queda;
                        }
                    }
                }
            }
        }
        if (in_ != -1) {
            for (j = l; j < n; j++) {
                if (estado[j] != NO_LIMITE_SUPERIOR) {
                    double dif_grad = gmax_n - G[j];
                    if (-G[j] >= gmax_n2) gmax_n2 = -G[j];
                    if (dif_grad > 0) {
                        double queda;
                        coef_quad = diagonal[in_ - l] + diagonal[j - l] - 2.0 * linha_in[j - l];
                        queda = (coef_quad > 0) ? -(dif_grad * dif_grad) / coef_quad
                                                : -(dif_grad * dif_grad) / TAU;
                        if (queda <= menor_queda) {
                            j_escolhido = j;
                            menor_queda = queda;
                        }
                    }
                }
            }
        }

        /* ===== 2) Critério de parada ===== */
        {
            double v1 = gmax_p + gmax_p2, v2 = gmax_n + gmax_n2;
            double violacao = v1 > v2 ? v1 : v2;
            if (violacao < tolerancia || j_escolhido == -1) {
                convergiu = 1;
                break;
            }
        }

        i = (j_escolhido < l) ? ip : in_;
        j = j_escolhido;
        iteracoes++;

        /* ===== 3) Passo analítico em (α_i, α_j): mesma soma, dentro de [0, C]. ===== */
        linha_i = K + (size_t)(i % l) * l;
        linha_j = K + (size_t)(j % l) * l;
        coef_quad = diagonal[i % l] + diagonal[j % l] - 2.0 * linha_i[j % l];
        if (coef_quad <= 0) coef_quad = TAU;
        delta = (G[i] - G[j]) / coef_quad;
        soma_ij = alfa[i] + alfa[j];
        antigo_i = alfa[i];
        antigo_j = alfa[j];
        alfa[i] -= delta;
        alfa[j] += delta;

        if (soma_ij > C) {
            if (alfa[i] > C) {
                alfa[i] = C;
                alfa[j] = soma_ij - C;
            }
        } else {
            if (alfa[j] < 0) {
                alfa[j] = 0.0;
                alfa[i] = soma_ij;
            }
        }
        if (soma_ij > C) {
            if (alfa[j] > C) {
                alfa[j] = C;
                alfa[i] = soma_ij - C;
            }
        } else {
            if (alfa[i] < 0) {
                alfa[i] = 0.0;
                alfa[j] = soma_ij;
            }
        }

        /* ===== 4) Atualiza o gradiente: G_k += Q_ik Δα_i + Q_jk Δα_j. =====
         * Como Q_ik = s_i s_k K, as duas metades de G recebem a mesma quantidade
         * com sinais opostos. */
        a_i = ((i < l) ? 1.0 : -1.0) * (alfa[i] - antigo_i);
        a_j = ((j < l) ? 1.0 : -1.0) * (alfa[j] - antigo_j);
        for (m = 0; m < l; m++) {
            double t = linha_i[m] * a_i + linha_j[m] * a_j;
            G[m] += t;
            G[m + l] -= t;
        }
        estado[i] = situacao_do_alfa(alfa[i], C);
        estado[j] = situacao_do_alfa(alfa[j], C);
    }

    /* ---- rho e epsilon (calculate_rho do Solver_NU da LIBSVM) ---- */
    {
        int livres1 = 0, livres2 = 0;
        double ub1 = INFINITY, ub2 = INFINITY, lb1 = -INFINITY, lb2 = -INFINITY;
        double soma_livres1 = 0.0, soma_livres2 = 0.0, r1, r2;
        for (k = 0; k < l; k++) {
            if (estado[k] == NO_LIMITE_SUPERIOR) { if (G[k] > lb1) lb1 = G[k]; }
            else if (estado[k] == NO_LIMITE_INFERIOR) { if (G[k] < ub1) ub1 = G[k]; }
            else { livres1++; soma_livres1 += G[k]; }
        }
        for (k = l; k < n; k++) {
            if (estado[k] == NO_LIMITE_SUPERIOR) { if (G[k] > lb2) lb2 = G[k]; }
            else if (estado[k] == NO_LIMITE_INFERIOR) { if (G[k] < ub2) ub2 = G[k]; }
            else { livres2++; soma_livres2 += G[k]; }
        }
        r1 = livres1 > 0 ? soma_livres1 / livres1 : (ub1 + lb1) / 2.0;
        r2 = livres2 > 0 ? soma_livres2 / livres2 : (ub2 + lb2) / 2.0;
        info->epsilon = -((r1 + r2) / 2.0);
        info->rho = (r1 - r2) / 2.0;
    }
    for (k = 0; k < l; k++) coef[k] = alfa[k] - alfa[k + l];
    info->iteracoes = iteracoes;
    info->convergiu = convergiu;

    free(alfa);
    free(G);
    free(estado);
    free(diagonal);
    free(Kdif);
}
