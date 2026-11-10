/*
 * opcoes.h - Opções de execução do CIBERSORT (preenchidas pela linha de comando).
 */
#ifndef CIBERSORT_OPCOES_H
#define CIBERSORT_OPCOES_H

#include <stdint.h>

typedef struct {
    int n_perm;                 /* nº de permutações (0 = sem valor-p) */
    int normalizar_por_quantis; /* 1 = normaliza a mistura por quantis */
    int absoluto;               /* 1 = modo absoluto */
    int metodo;                 /* METODO_RELATIVO ou METODO_SEM_SOMA1 (central.h) */
    uint64_t semente;           /* semente do sorteio das permutações */
    int tipo_p;                 /* P_ORIGINAL ou P_CLASSICO (valor_p.h) */
    double tolerancia;          /* critério de parada do ν-SVR (padrão 0,001, como no e1071) */
} Opcoes;

#endif /* CIBERSORT_OPCOES_H */
