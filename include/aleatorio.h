/*
 * aleatorio.h - Gerador de números aleatórios e sorteio das permutações.
 *
 * Gerador xoshiro256**, iniciado com splitmix64.
 * Atenção: o gerador do R é outro. Por isso o valor-p muda um pouco entre
 * sementes diferentes (e entre R e C), mas as proporções estimadas são
 * determinísticas e não dependem do sorteio.
 */
#ifndef CIBERSORT_ALEATORIO_H
#define CIBERSORT_ALEATORIO_H

#include <stdint.h>

typedef struct {
    uint64_t s[4];
} Gerador;

/* Inicia o gerador a partir de uma semente. */
void iniciar_gerador(Gerador *g, uint64_t semente);

/* Próximo número de 64 bits. */
uint64_t proximo_aleatorio(Gerador *g);

/* Inteiro uniforme em [0, limite), sem viés (rejeição). 'limite' deve ser > 0. */
uint64_t inteiro_ate(Gerador *g, uint64_t limite);

/*
 * Sorteia, para cada permutação, 'quantos' índices distintos entre 0 e total-1
 * (embaralhamento parcial de Fisher-Yates).
 * Resultado: vetor n_perm x quantos (linha por linha); liberar com free().
 */
int *sortear_indices(int n_perm, int quantos, int total, uint64_t semente);

#endif /* CIBERSORT_ALEATORIO_H */
