/*
 * aleatorio.c - xoshiro256** + sorteio de permutações. Veja aleatorio.h.
 */
#include "aleatorio.h"

#include <stdlib.h>

#include "util.h"

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void iniciar_gerador(Gerador *g, uint64_t semente)
{
    int i;
    for (i = 0; i < 4; i++) g->s[i] = splitmix64(&semente);
}

static uint64_t girar(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

uint64_t proximo_aleatorio(Gerador *g)
{
    uint64_t *s = g->s;
    uint64_t resultado = girar(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = girar(s[3], 45);
    return resultado;
}

uint64_t inteiro_ate(Gerador *g, uint64_t limite)
{
    uint64_t corte = (0ULL - limite) % limite, r;
    do {
        r = proximo_aleatorio(g);
    } while (r < corte);
    return r % limite;
}

int *sortear_indices(int n_perm, int quantos, int total, uint64_t semente)
{
    Gerador g;
    int *todos = (int *)alocar((size_t)total * sizeof(int));
    int *saida = (int *)alocar((size_t)n_perm * (size_t)quantos * sizeof(int));
    int p, i;
    iniciar_gerador(&g, semente);
    for (i = 0; i < total; i++) todos[i] = i;
    for (p = 0; p < n_perm; p++) {
        for (i = 0; i < quantos; i++) {
            int j = i + (int)inteiro_ate(&g, (uint64_t)(total - i));
            int tmp = todos[i];
            todos[i] = todos[j];
            todos[j] = tmp;
            saida[(size_t)p * quantos + i] = todos[i];
        }
    }
    free(todos);
    return saida;
}
