/*
 * util.h - Utilidades gerais: erro fatal e alocação de memória segura.
 */
#ifndef CIBERSORT_UTIL_H
#define CIBERSORT_UTIL_H

#include <stddef.h>

/* Escreve "ERRO: <mensagem>" em stderr (formato printf) e encerra o programa. */
void morrer(const char *formato, ...);

/* malloc que nunca devolve NULL: se faltar memória, encerra com erro. */
void *alocar(size_t bytes);

/* Cópia de um texto (strdup não existe no C99 puro). Liberar com free(). */
char *copiar_texto(const char *s);

#endif /* CIBERSORT_UTIL_H */
