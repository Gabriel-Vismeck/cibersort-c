/*
 * util.c - Utilidades gerais: erro fatal e alocação de memória segura.
 */
#include "util.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void morrer(const char *formato, ...)
{
    va_list args;
    va_start(args, formato);
    fprintf(stderr, "ERRO: ");
    vfprintf(stderr, formato, args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(1);
}

void *alocar(size_t bytes)
{
    void *p = malloc(bytes > 0 ? bytes : 1);
    if (!p) morrer("memória insuficiente (%zu bytes).", bytes);
    return p;
}

char *copiar_texto(const char *s)
{
    size_t n = strlen(s);
    char *c = (char *)alocar(n + 1);
    memcpy(c, s, n + 1);
    return c;
}
