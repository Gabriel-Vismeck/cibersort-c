/*
 * cli.h - Interpretação da linha de comando.
 */
#ifndef CIBERSORT_CLI_H
#define CIBERSORT_CLI_H

#include "opcoes.h"

typedef struct {
    const char *arquivo_assinatura;
    const char *arquivo_mistura;
    const char *arquivo_saida;
    Opcoes opcoes;
} Configuracao;

/* Imprime o texto de ajuda em stdout. */
void mostrar_ajuda(void);

/*
 * Interpreta os argumentos. Devolve 1 se o programa deve seguir normalmente,
 * 0 se só a ajuda foi pedida (o chamador deve encerrar com sucesso).
 * Argumentos inválidos encerram o programa com mensagem de erro.
 */
int ler_argumentos(int argc, char **argv, Configuracao *cfg);

#endif /* CIBERSORT_CLI_H */
