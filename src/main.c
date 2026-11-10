/*
 * main.c - Ponto de entrada: lê argumentos, carrega as tabelas e roda o pipeline.
 */
#include "cli.h"
#include "pipeline.h"
#include "tabela.h"

int main(int argc, char **argv)
{
    Configuracao cfg;
    Tabela assinatura, mistura;

    if (!ler_argumentos(argc, argv, &cfg)) return 0;

    assinatura = ler_tabela(cfg.arquivo_assinatura);
    mistura = ler_tabela(cfg.arquivo_mistura);
    executar_cibersort(&assinatura, &mistura, &cfg.opcoes, cfg.arquivo_saida);
    liberar_tabela(&assinatura);
    liberar_tabela(&mistura);
    return 0;
}
