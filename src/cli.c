/*
 * cli.c - Interpretação da linha de comando. Veja cli.h.
 */
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "central.h"
#include "util.h"
#include "valor_p.h"

void mostrar_ajuda(void)
{
    printf("CIBERSORT em C puro (em português)\n\n"
           "Uso: cibersort --assinatura LM22.txt --mistura expressao.txt --saida resultado.tsv [opções]\n\n"
           "  --assinatura ARQ        tabela da assinatura (genes x tipos celulares), escala linear\n"
           "  --expressao ARQ           tabela de expressão (genes x amostras)\n"
           "  --saida ARQ             arquivo de resultado (.tsv)\n"
           "  --perm N                número de permutações (padrão 100; 0 = sem valor-p)\n"
           "  --com-quantis           normaliza por quantis (recomendado para microarray)\n"
           "  --sem-quantis           não normaliza por quantis (padrão recomendado para RNA-seq)\n"
           "  --absoluto              modo absoluto\n"
           "  --metodo-absoluto M     sig.score (padrão) ou no.sumto1\n"
           "  --semente N             semente do sorteio das permutações (padrão 1)\n"
           "  --valor-p T             original (padrão, igual ao CIBERSORT) ou classico\n"
           "  --tolerancia X          tolerância do ν-SVR (padrão 0.001, a mesma do CIBERSORT original)\n"
           "  --threads N             número de núcleos (só se compilado com OpenMP)\n");
}

/* Valores padrão de todas as opções. */
static void definir_padroes(Configuracao *cfg)
{
    cfg->arquivo_assinatura = NULL;
    cfg->arquivo_mistura = NULL;
    cfg->arquivo_saida = NULL;
    cfg->opcoes.n_perm = 100;
    cfg->opcoes.normalizar_por_quantis = 0;
    cfg->opcoes.absoluto = 0;
    cfg->opcoes.metodo = METODO_RELATIVO;
    cfg->opcoes.semente = 1;
    cfg->opcoes.tipo_p = P_ORIGINAL;
    cfg->opcoes.tolerancia = 1e-3;
}

int ler_argumentos(int argc, char **argv, Configuracao *cfg)
{
    int i;
    Opcoes *op = &cfg->opcoes;

    definir_padroes(cfg);
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--ajuda") || !strcmp(a, "-h") || !strcmp(a, "--help")) {
            mostrar_ajuda();
            return 0;
        } 
        else if (!strcmp(a, "--sem-quantis")) op->normalizar_por_quantis = 0;
        else if (!strcmp(a, "--com-quantis")) op->normalizar_por_quantis = 1;
        else if (!strcmp(a, "--absoluto")) op->absoluto = 1;
        else if (i + 1 < argc && !strcmp(a, "--assinatura")) cfg->arquivo_assinatura = argv[++i];
        else if (i + 1 < argc && (!strcmp(a, "--expressao") || !strcmp(a, "--mistura"))) cfg->arquivo_mistura = argv[++i];
        else if (i + 1 < argc && !strcmp(a, "--saida")) cfg->arquivo_saida = argv[++i];
        else if (i + 1 < argc && !strcmp(a, "--perm")) op->n_perm = atoi(argv[++i]);
        else if (i + 1 < argc && !strcmp(a, "--tolerancia")) op->tolerancia = atof(argv[++i]);
        else if (i + 1 < argc && !strcmp(a, "--semente")) op->semente = (uint64_t)strtoull(argv[++i], NULL, 10);
        else if (i + 1 < argc && !strcmp(a, "--metodo-absoluto")) {
            const char *v = argv[++i];
            if (!strcmp(v, "sig.score")) op->metodo = METODO_RELATIVO;
            else if (!strcmp(v, "no.sumto1")) op->metodo = METODO_SEM_SOMA1;
            else morrer("--metodo-absoluto deve ser sig.score ou no.sumto1.");
        } else if (i + 1 < argc && !strcmp(a, "--valor-p")) {
            const char *v = argv[++i];
            if (!strcmp(v, "original")) op->tipo_p = P_ORIGINAL;
            else if (!strcmp(v, "classico")) op->tipo_p = P_CLASSICO;
            else morrer("--valor-p deve ser original ou classico.");
        } else if (i + 1 < argc && !strcmp(a, "--threads")) {
#ifdef _OPENMP
            omp_set_num_threads(atoi(argv[++i]));
#else
            ++i;
            fprintf(stderr, "Aviso: compilado sem OpenMP; --threads ignorado.\n");
#endif
        } else {
            mostrar_ajuda();
            morrer("opção desconhecida ou sem valor: %s", a);
        }
    }
    if (!cfg->arquivo_assinatura || !cfg->arquivo_mistura || !cfg->arquivo_saida) {
        mostrar_ajuda();
        morrer("faltam --assinatura, --mistura e/ou --saida.");
    }
    if (op->n_perm < 0) morrer("--perm não pode ser negativo.");
    return 1;
}
