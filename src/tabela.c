/*
 * tabela.c - Leitura de tabelas de texto (TSV/CSV). Veja tabela.h.
 */
#include "tabela.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

void liberar_tabela(Tabela *t)
{
    int i;
    for (i = 0; i < t->n_linhas; i++) free(t->nomes_linhas[i]);
    for (i = 0; i < t->n_colunas; i++) free(t->nomes_colunas[i]);
    free(t->nomes_linhas);
    free(t->nomes_colunas);
    free(t->dados);
}

/* Lê o arquivo inteiro para a memória (e remove a marca BOM do UTF-8, se houver). */
static char *ler_arquivo_inteiro(const char *caminho)
{
    FILE *f = fopen(caminho, "rb");
    long tamanho;
    char *buffer;
    if (!f) morrer("falha ao abrir o arquivo '%s'.", caminho);
    fseek(f, 0, SEEK_END);
    tamanho = ftell(f);
    fseek(f, 0, SEEK_SET);
    buffer = (char *)alocar((size_t)tamanho + 1);
    if (fread(buffer, 1, (size_t)tamanho, f) != (size_t)tamanho) morrer("erro ao ler '%s'.", caminho);
    buffer[tamanho] = '\0';
    fclose(f);
    if (tamanho >= 3 && (unsigned char)buffer[0] == 0xEF && (unsigned char)buffer[1] == 0xBB &&
        (unsigned char)buffer[2] == 0xBF)
        memmove(buffer, buffer + 3, (size_t)tamanho - 2);
    return buffer;
}

/* Remove espaços e aspas das pontas de um campo (altera o texto no lugar). */
static char *aparar(char *s)
{
    char *fim;
    while (*s == ' ' || *s == '\t' || *s == '"') s++;
    fim = s + strlen(s);
    while (fim > s && (fim[-1] == ' ' || fim[-1] == '"' || fim[-1] == '\r')) fim--;
    *fim = '\0';
    return s;
}

/* Divide uma linha em campos separados por 'sep'. Devolve quantos campos há. */
static int separar_campos(char *linha, char sep, char **campos, int maximo)
{
    int n = 0;
    char *inicio = linha;
    char *p;
    for (p = linha;; p++) {
        if (*p == sep || *p == '\0') {
            int fim_da_linha = (*p == '\0');
            *p = '\0';
            if (n < maximo) campos[n] = inicio;
            n++;
            if (fim_da_linha) break;
            inicio = p + 1;
        }
    }
    return n;
}

/* Conta quantos campos há numa linha, sem alterá-la. */
static int contar_campos(const char *linha, char sep)
{
    int n = 1;
    for (; *linha; linha++)
        if (*linha == sep) n++;
    return n;
}

static int linha_em_branco(const char *s)
{
    while (*s) {
        if (*s != ' ' && *s != '\t' && *s != '\r') return 0;
        s++;
    }
    return 1;
}

Tabela ler_tabela(const char *caminho)
{
    Tabela t;
    char *buffer = ler_arquivo_inteiro(caminho);
    char **linhas;
    int n_texto = 0, n_uteis = 0, i, j, n_campos_dados, n_campos_cab, max_campos;
    char *p;
    char sep;
    char **campos;

    /* Quebra o texto em linhas (troca '\n' por '\0'), ignorando linhas em branco. */
    for (p = buffer; *p; p++)
        if (*p == '\n') n_texto++;
    linhas = (char **)alocar(((size_t)n_texto + 2) * sizeof(char *));
    {
        char *inicio = buffer;
        for (p = buffer;; p++) {
            if (*p == '\n' || *p == '\0') {
                int fim = (*p == '\0');
                *p = '\0';
                if (!linha_em_branco(inicio)) linhas[n_uteis++] = inicio;
                if (fim) break;
                inicio = p + 1;
            }
        }
    }
    if (n_uteis < 2) morrer("o arquivo '%s' precisa ter cabeçalho e pelo menos uma linha de dados.", caminho);

    /* Separador: tab, ponto e vírgula ou vírgula (olhando o cabeçalho). */
    if (strchr(linhas[0], '\t')) sep = '\t';
    else if (strchr(linhas[0], ';') && !strchr(linhas[0], ',')) sep = ';';
    else sep = ',';

    /* Conta os campos do cabeçalho e da primeira linha de dados (sem alterar o texto). */
    n_campos_cab = contar_campos(linhas[0], sep);
    n_campos_dados = contar_campos(linhas[1], sep);
    max_campos = (n_campos_cab > n_campos_dados ? n_campos_cab : n_campos_dados) + 2;
    campos = (char **)alocar((size_t)max_campos * sizeof(char *));
    separar_campos(linhas[0], sep, campos, max_campos); /* corta o cabeçalho em campos */
    t.n_linhas = n_uteis - 1;
    if (n_campos_cab == n_campos_dados - 1) { /* formato do write.table do R */
        t.n_colunas = n_campos_cab;
        t.nomes_colunas = (char **)alocar((size_t)t.n_colunas * sizeof(char *));
        for (j = 0; j < t.n_colunas; j++) t.nomes_colunas[j] = copiar_texto(aparar(campos[j]));
    } else if (n_campos_cab == n_campos_dados) { /* cabeçalho com a célula do canto */
        t.n_colunas = n_campos_cab - 1;
        t.nomes_colunas = (char **)alocar((size_t)t.n_colunas * sizeof(char *));
        for (j = 0; j < t.n_colunas; j++) t.nomes_colunas[j] = copiar_texto(aparar(campos[j + 1]));
    } else {
        morrer("'%s': o cabeçalho tem %d campos e a primeira linha de dados tem %d.", caminho,
               n_campos_cab, n_campos_dados);
    }

    t.nomes_linhas = (char **)alocar((size_t)t.n_linhas * sizeof(char *));
    t.dados = (double *)alocar((size_t)t.n_linhas * (size_t)t.n_colunas * sizeof(double));
    for (i = 0; i < t.n_linhas; i++) {
        int n = separar_campos(linhas[i + 1], sep, campos, max_campos);
        if (n != t.n_colunas + 1)
            morrer("'%s', linha %d: %d campos (esperava %d).", caminho, i + 2, n, t.n_colunas + 1);
        t.nomes_linhas[i] = copiar_texto(aparar(campos[0]));
        for (j = 0; j < t.n_colunas; j++) {
            char *texto = aparar(campos[j + 1]);
            char *fim;
            double v = strtod(texto, &fim);
            if (fim == texto || *fim != '\0')
                morrer("valor não numérico '%s' no gene '%s', coluna '%s'.", texto, t.nomes_linhas[i],
                       t.nomes_colunas[j]);
            t.dados[(size_t)i * (size_t)t.n_colunas + (size_t)j] = v;
        }
    }
    free(campos);
    free(linhas);
    free(buffer);
    return t;
}
