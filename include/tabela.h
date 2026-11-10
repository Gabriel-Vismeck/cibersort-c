/*
 * tabela.h - Leitura de tabelas de texto (TSV/CSV) com genes nas linhas.
 *
 * Formato aceito (separador tab, ';' ou ','):
 *     Gene symbol <TAB> Amostra1 <TAB> Amostra2 ...
 *     ABCB4       <TAB> 12.5     <TAB> 30.1    ...
 * Também aceita o formato do write.table() do R, em que o cabeçalho tem uma
 * coluna a menos (sem a célula do canto).
 * Limitação: não trata separador dentro de aspas.
 */
#ifndef CIBERSORT_TABELA_H
#define CIBERSORT_TABELA_H

typedef struct {
    int n_linhas;          /* genes */
    int n_colunas;         /* amostras ou tipos celulares */
    char **nomes_linhas;   /* n_linhas nomes */
    char **nomes_colunas;  /* n_colunas nomes */
    double *dados;         /* n_linhas x n_colunas, linha por linha: dados[i*n_colunas + j] */
} Tabela;

/* Lê a tabela do arquivo. Encerra o programa com mensagem se o arquivo for inválido. */
Tabela ler_tabela(const char *caminho);

/* Libera toda a memória da tabela (não libera o ponteiro 't' em si). */
void liberar_tabela(Tabela *t);

#endif /* CIBERSORT_TABELA_H */
