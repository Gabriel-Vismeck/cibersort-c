/*
 * pipeline.h - Pipeline completo do CIBERSORT.
 *
 *   Matriz de expressão (genes x amostras)
 *     -> pré-processamento (anti-log, normalização por quantis)
 *     -> seleção de genes em comum + padronização da assinatura (ex.: LM22)
 *     -> distribuição nula por permutação (confiabilidade)
 *     -> ν-SVR em cada amostra (deconvolução)
 *     -> arquivo TSV com as proporções celulares estimadas
 */
#ifndef CIBERSORT_PIPELINE_H
#define CIBERSORT_PIPELINE_H

#include "opcoes.h"
#include "tabela.h"

/*
 * Executa o CIBERSORT e grava o resultado (TSV) em 'arquivo_saida'.
 * Colunas da saída: Mistura | uma por tipo celular | Valor_p | Correlacao | RMSE | [Escore_absoluto, se op->absoluto].
 */
void executar_cibersort(const Tabela *assinatura, const Tabela *mistura, const Opcoes *op, const char *arquivo_saida);

#endif /* CIBERSORT_PIPELINE_H */
