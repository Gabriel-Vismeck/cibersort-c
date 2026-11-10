/*
 * quantis.h - Normalização por quantis (Bolstad et al., 2003).
 */
#ifndef CIBERSORT_QUANTIS_H
#define CIBERSORT_QUANTIS_H

/*
 * Normaliza por quantis, no lugar. 'm' guarda uma amostra por vez, contígua:
 * m[amostra * n_genes + gene].
 *   1. ordena cada amostra.
 *   2. o "alvo" de cada posição é a média dessa posição entre as amostras.
 *   3. cada gene recebe o alvo do seu posto; empates recebem a média dos alvos das posições que ocupam (como o preprocessCore do R).
 */
void normalizar_quantis(double *m, int n_genes, int n_amostras);

#endif /* CIBERSORT_QUANTIS_H */
