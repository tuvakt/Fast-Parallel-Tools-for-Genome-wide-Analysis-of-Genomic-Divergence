/*
 * css.h
 *
 * Author: tuvakt
 */

#ifndef FILE_CSSS_SEEN
#define FILE_CSSS_SEEN

void compute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, int treshold, int runs,
		int drosophila, int mdsalg, double *scores, double *p);
double cluster_separation_scorer(double **distance, double *avals, double* bvals,  int *atracks,  int *btracks, int asize, int bsize, int npos, double **dissimilarity,
		int drosophila, int mdsalg, double **X,double **B, double **Z, double **tmp, double **L, double **Q, double **result);
void compare_all(double *avals, double *bvals, int asize, int bsize, int npos,  double **dissimilarity);
void compare_freq(double *avals, double *bvals, int npos, double **dissimilarity);
int fill_averages(double **dissimilarity, int m);
void allocate_matrix(double ***A, int m, int n);
void deallocate_matrix(double **A);
void matrix_mult(double **A, double **B, double **C, int m, int n, int p);
void matrix_squared(double **A, double **B, int m);
void matrix_sqrt(double **A, int m);
void setup_z_matrix(double **A, int m);
void cmds( double **compared, double **X, int dims, int m, double **B, double **Z, double **tmp, double **L, double **Q);
void calc_dist(double **A, double **distance, int m);
double css(double **distance,  int *atracks,  int *btracks, int a, int b);
void random_shuffle( int *elms, int n, unsigned short state[]);
double significance_treshold(double **distance,  int *tracks, int asize, int bsize, double score, int treshold, int runs, unsigned short state[]);
double smacof(double **dissimilarity, int m, int dims, double **X, double **Z, double **B, double **D, int max_iters, double epsilon);
void smacof_runs(double **dissimilarity, int m, int dims, double **X, double **Z, double **B, double **D, double **result, int max_iters, int n_init, double epsilon);


#endif
