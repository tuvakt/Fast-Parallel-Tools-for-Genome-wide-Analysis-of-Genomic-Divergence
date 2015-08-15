/*
 * cFisher.h
 *
 * Author: tuvakt
 */


#ifndef FILE_CSSS_SEEN
#define FILE_CSSS_SEEN

void compute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, 
	     double perc, double *scores, double *stddev);
void fisher_exact_test(double *results, double *avals, double *bvals, int asize, int bsize, int npos, int *f, int *tmp, double *samples,
		double *stdsamples, int nsamples, double *fetscores, unsigned short state[], double perc);
void fetcount(int *f, double *avals, double *bvals, int idx, int asize, int bsize);
unsigned long binomial(unsigned long n, unsigned long k);
double fet(int *counts, int *tmp);
double fet_p(int a, int b, int c, int d);
double std(double *a, int n);
double mean(double *a, int n);
void shift_table(int *f, int *tmp, int n);
void create_table(int *f, int n);
int min_idx(int *a, int n);
double calc_std(double *fetscores, double *samples, double *stdsamples, int nsamples, int npos, double perc, unsigned short state[]);
#endif
