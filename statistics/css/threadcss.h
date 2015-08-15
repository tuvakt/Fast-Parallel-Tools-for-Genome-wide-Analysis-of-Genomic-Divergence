/*
 * threadscss.h
 *
 * Author: tuvakt
 */

#ifndef THREADSCSS_H_
#define THREADSCSS_H_

#define FALSE 0
#define TRUE 1

/**
 * The data needed by the threads
 */
struct thread_data{
	int thread_id;
	int num_windows;
	int regend;
	int wsize;
	int wstep;
	int alen;
	int blen;
	int treshold;
	int runs;
	int drosophila;
	int mds;
	int *apos;
	int *bpos;
	double *avals;
	double *bvals;
	double *scores;
	double *p;
};

void threadcompute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen,
		int treshold, int runs, int drosophila, int mds, double *scores, double *p);
void mycompute(void *threadarg);

#endif /* THREADSCSS_H_ */
