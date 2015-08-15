/*
 * threadfisher.h
 *
 * Author: tuvakt
 */

#ifndef THREADFISHER_H_
#define THREADFISHER_H_

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
    double perc;
	int *apos;
	int *bpos;
	double *avals;
	double *bvals;
	double *scores;
	double *stddev;
};

void threadcompute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen,
		   double perc, double *scores, double *stddev);
void mycompute(void *threadarg);

#endif /* THREADFISHER_H_ */
