/*
 * threadcss.c
 *
 * Author: tuvakt
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "threadcss.h"
#include "css.h"
#include "comparative.h"

#define NUM_THREADS 64   							// number of threads
#define TASK_SIZE 100    							// number of windows per task

unsigned int task_id;    							// global task id
unsigned int num_tasks;      						// total number of tasks
pthread_mutex_t mutexTASK_ID;   					// mutex to protect 'task_id'
struct thread_data thread_data_array[NUM_THREADS];  // array for the thread data

/**
 * Computing the Cluster Separation Score (CSS) and the corresponding p-value for each window in a chromosome with threads.
 * The calculation is based on the article by Jones et al.
 * http://www.nature.com/nature/journal/v484/n7392/full/nature10944.html
 * and the master's thesis by Torkil Vederhus.
 *
 * @param avals a double array with values for group A, ordered by position and then individual id
 * @param bvals a double array with values for group B, ordered by position and then individual id
 * @param apos a double array with positions for group A, ordered by position and then individual id
 * @param bpos a double array with positions for group B, ordered by position and then individual id
 * @param regstart an integer region start
 * @param regend an integer, region end
 * @param wsize the size of the sliding window
 * @param wstep the step size of the sliding windows
 * @param asize the number of individuals in group A
 * @param bsize the number of individuals in group B
 * @param alen the length of the array avals and apos
 * @param blen the length of the array bvals and bpos
 * @param treshold the minimum amount of runs by the siginificance calculation
 * @param runs the maximum amount of runs by the significance calculation
 * @param drosophila an integer in [0,1] gives the choice of distance metric. 0: average of frequency, 1: count differences
 * @param mds an integer in [0,1,2] gives the choice of Multi-Dimensional Scaling (MDS) algorithm. 0: Classical MDS, 1: SMACOF, 2: SMACOF+CMDS
 * @param scores an double array with the resulting CSS scores for each window
 * @param p an double array with the resulting p-values for each window
 */
void threadcompute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen,
		   int treshold, int runs, int drosophila, int mds, double *scores, double *p) {
  
    struct timeval before, after;
    gettimeofday(&before, NULL);
    pthread_t threads[NUM_THREADS];
    
    int num_windows, rtn, i;
    num_windows = (regend/wstep)-3;
    printf("num windows %d\n", num_windows);
    num_tasks = num_windows/TASK_SIZE;
    printf("num tasks %d\n", num_tasks);
    
    // init mutex variable
    pthread_mutex_init(&mutexTASK_ID, NULL);
    task_id = 0;
    
    // create and start the threads
    for (i = 0; i < NUM_THREADS; i++) {
        thread_data_array[i].thread_id = i;
	thread_data_array[i].regend = regend;
	thread_data_array[i].num_windows = num_windows;
	thread_data_array[i].wsize = wsize;
	thread_data_array[i].wstep = wstep;
	thread_data_array[i].apos = apos;
	thread_data_array[i].bpos = bpos;
	thread_data_array[i].avals = avals;
	thread_data_array[i].bvals = bvals;
	thread_data_array[i].alen = alen;
	thread_data_array[i].blen = blen;
	thread_data_array[i].treshold = treshold;
	thread_data_array[i].runs = runs;
	thread_data_array[i].drosophila = drosophila;
	thread_data_array[i].mds = mds;
	thread_data_array[i].scores = scores;
	thread_data_array[i].p = p;
	
	if ((rtn = pthread_create(&threads[i], NULL, (void*)mycompute, (void*)&thread_data_array[i]))) {
	    printf("Error: pthread_create %d, %s\n", i, strerror(rtn));
	    exit(-1);
	}
	
    }
    
    // join the threads
    for (i = 0; i < NUM_THREADS; i++) {
        if ((rtn = pthread_join(threads[i], NULL))) {
	    printf("Error: pthread_join %d, %s\n", i, strerror(rtn));
	    exit(-1);
	}
    }
    
	// destroy mutex variable
    pthread_mutex_destroy(&mutexTASK_ID);
    gettimeofday(&after, NULL);
    printf("regend %d time %g\n", regend, time_ddiff(before, after));
    
}

/*
 * Get the start and stop position of the regionf the given taskid.
 */
inline void get_positions(int wsize, int wstep, int taskid, int *start, int *stop, int regend) {
    int l, r;
    
    /* start value, this ok */
    l = (taskid*TASK_SIZE)*wstep;
    *start = l;
    
    /* end value */
    r = ((taskid+1)*TASK_SIZE)*wstep + (wsize - wstep);
    *stop = r;
}

/**
 * Computing the Cluster Separation Score (CSS) and the corresponding p-value for a given thread.
 * @param threadarg: the arguments to the thread: a thread_data struct
 */
void mycompute(void *threadarg) {
    srand48(time(NULL));  // seed for PRNG
    struct thread_data *my_data;
    
    /* thread vars */
    int tid, my_task_id, wcount = 0;
    
    /* program vars */
    int start, stop, regend, num_windows, wsize, wstep, alen, blen, asize, bsize, treshold, runs, drosophila, mds;
    int wstart, wstop, i, m, dims, npos, idx;
    double result;
    int *apos;
    int *bpos;
    double *avals;
    double *bvals;
    double *scores;
    double *p;
    
    int *aidx;
    int *bidx;
    int *atracks;
    int *btracks;
    int *signtracks;
    double **dissimilarity;
    double **distance;
    double **X;
    double **B;
    double **Z;
    double **tmp;
    double **L;
    double **Q;
    double **r;
    
    // get data from threadarg
    my_data = (struct thread_data *) threadarg;
    tid = my_data->thread_id;
    regend = my_data->regend;
    num_windows = my_data->num_windows;
    wsize = my_data->wsize;
    wstep = my_data->wstep;
    apos = my_data->apos;
    bpos = my_data->bpos;
    avals = my_data->avals;
    bvals = my_data->bvals;
    alen = my_data->alen;
    blen = my_data->blen;
    treshold = my_data->treshold;
    runs = my_data->runs;
    drosophila = my_data->drosophila;
    mds = my_data->mds;
    scores = my_data->scores;
    p = my_data->p;
    
    // initialize array indices
    aidx = (int*)malloc(2*sizeof(int));
    bidx = (int*)malloc(2*sizeof(int));
    
    // get the size of each population
    asize = get_population_size(apos);
    bsize = get_population_size(bpos);
    
    m = asize+bsize;
    dims = 2;
    
    // initialize the arrays with individual ids
    // 0 -> asize-1 for group a,
    // asize -> (asize+bsize) for group b
    atracks = (int*)malloc(asize*sizeof(int));
    btracks = (int*)malloc(bsize*sizeof(int));
    signtracks = (int*)malloc(m*sizeof(int));
    for (i = m; i--;) {
      if (i < asize) atracks[i] = i;
      if (i >= asize) btracks[i-asize] = i;
      signtracks[i] = i;
    }
    
    /* allocate matrixes (2d arrays) */
    allocate_matrix(&dissimilarity, m, m);
    allocate_matrix(&distance, m, m);
    allocate_matrix(&X, m, dims);
    allocate_matrix(&B, m, m);
    allocate_matrix(&Z, m, m);
    allocate_matrix(&tmp, m, m);
    allocate_matrix(&Q, m, dims);
    allocate_matrix(&L, dims, dims);
    allocate_matrix(&r, m, dims);
    
    idx = 0;
    my_task_id = 0;
    
    /* each thread has it's own state - they will generate different streams of pseudo-random numbers */
    unsigned short state[3] = {0,0,0};
    unsigned short seed = time(NULL) + (unsigned short) pthread_self();
    memcpy(state, &seed, sizeof(seed));
    
	// fetch new tasks and calculate css
    while (my_task_id < num_tasks) {
      
        // try to fetch a new task
        pthread_mutex_lock (&mutexTASK_ID);
	my_task_id= task_id;
	task_id++;
	pthread_mutex_unlock(&mutexTASK_ID);
	
	// if no more tasks, exit
	if (my_task_id > num_tasks)
	    break;
	
	// calculate chromosome positions for the given task
	get_positions(wsize, wstep, my_task_id, &start, &stop, regend);
	
	// to include the last (possible smaller) window
	if (stop >= regend) {
	    stop = regend + wstep;
	}

	wstart = start;
	wstop = start+wsize;
	
	aidx[0] = 0; aidx[1] = 0;
	bidx[0] = 0; bidx[1] = 0;
	
	// calculate css for each window
	while (wstart + wsize <= stop) {
	  
	    slide_right(aidx, apos, wstart, wstop, alen);
	    slide_right(bidx, bpos, wstart, wstop, blen);
	    
	    npos = (aidx[1] - aidx[0])/asize;
	    
	    if (npos > 0) {
	      
	        idx = wstart/wstep;
		
		result = cluster_separation_scorer(distance, &(avals[aidx[0]]), &(bvals[bidx[0]]), atracks, btracks, asize, bsize, npos, dissimilarity,
						   drosophila, mds, X, B, Z, tmp, L, Q, r);
		
		if (result != -1) {
		    scores[idx] = result;
		    p[idx] = significance_treshold(distance, signtracks, asize, bsize, result, treshold, runs,state);
		}
	    }
	    wstart += wstep;
	    wstop += wstep;
	    wcount++;
	}
    }

    /* cleanup */
    free(aidx);
    free(bidx);
    free(atracks);
    free(btracks);
    free(signtracks);
    deallocate_matrix(dissimilarity);
    deallocate_matrix(distance);
    deallocate_matrix(X);
    deallocate_matrix(B);
    deallocate_matrix(Z);
    deallocate_matrix(tmp);
    deallocate_matrix(Q);
    deallocate_matrix(L);
    deallocate_matrix(r);
}


