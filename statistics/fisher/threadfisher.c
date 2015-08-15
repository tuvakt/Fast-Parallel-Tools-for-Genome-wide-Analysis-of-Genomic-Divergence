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
#include "threadfisher.h"
#include "cFisher.h"
#include "comparative.h"

#define NUM_THREADS 64      // number of threads
#define TASK_SIZE 100       // number of windows per task

unsigned int task_id;       // global task id
unsigned int num_tasks;     // number of tasks
pthread_mutex_t mutexTASK_ID;    // mutex to protect 'task_id'
struct thread_data thread_data_array[NUM_THREADS];    // array with data for threads

/**
 * Compute the Fisher Exact Test (FET) described by Burke et. al. for the entire chromosome with threads
 * Burke et al: http://www.nature.com/nature/journal/v467/n7315/abs/nature09352.html
 *
 * @param avals a double array with values for group a, ordered by position and then individual id
 * @param bvals a double array with values for group b, ordered by position and then individual id
 * @param apos a double array with positions for group a, ordered by position and then individual id
 * @param bpos a double array with positions for group b, ordered by position and then inidivudal id
 * @param regstart an integer region start
 * @param regend an integer, region end
 * @param wsize the size of the sliding window
 * @param wstep the step size of the sliding windows
 * @param asize the number of individuals in group a
 * @param bsize the number of individuals in group b
 * @param alen the length of the array avals and apos
 * @param blen the length of the array bvals and bpos
 * @param perc the percentile of the FET score to represent the window
 * @param scores an double array with the resulting FET scores for each window
 * @param stddev an double array with the resulting standard deviations for each window
 */
void threadcompute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen,
		   double perc, double *scores, double *stddev) {
	
    struct timeval before, after;
    gettimeofday(&before, NULL);
    pthread_t threads[NUM_THREADS];
    
    int num_windows, rtn;
    num_windows = (regend/wstep)-3;
    printf("num windows %d\n", num_windows);
    //printf("windows per thread %d\n", windows_per_thread);
    num_tasks = num_windows/TASK_SIZE;
    printf("num tasks %d\n", num_tasks);
    
    pthread_mutex_init(&mutexTASK_ID, NULL);    // init mutex
    int i;
    task_id = 0;
    
    // create and start threads
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
	thread_data_array[i].perc = perc;
	thread_data_array[i].scores = scores;
	thread_data_array[i].stddev = stddev;
	
	if ((rtn = pthread_create(&threads[i], NULL, (void*)mycompute, (void*)&thread_data_array[i]))) {
	    printf("Error: pthread_create %d, %s\n", i, strerror(rtn));
	    exit(-1);
	}
    }
    
    // join threads
    for (i = 0; i < NUM_THREADS; i++) {
        if ((rtn = pthread_join(threads[i], NULL))) {
	    printf("Error: pthread_join %d, %s\n", i, strerror(rtn));
	    exit(-1);
	}
    }
    
    // destroy mutex
    pthread_mutex_destroy(&mutexTASK_ID);
    gettimeofday(&after, NULL);
    printf("regend %d time %g\n", regend, time_ddiff(before, after));
}


/*
 * Get the start and stop position of the regionf the given taskid.
 */
inline void get_positions(int wsize, int wstep, int tid, int *start, int *stop, int regend) {
    int l, r;
    
    /* start value, this ok */
    l = (tid*TASK_SIZE)*wstep;
    *start = l;
    
    /* end value */
    r = ((tid+1)*TASK_SIZE)*wstep + (wsize - wstep);
    *stop = r;
}

/**
 * Computing the Fisher's Exact Test (FET) and the corresponding standard deviation for a given thread.
 * @param threadarg: the arguments to the thread: a thread_data struct
 */
void mycompute(void *threadarg) {
    struct thread_data *my_data;
    int tid;
    int start, stop, regend, num_windows, wsize, wstep, alen, blen, asize, bsize, totalpos, npos, nsamples, my_task_id = 0;
    double perc;
    int *apos;
    int *bpos;
    double *avals;
    double *bvals;
    double *scores;
    double *stddev;
    
    /* fetch data from threadarg */
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
    perc = my_data->perc;
    scores = my_data->scores;
    stddev = my_data->stddev;
    
    int wcount = 0;
    int wstart;
    int wstop;
    
    int *aidx;
    int *bidx;
    int *f;
    int *tmp;
    double *results;
    double *fetscores = NULL;
    double *samples = NULL;
    double *stdsamples;
    
    /* get size of population a and b*/
    asize = get_population_size(apos);
    bsize = get_population_size(bpos);
    
    wcount = 0;
    start = 0;
    stop = wsize;
    nsamples = 100;
    
    /* setup idx arrays */
    aidx = (int*)malloc(2*sizeof(int));
    bidx = (int*)malloc(2*sizeof(int));
    
    /* setup arrays for calc. of FET */
    f = (int*)malloc(4*sizeof(int));
    tmp = (int*)malloc(4*sizeof(int));
    results = (double*)malloc(2*sizeof(double));
    stdsamples = (double*)malloc(nsamples*sizeof(double));
    
    int idx = 0;
    
    /* generate seed for the PRNG */
    unsigned short state[3] = {0,0,0};
    unsigned short seed = time(NULL) + (unsigned short) pthread_self();
    memcpy(state, &seed, sizeof(seed));
    
    // fetch new tasks and calculate FET
    while (my_task_id < num_tasks) {
      
        // try to fetch new task
        pthread_mutex_lock (&mutexTASK_ID);
	my_task_id= task_id;
	task_id++;
	pthread_mutex_unlock(&mutexTASK_ID);
	
	// if no more tasks, break
	if (my_task_id > num_tasks)
	    break;
	
	// get chromosome position of the current task
	get_positions(wsize, wstep, my_task_id, &start, &stop, regend);
	
	// to include the last (possibly smaller) window
	if (stop >= regend) {
	    stop = regend + wstep;
	}
	
	wstart = start;
	wstop = start+wsize;
	
	aidx[0] = 0; aidx[1] = 0;
	bidx[0] = 0; bidx[1] = 0;
	
	// calculate FET for each window
	while (wstart + wsize <= stop) {
	  
	    slide_right(aidx, apos, wstart, wstop, alen);
	    slide_right(bidx, bpos, wstart, wstop, blen);
	    
	    npos = (aidx[1] - aidx[0])/asize;
	    
	    if (npos > 0) {
	      
	        fetscores = (double*)realloc(fetscores, npos*sizeof(double));
		samples = (double*)realloc(samples, npos*sizeof(double));
		idx = wstart/wstep;
		
		fisher_exact_test(results, &(avals[aidx[0]]), &(bvals[bidx[0]]), asize, bsize, npos, f, tmp, samples, stdsamples, nsamples, fetscores, state, perc);
		scores[idx] = results[0];
		stddev[idx] = results[1]; 
	    }
	    wstart += wstep;
	    wstop += wstep;
	    wcount++;
	}
    }
    
    
    /* cleanup */
    free(aidx);
    free(bidx);
    free(f);
    free(tmp);
    free(results);
    free(samples);
    free(stdsamples);
    free(fetscores);
}

