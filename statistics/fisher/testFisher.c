#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include "cFisher.h"
#include "comparative.h"
#include "threadfisher.h"

#define FALSE 0
#define TRUE 1

/* Tests 
   To run the test the of the complete C program, you need the following inputs
   - infile1: GTrack file for population 1
   - infile2: GTrack file for population 2
   - outfile: the file to write to
   - file with correct answers

   Notice that the C code can only be run for one chromosome at a time. 
*/

/* compare the serial and to the parallel results*/
void test_thread(int argc, char *argv[]) {

    int i, serial, threads, scount, tcount, length, chunk, ant;
    char buf[32];
    char *pos_end;
    int *spos;
    int *tpos;
    double *sresults;
    double *tresults;
    
    if (argc != 5) {
        printf("usage: %s <infile1> <infile2> <outfile> <file with correct answers>\n", argv[0]);
	exit(-1);
    }
    
    serial = open(argv[4], O_RDONLY);
    if (serial < 0) {
        printf ("error %d opening %s: %s\n", errno, argv[4], strerror(errno));
	exit(-1);
    }
    
    threads = open(argv[3], O_RDONLY);
    if (threads < 0) {
        printf ("error %d opening %s: %s\n", errno, argv[3], strerror(errno));
	exit(-1);
    }
    
    length = 1000000;
    chunk = length;
    spos = (int*)malloc(length*sizeof(int));
    tpos = (int*)malloc(length*sizeof(int));
    sresults = (double*)malloc(length*sizeof(double));
    tresults = (double*)malloc(length*sizeof(double));
    
    printf("reading serial results\n");
    scount = 0;
    length = chunk;
    while (TRUE) {
        ant = readln(serial, buf, 32);
	if (ant <= 0) {
	    break;
	}
	
	/* if the array is too small, realloc
	 * increase with 'chunk' */
	if (scount >= length) {
	    length += chunk;
	    spos = (int*)realloc(spos, length*sizeof(int));
	    sresults = (double*)realloc(sresults, length*sizeof(double));
	}
	
	spos[scount] = (int)strtol(&buf[0], &pos_end, 10);
	sresults[scount] = strtod(pos_end, &pos_end);
	
	scount++;
    }
    
    printf("reading parallel results\n");
    tcount = 0;
    length = chunk;
    while (TRUE) {
        ant = readln(threads, buf, 32);
	if (ant <= 0) {
	  break;
	}
	
	/* if the array is too small, realloc
	 * increase with 'chunk' */
	if (tcount >= length) {
	    length += chunk;
	    tpos = (int*)realloc(tpos, length*sizeof(int));
	    tresults = (double*)realloc(tresults, length*sizeof(double));
	}
	
	tpos[tcount] = (int)strtol(&buf[0], &pos_end, 10);
	tresults[tcount] = strtod(pos_end, &pos_end);
	tcount++;
    }
    
    if (scount != tcount) {
        printf("ERROR: scount was %d but tcount was %d\n", scount, tcount);
	exit(-1);
    }
    
    int count = 0;
    for (i = 0; i < scount; i++) {
        if (spos[i] != tpos[i]) {
	    printf("s: %d t: %d\n", spos[i], tpos[i]);
	    count++;
	}
	
	if (sresults[i]  != tresults[i]) {
	    printf("s: %g t: %g\n", sresults[i], tresults[i]);
	    count++;
	}
    }
    
    if (count == 0) {
        printf("SUCCESS, same result\n");
    } else {
        printf("ERROR, not same result\n");
    }
    
    close(serial);
    close(threads);
    free(sresults);
    free(tresults);
    free(spos);
    free(tpos);
}

void test_all(int argc, char *argv[]) {
  
    /* SNP position (5)
     * val (9)
     * ind (12)
     */
  
    int i, acount, bcount, ant, fileA, fileB;
    FILE *outfile;
    char buf[32];
    char *pos_end;
    int chunk = 1000000;
    int length = chunk;
    
    int *apos = (int*)malloc(length*sizeof(int));
    double *avals = (double*)malloc(length*sizeof(double));
    int *bpos = (int*)malloc(length*sizeof(int));
    double *bvals = (double*)malloc(length*sizeof(double));
    
    double *scores;
    double *stddev;
    int totalpos;
    
    int asize = 11; // cheating
    int bsize = 10; // cheating
    
    int wsize = 2500;
    int wstep = 500;
    
    /* how many positions should we read ? */
    int npos = 100;
    
    if (argc < 4) {
        printf("Usage: %s <file group A> <file group B> <outfile> \n", argv[0]);
	exit(-1);
    }
    
    fileA = open(argv[1], O_RDONLY);
    if (fileA < 0) {
        printf ("error %d opening %s: %s\n", errno, argv[1], strerror(errno));
	exit(-1);
    }
    
    fileB = open(argv[2], O_RDONLY);
    if (fileB < 0) {
        printf ("error %d opening %s: %s\n", errno, argv[2], strerror(errno));
	exit(-1);
    }
    
    outfile = fopen(argv[3], "w+");
    if (outfile < 0) {
        printf("error %d openind %s: %s\n", errno, argv[3], strerror(errno));
    }
    
    /* reading in file A */
    for (i = 0; i < 5; i++) {
        ant = readln(fileA, buf, 32);
	
	if (ant <= 0) {
	  break;  //END OF FILE
	}
    }
    
    acount = 0;
    while(TRUE) {
        ant = readln(fileA, buf, 32);
	if (ant <= 0) {
	  break;
	}
	
	char c = buf[0];
	i = 0;
	while (c != '\t') {
	    c = buf[++i];
	}
	
	/* if the array is too small, realloc
	 * increase with 'chunk' */
	if (acount >= length) {
	    length += chunk;
	    apos = (int*)realloc(apos, length*sizeof(int));
	    avals = (double*)realloc(avals, length*sizeof(double));
	}
	
	apos[acount] = (int)strtol(&buf[i], &pos_end, 10);  // atoi needs a pointer
	avals[acount] = strtod(pos_end, &pos_end);
	
	//printf("position %d val %g \n", apos[acount], avals[acount]);
	acount++;
    }
    
    /* read in file B */
    for (i = 0; i < 5; i++) {
        ant = readln(fileB, buf, 32);
	
	if (ant <= 0) {
	    break;  //END OF FILE
	}
    }
    
    bcount = 0;
    length = chunk;
    //while (bcount < npos*bsize) {
    while(TRUE) {
        ant = readln(fileB, buf, 32);
	if (ant <= 0) {
	    break;
	}
	
	/* if the array is too small, realloc
	 * increase with 'chunk' */
	if (bcount >= length) {
	    length += chunk;
	    bpos = (int*)realloc(bpos, length*sizeof(int));
	    bvals = (double*)realloc(bvals, length*sizeof(double));
	}
	
	char c = buf[0];
	i = 0;
	while (c != '\t') {
	     c = buf[++i];
	}
	
	bpos[bcount] = (int)strtol(&buf[i], &pos_end, 10);  // atoi needs a pointer
	bvals[bcount] = strtod(pos_end, &pos_end);

	bcount++;
    }
    
    
    // TODO realloc() to 'count' size
    apos = (int*)realloc(apos, acount*sizeof(int));
    avals = (double*)realloc(avals, acount*sizeof(double));
    bpos = (int*)realloc(bpos, bcount*sizeof(int));
    bvals = (double*)realloc(bvals, bcount*sizeof(double));
    
    close(fileA);
    close(fileB);
    
    totalpos = (apos[acount-1]+1)/wstep;
    scores = (double*)malloc(totalpos*sizeof(double));
    stddev = (double*)malloc(totalpos*sizeof(double));
    
    memset(scores, 0, sizeof(double)*totalpos);
    memset(stddev, 0, sizeof(double)*totalpos);
    
    // testing.. scary!
    
    double perc = 0.95;
    struct timeval before, after;
    
    gettimeofday(&before, NULL);
    /* for serial program: call method 'compute'. */
    threadcompute(avals, bvals, apos, bpos, 0, (apos[acount-1]+1), wsize, wstep, acount, bcount, perc, scores, stddev);
    gettimeofday(&after, NULL);
    
    for (i = 0; i < totalpos; i++) {
        if (scores[i] != 0) {
	  // printf("pos %d score %f std %f\n", i*wstep, scores[i], stddev[i]);
	  fprintf(outfile, "%d\t%g\t%g\n", i*wstep, scores[i], stddev[i]);
	  //printf("std[%d]: %g\n", i, stddev[i]);
	}
    }   
    fclose(outfile);
    
    printf("last position %d \n", apos[acount-1]);
    printf("Time in seconds %f \n", time_ddiff(before, after));
    
    free(apos);
    free(avals);
    free(bpos);
    free(bvals);
    free(scores);
    free(stddev);
}


void test_compute_and_slide_right() {
    int asize, bsize, regstart, regend, wsize, wstep, alen, blen, npos;
    
    asize = 3;
    bsize = 3;
    regstart = 0;
    regend = 10;
    wsize = 5;
    wstep = 2;
    
    npos = 6;
    
    double perc = 0.95;
    
    alen = asize*npos;
    blen = bsize*npos;
    
    double avals[18] = {-3,3,-3,3,-3,2,2,2,2,2,-3,-3,-3,-3,-3,3,3,3};
    double bvals[18] = {2,2,2,2,2,-3,3,-3,3,-3,3,3,3,3,3,-3,-3,-3};
    
    int apos[18] = {2,2,2,4,4,4,5,5,5,6,6,6,8,8,8,10,10,10};
    int bpos[18] = {2,2,2,4,4,4,5,5,5,6,6,6,8,8,8,10,10,10};
    
    double *scores = (double*)malloc(18*sizeof(double));
    double *stddev = (double*)malloc(18*sizeof(double));
    
    compute(avals, bvals, apos, bpos, regstart, regend, wsize, wstep, alen, blen, perc, scores, stddev);
    
}


void test_mean() {
    int i, n;
    double *vals;
    double res;

    n = 10;
    vals = (double*)malloc(n*sizeof(double));

    for (i = 0; i < n; i++) {
    	vals[i] = i+1;
    }

    for (i = 0; i < n; i++) {
    	printf(" %g", i, vals[i]);
	}
    printf("\n");
    res = mean(vals, n);
    printf("MEAN: Expected %g got %g\n", 5.5, res);

    /* clean up */
    free(vals);

}
//
void test_std() {
    int i, n;
    double *vals;
    double res;

    n = 10;
    vals = (double*)malloc(n*sizeof(double));

    for (i = 0; i < n; i++) {
    	vals[i] = i+1;
    }

    res = std(vals, n);
    printf("STD: Expected %g got %g\n", 2.8722813232690143, res);

    /* clean up */
    free(vals);
}

void test_percentile() {
    int n;
    double res;
    double vals[25] = {43, 54, 56, 61, 62, 66, 68, 69, 69, 70, 71, 72, 77, 78, 79, 85, 87, 88, 89, 93, 95, 96, 98, 99, 99};
    res = percentile(vals, 25, 0.90);
    printf("PERCENTILE: Expected %g got %g\n", 97.2, res);

    double vals2[10] = {0,1,2,3,4,5,6,7,8,9};
    res = percentile(vals2, 10, 0.5);
    printf("PERCENTILE: Expected %g got %g\n", 4.5, res);

}

void test_min_idx() {
    int idx;
    int a[10] = {2, 5, 7, 1, 11, 2, 76, 3, 6, 10};
    idx = min_idx((int*)a, 10);
    printf("MIN IDX: Expected 1 got %d\n", a[idx]);

}

void test_binomial() {
    printf("expected %lu got %lu\n", 10, binomial(5, 3));
    printf("expected %lu got %lu\n", 131282408400, binomial(40, 19));
    printf("expected %lu got %lu\n", 11923179284862717872, binomial(67, 31));
}

void test_fetcount() {
    int i, asize, bsize;
    int *f;

    asize = 10;
    bsize = 10;

    f = (int*)malloc(4*sizeof(int));

    double avals[10*2] = {3, 3, 0, -10000, 3, -3, 0, 0, 0, 3, 3, 3, 0, -10000, 3, -3, 0, 0, 0, 3};
    double bvals[10*2] = {-3,-3, 0, -10000, -3, -3, 3, 3, 3, 0, -3,-3, 0, -10000, -3, -3, 3, 3, 3, 0};
    printf("All ok so far\n");

    printf("group a: ");
    for (i = 0; i < asize*2; i++) {
    	printf("%g ", avals[i]);
    }
    printf("\n");

    printf("group b: ");
    for (i = 0; i < bsize*2; i++) {
    	printf("%g ", bvals[i]);
    }
    printf("\n");

    fetcount(f, avals, bvals, 0, asize, bsize);

    printf("frequencies:\nSHOULD BE:\n4 | 1\n3 | 4\n");
    printf("IS:\n%d | %d\n%d | %d\n", f[0], f[1], f[2], f[3]);

    free(f);

}

void test_shift_table() {

    int n = 4;
    int *f = (int*)malloc(n*sizeof(int));
    int *tmp = (int*)malloc(n*sizeof(int));
    
    f[0] = 3; f[1] = 5; f[2] = 2; f[3] = 7;
    
    printf("F:\n%d | %d\n%d | %d\n", f[0], f[1], f[2], f[3]);
    shift_table(f, tmp, n);
    printf("F SHOULD BE:\n2 | 3\n7 | 5\n");
    printf("F IS:\n%d | %d\n%d | %d\n", f[0], f[1], f[2], f[3]);
    
    free(f);
    free(tmp);
}

void test_fet_p() {
    printf("Expected %g got %g\n", 0.001346076, fet_p(1, 9, 11, 3));
    printf("Expected %g got %g\n", 0.000033652, fet_p(0, 10, 12, 2));
    printf("Expected %g got %g\n", 0.0166, fet_p(9, 3, 1, 6));
    printf("\n");
}

void test_fet() {
    double res;
    int *tmp;

    tmp = (int*)malloc(4*sizeof(int));

    int a[4] = {2, 7, 8, 2};
    res = fet((int*)a, tmp);
    printf("TEST FET: 2782 Expected %g got %g\n", 0.0230141, res);

    int b[4] = {2, 3, 6, 4};
    res = fet((int*)b, tmp);
    printf("TEST FET: 2364 Expected %g got %g\n", 0.6083916, res);

    //double size, goes over 1
    int c[4] = {2, 2, 3, 3};
    res = fet((int*)c, tmp);
    printf("TEST FET: 2233 Expected %d got %g\n", 1, res);

    // goes over 1
    int d[4] = {1, 3, 2, 3};
    res = fet((int*)d, tmp);
    printf("TEST FET: 1323 Expected %d got %g\n", 1, res);

    free(tmp);

}



/* tests for inlined functions. */
void test_boostrap_sample() {

    double fetscores[5] = {1, 0.9, 0.6, 0.5, 0.3};
    double *samples;
    int i, n, npos;
    
    unsigned short state[3] = {0,0,0};
    unsigned short seed = time(NULL);
    printf("SEED %d \n", seed);
    memcpy(state, &seed, sizeof(seed));
    
    n = 10;
    npos = 5;
    
    samples = (double*)malloc(n*sizeof(double));
    
    printf("FETSCORES: ");
    for (i = 0; i < npos; i++) {
        printf("%g ", fetscores[i]);
    }
    printf("\n");
    bootstrap_sample(fetscores, samples, n, npos, state);
    printf("SAMPLES 1: ");
    for (i = 0; i < n; i++) {
        printf("%g ", samples[i]);
    }
    printf("\n");
    
    bootstrap_sample(fetscores, samples, n, npos, state);
    printf("SAMPLES 2: ");
    for (i = 0; i < n; i++) {
        printf("%g ", samples[i]);
    }
    printf("\n");
    
    bootstrap_sample(fetscores, samples, n, npos, state);
    printf("SAMPLES 3: ");
    for (i = 0; i < n; i++) {
        printf("%g ", samples[i]);
    }
    printf("\n");
    
    free(samples);
}

/* inlined functions */
/*void test_sort() {
  
    double fetscores[10] = {0.3, 0.6, 0.1, 0.2, 0.8, 0.3, 0.2, 0.5, 1, 0.01};
    int i, n = 10;
    qsort((double*)fetscores, n, sizeof(double), compare_doubles);
    
    printf("fetscores IS:\n");
    for ( i= 0; i < n; i++) {
        printf("%g ", fetscores[i]);
    }
    printf("\n");
}*/

void test_swap() {
    int a, b;
    a = 10;
    b = 5;
    printf("Before: a = %d, b = %d\n", a, b);
    swap(&a, &b);
    printf("After: a = %d, b = %d\n", a, b);
}

int main(int argc, char *argv[]) {
 
    /* inlined functions, test commented out */
    //printf("\n*** Testing bootstrap sample ***\n");
    //test_boostrap_sample();
    //printf("\n*** Testing compare doubles and qsort\n");
    //test_sort();
    //printf("\n*** Testing swap ***\n");
    //test_swap();
     
    /* regular functions */
    printf("\n*** Testing percentile ***\n");
    test_percentile();
    printf("\n*** Testing min idx ***\n");
    test_min_idx();
    printf("\n*** Testing mean ***\n");
    test_mean();
    printf("\n*** Testing std dev ***\n");
    test_std();
    printf("\n*** Testing binomial ***\n");
    test_binomial();
    printf("\n*** Testing fetcount ***\n");
    test_fetcount();
    printf("\n*** Testing shift table ***\n");
    test_shift_table();
    printf("\n*** Testing fet p0 ***\n");
    test_fet_p();
    printf("\n*** Testing fet ***\n");
    test_fet();

    test_all(argc, argv);
    test_thread(argc, argv);
    
    return 0;
    
}
