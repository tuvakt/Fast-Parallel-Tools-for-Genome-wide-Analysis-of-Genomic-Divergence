#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "css.h"
#include "comparative.h"
#include "threadcss.h"
#include "css_smacof.h"

#define FALSE 0
#define TRUE 1

void test_thread(int argc, char *argv[]) {

	int i, serial, threads, scount, tcount, length, chunk, ant;
	char buf[32];
	char *pos_end;
	int *spos;
	int *tpos;
	double *sresults;
	double *tresults;

	if (argc != 7) {
		printf("usage: %s 'x' mdsalg <infile1> <infile2> <outfile> <file with correct answers>\n", argv[0]);
		exit(-1);
	}

	serial = open(argv[6], O_RDONLY);
	if (serial < 0) {
		printf ("error %d opening %s: %s\n", errno, argv[1], strerror(errno));
		exit(-1);
	}

	threads = open(argv[5], O_RDONLY);
	if (threads < 0) {
		printf ("error %d opening %s: %s\n", errno, argv[2], strerror(errno));
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
			//printf("increase with %d too %d\n", chunk, length);
		}

		tpos[tcount] = (int)strtol(&buf[0], &pos_end, 10);
		tresults[tcount] = strtod(pos_end, &pos_end);
		//printf("%g\n", tresults[tcount]);
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

	int i, acount, bcount, ant, fileA, fileB, pcount, drosophila, mds;
	FILE *outfile;
	double val;
	char buf[32];
	char *pos_end;
	int chunk = 1000000;
	int length = chunk;

	int *apos = (int*)malloc(length*sizeof(int));
	double *avals = (double*)malloc(length*sizeof(double));
	int *bpos = (int*)malloc(length*sizeof(int));
	double *bvals = (double*)malloc(length*sizeof(double));

	double *scores;
	double *p;
	int totalpos;

	int asize = 11; // cheating
	int bsize = 10; // cheating

	int wsize = 2500;
	int wstep = 500;

/* how many positions should we read ? */
	int npos = 100;

	printf("Testing\n");
//
	if (argc < 6) {
		printf("Usage: %s 'x' mdsalg <file group A> <file group B> <outfile> \n", argv[0]);
		exit(-1);
	}

	// drosophila = 0: stickleback
	// drosophila = 1: drosophila
	drosophila = (int)strtol(argv[1], &pos_end,10);
	printf("drosophila %d\n", drosophila);
	mds = (int)strtol(argv[2], &pos_end, 10);
	printf("mdsalg %d\n", mds);

	printf("cprog %s drosophila %s mdsalg %s infileA %s infileB %s outfile%s\n", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);

	fileA = open(argv[3], O_RDONLY);
	if (fileA < 0) {
		printf ("error %d opening %s: %s\n", errno, argv[2], strerror(errno));
		exit(-1);
	}

	fileB = open(argv[4], O_RDONLY);
	if (fileB < 0) {
		printf ("error %d opening %s: %s\n", errno, argv[3], strerror(errno));
		exit(-1);
	}

	/*outfile = fopen(argv[5], "w+");
	if (outfile < 0) {
		printf("error %d openind %s: %s\n", errno, argv[3], strerror(errno));
		}*/

	printf("reading file A\n");
	/* reading in file A */
	for (i = 0; i < 5; i++) {
		ant = readln(fileA, buf, 32);

		if (ant <= 0) {
			break;  //END OF FILE
		}
	}

	acount = 0;
	//while (acount < npos*asize) {
	while(TRUE) {
		ant = readln(fileA, buf, 32);
		if (ant <= 0) {
			break;
		}

		/* if the array is too small, realloc
		 * increase with 'chunk' */
		if (acount >= length) {
			length += chunk;
			apos = (int*)realloc(apos, length*sizeof(int));
			avals = (double*)realloc(avals, length*sizeof(double));

			//printf("increase with %d too %d\n", chunk, length);
		}

		char c = buf[0];
		i = 0;
		while (c != '\t') {
			c = buf[++i];
		}
		apos[acount] = (int)strtol(&buf[i], &pos_end, 10);  // atoi needs a pointer
		val = strtod(pos_end, &pos_end);
		avals[acount] = val;

		//printf("position %d val %g \n", apos[acount], avals[acount]);
		acount++;
	}

	printf("reading file B\n");
	/* read in file B */
	for (i = 0; i < 5; i++) {
		ant = readln(fileB, buf, 32);

		if (ant <= 0) {
			break;  //END OF FILE
		}
	}

	bcount = 0;
	length = chunk;
	//while (bcount < 10) {
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

			//printf("increase with %d too %d\n", chunk, length);
		}

		char c = buf[0];
		i = 0;
		while (c != '\t') {
			c = buf[++i];
		}

		bpos[bcount] = (int)strtol(&buf[i], &pos_end, 10);  // atoi needs a pointer

		val = strtod(pos_end, &pos_end);
		bvals[bcount] = val;

		//printf("position %d val %g \n", bpos[bcount], bvals[bcount]);
		bcount++;
	}

	// TODO realloc() to 'count' size
	apos = (int*)realloc(apos, acount*sizeof(int));
	avals = (double*)realloc(avals, acount*sizeof(double));
	bpos = (int*)realloc(bpos, bcount*sizeof(int));
	bvals = (double*)realloc(bvals, bcount*sizeof(double));
	printf("Allocated apos %d + avals %d \n", (acount*sizeof(int)), (acount*sizeof(double)));
	printf("Allocated bpos %d + bvals %d \n", (bcount*sizeof(int)), (bcount*sizeof(double)));
	
	close(fileA);
	close(fileB);

	totalpos = (apos[acount-1]+1)/wstep;
	scores = (double*)malloc(totalpos*sizeof(double));
	p = (double*)malloc(totalpos*sizeof(double));
	printf("Allocated scores and p: %d\n", (2*totalpos*sizeof(double)));

	printf("totalpos %d\n", totalpos);
	memset(scores, 0, sizeof(double)*totalpos);
	memset(p, 0, sizeof(double)*totalpos);

	// testing.. scary!

	struct timeval before, after;

	int treshold = 10;
	int runs = 200000;
	printf("into the fire\n");
	gettimeofday(&before, NULL);
	compute(avals, bvals, apos, bpos, 0, (apos[acount-1]+1), wsize, wstep, acount, bcount, treshold, runs, drosophila, mds, scores, p);
	gettimeofday(&after, NULL);
	printf("Time in seconds %f \n", time_ddiff(before, after));
	/*for (i = 0; i < totalpos; i++) {
		  if (scores[i] != 0) {
			  //printf("pos %d score %f p %f\n", i*wstep, scores[i], p[i]);
			  fprintf(outfile, "%d\t%g\t%g\n", i*wstep, scores[i], p[i]);
		  }
	}

	fclose(outfile);*/
	

	printf("last position %d \n", apos[acount-1]);

	free(apos);
	free(avals);
	free(bpos);
	free(bvals);
	free(scores);
	free(p);
}

//void test_get_positions() {
//
//	int wsize = 2500;
//	int wstep = 500;
//	int num_threads = 3;  // 4;
//	int regstart = 0;
//	int regend = 7000;
//	int num_windows = ((regend/wstep)-4);
//	//int windows_per_thread = ((regend/wstep)-4)/num_threads;
//	printf("num_windows %d num_threads %d\n", num_windows, num_threads);
//
//	int start[num_threads];
//	int stop[num_threads];
//	int tid;
//	for (tid = 0; tid < num_threads; tid++) {
//		get_positions(wsize, wstep, (long)tid, &start[tid], &stop[tid], regend);
//	}
//
//	for (tid = 0; tid < num_threads; tid++) {
//		printf("tid %d:\n", tid);
//		int nw = ((stop[tid]-start[tid])/wstep)-4;
//		printf("start %d stop %d number of windows: %d\n", start[tid], stop[tid], nw);
//	}
//}

void test_compare_freq() {
	int i, asize, bsize, npos, m;
	double *avals;
	double *bvals;
	double **compared;

	asize = 1;
	bsize = 1;
	m = asize + bsize;
	npos = 3;

	avals = (double*)malloc(asize*npos*sizeof(double));
	bvals = (double*)malloc(bsize*npos*sizeof(double));

	avals[0] = 0.75; avals[1] = 0.3; avals[2] = 0.2;
	bvals[0] = 0.2; bvals[1] = 0.5; bvals[2] = 0.75;

	allocate_matrix(&compared, m, m);
	for (i = 0; i < m; i++) {
		memset(compared[i], 0, m*sizeof(double));
	}

	compare_freq(avals, bvals, npos, compared);

	printf("SHOULD BE\n0\t0.433333\n0.433333\t0\n");
	printf("IS:\n%g\t%g\n%g\t%g\n", compared[0][0], compared[0][1], compared[1][0], compared[1][1]);

	free(avals);
	free(bvals);
	deallocate_matrix(compared);
}

void test_compare_all() {
    int i, j, asize, bsize, npos, m;
    double **d;
    double *avals;
    double *bvals;

    asize = 2;
    bsize = 2;
    npos = 2;

    avals = (double*)malloc(asize*npos*sizeof(double));
    bvals = (double*)malloc(bsize*npos*sizeof(double));

    avals[0] = -3; avals[1] = 3; avals[2] = 0; avals[3] = 3;
    bvals[0] = 3; bvals[1] = -3; bvals[2] =- 10000; bvals[3] = 0;

    for ( i= 0; i < asize*npos; i++) {
    	printf("%g ", avals[i]);
    }
    printf("\n");

    for ( i= 0; i < bsize*npos; i++) {
    	printf("%g ", bvals[i]);
    }
    printf("\n");

    m = asize + bsize;
    printf("m %d\n", m);

    allocate_matrix(&d, m, m);
    for ( i= 0; i < m; i++) {
    	for (j = 0; j < m; j++) {
    		d[i][j] = 0;
    	}
    }

    compare_all(avals, bvals, asize, bsize, npos, d);

    printf("done!\n");

    printf("*** d *** \n");
    for (i = 0; i < m; i++) {
    	for (j = 0; j < m; j++) {
    		printf("%g ", d[i][j]);
    	}
    	printf("\n");
    }


    /* clean up */
    free(avals);
    free(bvals);
    deallocate_matrix(d);
}

void test_fill_averages() {

	int i, j, m;
	double **compared;

	m = 5;
	allocate_matrix(&compared, m, m);

	for (i = 0; i < m; i++) {
		compared[i][i] = 0;
		for (j = i+1; j < m; j++) {
			compared[i][j] = i+j;
			compared[j][i] = i+j;
		}
	}

	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g\t", compared[i][j]);
		}
		printf("\n");
	}


	int res = fill_averages(compared, m);
	printf("Result: SHOULD BE %d is %d\n", 1, res);
	printf("compared[0][0]: SHOULD BE 3.2 is %g\n", compared[0][0]);
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", compared[i][j]);
		}
		printf("\n");
	}

	for (i = 0; i < m; i++) {
		compared[i][i] = 0;
		compared[i][0] = 0;
		compared[i][m-1] = 0;
		compared[0][i] = 0;
	}

	res = fill_averages(compared, m);
	printf("Result: SHOULD BE %d is %d\n", 0, res);
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", compared[i][j]);
		}
		printf("\n");
	}

	deallocate_matrix(compared);
}

void test_calc_dist() {

    int i, j, m, dims, count;
    double  **A;
    double **distance;

    m = 2;
    dims = 2;

    allocate_matrix(&A, m, dims);
    allocate_matrix(&distance, m, m);

    // init A
    for (i = 0; i < m; i++) {
    	for (j = 0; j < dims; j++) {
    		A[i][j] = i+j;
    		printf("%g ", A[i][j]);
    	}
    	printf("\n");
    }

    double expected[4] = {0, 1.41421356, 1.41421356, 0};
    count = 0;
    calc_dist(A, distance, m);
    for (i = 0; i < m; i++) {
    	for (j = 0; j < m; j++) {
    		printf("expected: %g distance %g\n", expected[count], distance[i][j]);
    		count++;
    	}
    }

    /* clean up */
    deallocate_matrix(A);
    deallocate_matrix(distance);

}


void test_matrix_mult() {

	int i, j, m, n, p;

	double **A;
	double **B;
	double **C;

	m = 4;
	n = 3;
	p = 2;

	allocate_matrix(&A, m, n);
	allocate_matrix(&B, n, p);
	allocate_matrix(&C, m, p);

	double a[4*3] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

	printf("A:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			A[i][j] = a[n * i + j];
			printf(" %g", A[i][j]);
	    }
		printf("\n");
	}


	printf("B:\n");
	double b[3*2] = {6, 5, 4, 3, 2, 1};
	for (i = 0; i < n; i++) {
		for (j = 0; j < p; j++) {
			B[i][j] = b[p * i + j];
			printf(" %g", B[i][j]);
	    }
		printf("\n");
	}

	printf("Multiplying a %dx%d matrix with a %dx%d matrix: result should be %dx%d\n", m, n, n, p, m, p);
	matrix_mult(A, B, C, m, n, p);

	printf("m %d p %d\n",m, p);
	printf("C should be:\n20 14\n56 41\n92 68\n128 95\n");

	printf("C is:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < p; j++) {
			printf("%g ", C[i][j]);
		}
		printf("\n");
	}

	deallocate_matrix(A);
	deallocate_matrix(B);
	deallocate_matrix(C);
}

void test_cmds() {
	int i, j, n, m, dims;
	double **compared;
	double **X;
	double **B;
	double **Z;
	double **tmp;
	double **L;
	double **Q;

	n = 4;
	dims = 2;
	m = n;

	allocate_matrix(&compared, m, m);
	allocate_matrix(&X, m, dims);
	allocate_matrix(&B, m, m);
	allocate_matrix(&Z, m, m);
	allocate_matrix(&tmp, m, m);
	allocate_matrix(&Q, m, dims);
	allocate_matrix(&L, dims, dims);

    /* Test 1 */
    double a[4*4] = {0, 4.05, 8.25, 5.57, 4.05, 0, 2.54, 2.69, 8.25, 2.54, 0, 2.11, 5.57, 2.69, 2.11, 0};

    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
    		compared[i][j] = a[n * i + j];
    	}
	}

	printf("compared:\n");
    for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
    		printf("%g ", compared[i][j]);
    	}
    	printf("\n");
    }
    printf("\n");

    cmds(compared, X, dims,  m,B, Z,tmp, L, Q);

    printf("X SHOULD BE:\n");
    printf("4.62 0.07\n0.09 -1.11\n-3.63 -0.34\n-1.08 1.38\n\n");
    printf("X IS:\n");
    for (i = 0; i < n; i++) {
    	for (j = 0; j < dims; j++) {
    		printf("%g ", X[i][j]);
    	}
		printf("\n");
    }

    /* cleanup */
    deallocate_matrix(compared);
    deallocate_matrix(X);
    deallocate_matrix(B);
    deallocate_matrix(Z);
    deallocate_matrix(tmp);
    deallocate_matrix(Q);
    deallocate_matrix(L);

}

void test_smacof() {
	int i, j, m = 4, dims = 2, max_iters = 300, n_init = 1;
	double epsilon = 0.000001; // 1e-6
	double **compared;
	double **X;
	double **Z;
	double **D;
	double **B;
	double **result;

	allocate_matrix(&compared, m, m);
	allocate_matrix(&X, m, dims);
	allocate_matrix(&Z, m, dims);
	allocate_matrix(&D, m, m);
	allocate_matrix(&B, m, m);
	allocate_matrix(&result, m, dims);


	double a[4*4] = {0,5,3,4,5,0,2,2,3,2,0,1,4,2,1,0};

	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			compared[i][j] = a[m * i + j];
	   	}
	}

	printf("compared:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", compared[i][j]);
	    }
	    printf("\n");
	}

    double b[4*2] = {-0.266, -0.539, 0.451, 0.252, 0.016, -0.238, -0.200, 0.524};
    for (i = 0; i < m; i++) {
    	for (j = 0; j < dims; j++) {
    		X[i][j] = b[dims * i + j];
    	}
    }

    printf("X:\n");
    for (i = 0; i < m; i++) {
    	for (j = 0; j < dims; j++) {
    		printf("%g ", X[i][j]);
    	}
    	printf("\n");
    }

    smacof(compared, m, dims, X, Z, D, B, max_iters, epsilon);

    printf("X SHOULD BE:\n");
    printf("-1.457 -2.575\n1.730 1.23\n-0.028 0.16\n-0.245 1.185\n\n");
    printf("X IS:\n");
    for (i = 0; i < m; i++) {
    	for (j = 0; j < dims; j++) {
    		printf("%g ", X[i][j]);
        }
    	printf("\n");
    }

//    smacof_runs(compared, m, dims, X, Z, D, B, result, max_iters, n_init, epsilon);
//    printf("X with X0 RANDOM IS:\n");
//    for (i = 0; i < m; i++) {
//    	for (j = 0; j < dims; j++) {
//    		printf("%g ", X[i][j]);
//        }
//       	printf("\n");
//    }

	deallocate_matrix(compared);
	deallocate_matrix(X);
	deallocate_matrix(Z);
	deallocate_matrix(B);
	deallocate_matrix(D);
	deallocate_matrix(result);
}

void test_css() {
    int i, j, m, dims, a, b;

    double result;
    double **X;
    double **distance;

    m = 100;
    dims = 2;
    allocate_matrix(&distance, m, m);
    allocate_matrix(&X, m, dims);

    // fill A
    for (i = 0; i < m; i++) {
    	for (j = 0; j < dims; j++) {
    		X[i][j] = i+j;
    		//printf(" %g ", X[i][j]);
    	}
    	//printf("\n");
    }

    // testing 'css'
    a = 50;
    b = 50;

    int group_a[a];
    int group_b[b];

    for (i = 0; i < a; i++) {
    	group_a[i] = i;
    }

    for (j = 0; j < b; j++) {
    	group_b[j] = j+a;
    }

    // calculating the distance first
    calc_dist(X, distance, m);
    result = css(distance, group_a, group_b, a, b);
    double expected = 70.5975410337;

    double bools = expected - result;
    if (bools <= 0.00001) {
    	printf("TRUE: ");
    } else {
    	printf("FALSE: ");
    }
    printf("result of css: %f expected: %f\n", result, expected);

    /* clean up */
    deallocate_matrix(distance);
    deallocate_matrix(X);
}


void test_swap() {
    int a, b;
    a = 10;
    b = 5;
    printf("Before: a = %d, b = %d\n", a, b);
    swap(&a, &b);
    printf("After: a = %d, b = %d\n", a, b);
}

void test_copy_matrix() {
	int i, j;
	int m;

	double **Z;
	double **X;

	m = 3;
	allocate_matrix(&X, m, 2);
	allocate_matrix(&Z, m, 2);


	printf("X IS: \n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < 2; j++) {
			X[i][j] = i+j;
			printf("%g ", X[i][j]);
		}
		printf("\n");
	}

	copy_matrix(Z, X, m);

	X[0][0] = 4;

	printf("Z IS: \n");
	for (i = 0; i < m; i++) {
		printf("%g %g\n", Z[i][0], Z[i][1]);
	}

	deallocate_matrix(Z);
	deallocate_matrix(X);
}

void test_stress() {

	int i, j, m;

	double **delta;
	double **D;

	m = 4;

	allocate_matrix(&delta, m, m);
	allocate_matrix(&D, m, m);

	double a[4*4] = {0,5,3,4,5,0,2,2,3,2,0,1,4,2,1,0};
	double b[4*4]= {0, 1.068, 0.412, 1.065, 1.068, 0, 0.655, 0.706, 0.412, 0.655, 0, 0.792, 1.065, 0.706, 0.792, 0};

	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			delta[i][j] = a[m * i + j];
			D[i][j] = b[m * i + j];
		}
	}

	printf("delta:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", delta[i][j]);
		}
		printf("\n");
	}
	printf("D:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", D[i][j]);
		}
		printf("\n");
	}

	double sigma = stress(delta, D, m);
	printf("\nSTRESS SHOULD BE %g IS %g\n", 34.29899413, sigma);

	deallocate_matrix(delta);
	deallocate_matrix(D);
}

void test_random_shuffle() {

	int i, n;
	int *elms;

	n = 5;

	unsigned short state[3] = {0,0,0};
	unsigned short seed = time(NULL);
	printf("SEED %d \n", seed);
	memcpy(state, &seed, sizeof(seed));

	elms = malloc(n*sizeof(double));
	for (i = 0; i < n; i++) {
		elms[i] = i+1;
	}

	printf("ELMS: ");
	for (i = 0; i < n; i++) {
		printf("%d ", elms[i]);
	}
	printf("\n");

	random_shuffle(elms, n, state);
	printf("SHUFFLE 1: ");
	for (i = 0; i < n; i++) {
		printf("%d ", elms[i]);
	}
	printf("\n");

	random_shuffle(elms, n, state);
	printf("SHUFFLE 2: ");
	for (i = 0; i < n; i++) {
		printf("%d ", elms[i]);
	}
	printf("\n");

	random_shuffle(elms, n, state);
	printf("SHUFFLE 3: ");
	for (i = 0; i < n; i++) {
		printf("%d ", elms[i]);
	}
	printf("\n");

	free(elms);

}

void test_setup_z_matrix() {

	int i, j, m;
	double **A;

	m = 4;
	allocate_matrix(&A, m, m);

	setup_z_matrix(A, m);

	printf("Z should be:\n0.75 on the diagonal, -0.25 elsewhere\n");
	printf("Z IS:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", A[i][j]);
		}
		printf("\n");
	}

	deallocate_matrix(A);
}

void test_matrix_sqrt() {

	int i, j, m;
	double **A;

	m = 2;
	allocate_matrix(&A, m, m);

	A[0][0] = 100; A[0][1] = 100; A[1][0] = 100; A[1][1] = 100;

	printf("A SHOULD BE:\n");
	printf("10 10\n10 10\n");

	matrix_sqrt(A, m);

	printf("A IS:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", A[i][j]);
		}
		printf("\n");
	}

	deallocate_matrix(A);
}

void test_matrix_squared() {
	int i, j, m;
	double **A;
	double **B;

	m = 2;
	allocate_matrix(&A, m, m);
	allocate_matrix(&B, m, m);

	A[0][0] = 10; A[0][1] = 10; A[1][0] = 10; A[1][1] = 10;

	printf("B SHOULD BE:\n");
	printf("100 100\n100 100\n");

	matrix_squared(A,B, m);

	printf("B IS:\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) {
			printf("%g ", B[i][j]);
		}
		printf("\n");
	}

	deallocate_matrix(A);
	deallocate_matrix(B);
}

void test_sign_tresh() {

	int i, j, n, m, asize, bsize;

	unsigned short state[3] = {0,0,0};
	unsigned short seed = time(NULL);
	memcpy(state, &seed, sizeof(seed));

	int tracks[10] = {0,1, 2,3,4,5,6,7,8,9};
	int *atracks;
	int *btracks;

	asize = 6;
	bsize = 4;
	m = asize+bsize;
	n = 5;

	printf("tracks: ");
	for (i = 0; i < m; i++) {
		printf("%d ", tracks[i]);
	}
	printf("\n");

	for (j = 0; j < n; j++) {
		random_shuffle(tracks, m, state);
		atracks = &(tracks[0]);
		btracks = &(tracks[asize]);

		printf("atracks: ");
		for (i = 0; i < asize; i++) {
			printf("%d ", atracks[i]);
		}
		printf("\n");

		printf("btracks: ");
		for (i = 0; i < bsize; i++) {
			printf("%d ", btracks[i]);
		}
		printf("\n\n");
	}
}

void test_slide_right() {

	int alen = 20;
	int asize = 4;
	int npos = 0;
	int regstart = 0;
	int regend = 1000;
	int wsize = 400;
	int wstep = 300;

	int start = 0;
	int stop = wsize;
	int pos[20] = {262,262,262,262,300,300,300,300,500,500,500,500,667,667,667,667,999,999,999,999};
	int idx[2] = {0,0} ;

	 while (start + wsize <= regend + wstep) {
	    slide_right(idx, pos, start, stop, alen);
	    npos = (idx[1] - idx[0])/asize;
	    printf("start %d stop %d\n", start, stop);
	    printf("\t idx: %d %d npos: %d\n", idx[0], idx[1], npos);

	    start += wstep;
	   	stop += wstep;
	 }


}

int main(int argc, char **argv) {
    

//	printf("\n***Testing get position ***\n");
//	test_get_positions();

  printf("\n***Testing compare freq ***\n");
	test_compare_freq();
	printf("\n***Testing fill average ***\n");
	test_fill_averages();
	printf("\n***Testing compare all***\n");
	test_compare_all();
	printf("\n***Testing calculate distance***\n");
	test_calc_dist();
	printf("\n***Testing matrix mult***\n");
	test_matrix_mult();
	printf("\n***Testing cmds ***\n");
	test_cmds();
	printf("\n***Testing smacof ***\n");
	test_smacof();
	printf("\n***Testing css***\n");
	test_css();
	//printf("\n***Testing swap***\n");
	//test_swap();
	//printf("\n***Testing copy matrix ***\n");
	//test_copy_matrix();
	printf("\n***Testing stress ***\n");
	test_stress();
	printf("\n***Testing random shuffle ***\n");
	test_random_shuffle();
	printf("\n*** Testing setup z matrix ***\n");
	test_setup_z_matrix();
	printf("\n*** Testing matrix sqrt ***\n");
	test_matrix_sqrt();
	printf("\n*** Testing matrix squared ***\n");
	test_matrix_squared();
	printf("\n*** Testing part of sign tresh *** \n");
	test_sign_tresh();
	printf("\n*** Testing slide right ***\n");
	test_slide_right();

	//to run the full parallel test: run with Galaxy37 Galaxy38 parallel-outfile chrV-serial-results.txt
  test_all(argc, argv);
  test_thread(argc, argv);
    return EXIT_SUCCESS;
}


/* run:

gcc -c CSSS.c CSSS.h testCSSS.c
gcc -o css CSSS.o testCSSS.o -lgsl -lgslcblas -lm
./css
*/
