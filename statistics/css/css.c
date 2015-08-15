/*
 * css.c
 *
 * Author: tuvakt
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include "css.h"
#include "comparative.h"

/**
 * Computing the Cluster Separation Score (CSS) and the corresponding p-value for each window in a chromosome.
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

void compute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, int treshold, int runs,
		int drosophila, int mds, double *scores, double *p) {

	int i, start, stop, npos, asize, bsize, m, dims, wcount;
	double result;
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

	/* get number of individuals in group a and b */
	asize = get_population_size(apos);
	bsize = get_population_size(bpos);
	printf("asize %d\n", asize);
	printf("bsize %d\n", bsize);

	m = asize + bsize;
	dims = 2;
	wcount = 0; 	// adress in array is window_start/wstep
	start = 0;
	stop = wsize;
	result = 0;
	npos = 0;
	// to control number of calcs, only for debug
	//regend = start + wsize + 10*wstep;

	/* setup idx arrays */
	aidx = (int*)malloc(2*sizeof(int));
	bidx = (int*)malloc(2*sizeof(int));
	aidx[0] = 0; aidx[1] = 0;
	bidx[0] = 0; bidx[1] = 0;

	/* setup track arrays */
	atracks = (int*)malloc(asize*sizeof( int));
	btracks = (int*)malloc(bsize*sizeof( int));
	signtracks = ( int*)malloc(m*sizeof( int));
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

	unsigned short state[3] = {0,0,0};
	unsigned short seed = time(NULL);
	memcpy(state, &seed, sizeof(seed));
	srand48(time(NULL));

  while (start + wsize <= regend + wstep) {

    slide_right(aidx, apos, start, stop, alen);
    slide_right(bidx, bpos, start, stop, blen);
    npos = (aidx[1] - aidx[0])/asize;

    if (npos > 0) {
    	result = cluster_separation_scorer(distance, &(avals[aidx[0]]), &(bvals[bidx[0]]), atracks, btracks, asize, bsize, npos, dissimilarity,
    			drosophila, mds, X, B, Z, tmp, L, Q, r);

    	if (result != -1) {
    		scores[wcount] = result;
    		// not eligible for drosophila
    		p[wcount] = significance_treshold(distance, signtracks, asize, bsize, scores[wcount], treshold, runs, state);
    	}
    }
    start += wstep;
    stop += wstep;
    wcount++;

  }

  printf("npos in last window %d\n", npos);
  printf("wcount %d\n", wcount);

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

/**
 * Calculate a Cluster Separation Score (CSS) for one window in the chromosome
 * This calculation consists of three steps:
 * 1: Pairwise compare all individuals with a distance measure to get the dissimilarity matrix
 * 2: Use Multi-Dimensional Scaling (MDS) to scale down to two dimensions
 * 3: Calculate the cluster separation score
 * The calculations are based on Jones et al. and Vederhus.
 *
 * @param distance a (asize+bsize) x (asize+bsize) matrix with the distances ...
 * @param avals the values for group A in the window
 * @param bvals the values for group B in the window
 * @param atracks the ids for group A, in [0, asize-1]
 * @param btracks the ids for group B, in [asize, asize+bsize-1]
 * @param asize the size of group A
 * @param bsize the size of group B
 * @param npos number of SNP positions in the window
 * @param dissimilarity a (asize+bsize) x (asize+bsize) matrix with the dissimilarites between all the individuals
 * @param drosophila the choice of distance metric
 * @param mds the choice of MDS algorithm
 * @param X the init and result matrix from the MDS method
 * @param B, Z, tmp, L, Q, result: helper values for the MDS method
 * @return the CSS value for the window
 */
double cluster_separation_scorer(double **distance, double *avals, double* bvals, int *atracks, int *btracks, int asize, int bsize, int npos, double **dissimilarity,
		int drosophila, int mds, double **X, double **B, double **Z, double **tmp, double **L, double **Q, double **result) {
	int i, m, dims, valcount;

	m = asize + bsize; /* m: number of individuals */
	dims = 2;   /* dimension of MDS space */

	/* make sure diagonal is 0 */
    for (i = 0; i < m; i++) {
    	dissimilarity[i][i] = 0;
    }

	/* 1: Pairwise compare all individuals with a distance measure
	 * to get the dissimilarity matrix */
	if (drosophila) {
		compare_freq(avals, bvals, npos, dissimilarity);
	} else {
		compare_all(avals, bvals, asize, bsize, npos, dissimilarity);
	}

	/* fill avg and check if window should be discarded */
	if (!fill_averages(dissimilarity, m)) {
		return -1;
	}

	/* 2: Use MDS to scale down to two dimensions */
	if (mds == 0) {
		/* Classical MDS */
		cmds(dissimilarity, X, dims, m, B, Z, tmp, L, Q);
	} else if (mds == 1) {
		/* SMACOF */
		smacof_runs(dissimilarity, m, dims, X, Q, B, Z, result, 300, 4, 0.000001);
	} else {
		/* SMACOF + CMDS */
		cmds(dissimilarity, X, dims, m, B, Z, tmp, L, Q);
		smacof(dissimilarity, m, dims, X, Q, B, Z, 300, 0.000001);
	}

	/* 3: Calculate the cluster separation score - precalc. the distance */
	calc_dist(X, distance, m);
	return css(distance, atracks, btracks, asize, bsize);
}

/* gives the abs value in double
 * inlined for speed
 */
inline double dabs(double number) {
	if (number < 0) {
		return -1*number;
	}
	return number;
}

/**
 * Compute the distance between all individuals with distance metric
 * = the average of the frequency of the minor allele for each population
 * to create a dissimilarity matrix.
 * Used for Drosophila Melanogaster (Burke et al.)
 *
 * @param avals values for group A
 * @param bvals values for group B
 * @param compared the matrix to be filled with dissimilarity measures
 */
void compare_freq(double *avals, double *bvals, int npos, double **dissimilarity) {

	/* We have frequency of the minor allele for each population,
	 * so dissim. is a 2x2 matrix.
	 * The diagonal is 0 (a is equal to a)
	 * The rest contains the average of the (absolute value of the) difference
	 * between a and b
	 */
	int i;
	double avg = 0;

	for (i = npos; i--; ) {
		avg += dabs(avals[i] - bvals[i]);
	}

	avg = avg/npos;

	dissimilarity[0][1] = avg;
	dissimilarity[1][0] = avg;
}

/**
 * Compare all individuals with each other with distance metric = count differences
 * to create a dissimilarity matrix
 *
 * @param avals the values for group A
 * @param bvals the values for group B
 * @param asize the size of group A
 * @param bsize the size of group B
 * @param npos number of SNP positions in the window
 * @param dissimilarity the matrix to be filled with dissimilarity measures
 */
void compare_all(double *avals, double *bvals, int asize, int bsize, int npos, double **dissimilarity) {
  int i, j, k;
  double count;

  /* we skip i == j, because the dist(i,i) = 0
   * to increase speed, this function is longer than it should be
   * (we avoid unnecessary method calls)
   */

  /* for comparing avals with it self */
  for (i = asize; i--; ) {
	  for (j = i; j--; ) {
		  count = 0;
		  for (k = npos; k--; ) {
		    if (((int)avals[k*asize + i]*avals[k*asize + j]) == -9) {
				  count++;
			  }
		  }
		  dissimilarity[i][j] = count;
		  dissimilarity[j][i] = count;
	  }
  }

  /* bvals */
  for (i = bsize; i--; ) {
	  for (j = i; j--; ) {
		  count = 0;
		  for (k = npos; k--; ) {
		    if (((int)bvals[k*bsize + i]*bvals[k*bsize+j]) == -9) {
				  count++;
			  }
		  }
		  dissimilarity[i+asize][j+asize] = count;
		  dissimilarity[j+asize][i+asize] = count;
	  }
  }

  /* avals and bvals */
  for (i = asize; i--; ) {
	  for (j = bsize; j--; ) {
		  count = 0;
		  for (k = npos; k--; ) {
		    if (((int)avals[k*asize + i]*bvals[k*bsize+j]) == -9) {
				  count++;
			  }
		  }
		  dissimilarity[i][j+asize] = count;
		  dissimilarity[j+asize][i] = count;
	  }
  }
}

/**
 * Fill the empty positions in the dissimilarity matrix with the average value.
 * If more than half the elements in the matrix is 0, we discard the window
 *
 * @param dissimilarity a matrix to be filled with dissimilarity values
 * @param m the dimension of the matrix
 * @return 0 if the window is discarded, 1 otherwise
 */
int fill_averages(double **dissimilarity, int m) {
	int i, j, unval_pos = 0, totpos = m*m;
	double avg = 0;

	for (i = m; i--; ) {
		for (j = m; j--; ) {
			if (dissimilarity[i][j] < 0.00001) {
				unval_pos++;
			} else {
				avg += dissimilarity[i][j];
			}
		}
	}


	avg = avg/totpos;

	// if too many were 0, discard window
	if (unval_pos > totpos/2) {
	  return 0;
	}

	// else, keep on filling the blanks with averages
	for (i = m; i--; ) {
		for (j = m; j--; ) {
			if (dissimilarity[i][j] < 0.00001) {
				dissimilarity[i][j] = avg;
			}
		}
	}
	return 1;
}

/**
 * Allocate a m x n double matrix
 *
 * @param A a pointer to the double 2D matrix to be allocated
 * @param m number of rows
 * @param n number of columns
 */
void allocate_matrix(double ***A, int m, int n) {
	int i;
	double **data = (double**)malloc(m*sizeof(double*));
	if (data == NULL) {
		printf("'data' is NULL, something went wrong when mallocing\n");
		exit(1);
	}

	data[0] = (double*)malloc(m*n*sizeof(double));
	if (data[0] == NULL) {
		printf("'data[0]' is NULL, something went wrong when mallocing\n");
		exit(1);
	}

	for (i = 0; i < m; i++) {
		data[i] = &(data[0][i*n]);
	}

	*A = data;
}

/**
 * Deallocate a double 2D matrix
 *
 * @param A the matrix to be deallocated
 */
void deallocate_matrix(double **A) {
	free(A[0]);
	free(A);
}

/**
 * Multiply two matrices.
 * Computes C = A*B, where
 *
 * A is m x n,
 * B is n x p
 * and C is m x p
 *
 * @param A a m x n matrix to be multiplied
 * @param B a n x p matrix to be multiplied
 * @param C a m x p matrix where the result is stored
 * @param m the number of rows in matrix A
 * @param n the number of rows in matrix B and columns in matrix A
 * @param p the number of columns in matrix B
 */
void matrix_mult(double **A, double **B, double **C, int m, int n, int p) {
	gsl_matrix_view a = gsl_matrix_view_array(A[0], m, n);
	gsl_matrix_view b = gsl_matrix_view_array(B[0], n, p);
	gsl_matrix_view c = gsl_matrix_view_array(C[0], m, p);

	/* Compute C = A B */
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			1.0, &a.matrix, &b.matrix,
			0.0, &c.matrix);
}

/**
 * Square all the elements in a matrix
 * Computes B = A^2
 *
 * @param A the matrix to the squared
 * @param B the resulting matrix
 * @param m the dimension of the matrices
 */
void matrix_squared(double **A, double **B, int m) {

    int i, j;
    for (i = m; i--; ) {
    	for (j = m; j--; ) {
    		B[i][j] = A[i][j]*A[i][j];
       }
   }  
}

/**
 * Compute the square root of the matrix
 *
 * Computes Aij = sqrt(Aij), for all i,j
 *
 * @param A the matrix to be square rooted
 * @param m the dimension of the matrix
 */
void matrix_sqrt(double **A, int m) {
	int i, j;

	for (i = m; i--; ) {
		for (j = m; j--; ) {
			A[i][j] = sqrt(A[i][j]);
       }
   }  
}

/**
 * Set up the Z matrix needed in Classical MDS
 *
 * Computes Z = I - (1\m)*11^T
 *
 * where I is the identity matrix and 1 is a column vector with m ones.
 * m is the dimension of the dissimilarity matrix.
 *
 * @param A the matrix to be filled
 * @param m the dimension of the matrix
 */
void setup_z_matrix(double **A, int m) {
    int i, j;

    for (i = m; i--; ) {
    	for (j = m; j--; ) {
    		if (i != j) A[i][j] = -1.0/m;
    		else A[i][j] = (m-1)/(m*1.0);
    	}
    }
}

/**
 * Compute the Classical MDS (CMDS) from the dissimilarity matrix 'compared'.
 * The algorithm consist of several steps:
 * 1: Square the dissimilarity data
 * 2: Convert the squared dissimilarity matrix to scalar products through double centering
 * 3: Compute the eigen-decomposition B_delta = QLQ'
 * 4: Store the first m eigenvalues greater than 0 in L and the corresponding vectors in Q.
 * The solution of CMDS is X = QL^(1/2)
 *
 * @param dissimilarity the m x m dissimilarity matrix
 * @param X the matrix where the solution is stored
 * @param dims the dimension we are going to scale down to
 * @param m the dimension of the dissimilarity matrix
 * @param B,Z,tmp,L,Q helper matrices for the algorithm
 */
void cmds( double **dissimilarity, double **X, int dims, int m, double **B, double **Z, double **tmp, double **L, double **Q) {
    /* Classical MDS: Torgerson scaling
    Used if "one wants to assume or if one can prove that the dissimilarity data (...) are Euclidean distances". Here, I assume. 
       from http://download.springer.com/static/pdf/24/chp%253A10.1007%252F978-3-642-31848-1_8.pdf?auth66=1425904830_7c1a5d6c4eedd1ec6cfcd8b77716366a&ext=.pdf */
    int i, j, n;

    /* init */
    n = m;  /* distances square */

    /* 1: Square the dissimilarity data, saving in B */
    matrix_squared(dissimilarity, B, n);

    /* 2: Convert squared dissim. to scalar products through double centering */
    setup_z_matrix(Z, n);

    matrix_mult(Z, B, tmp, n, n, n);
    matrix_mult(tmp, Z, B, n, n, n);

    for (i = n; i--; ) {
    	for (j = n; j--; ) {
    		B[i][j] *= -0.5;
    	}
    }

    /* 3: Compute the eigen-decomposition B_delta = QLQ' 
     I only keep the 'dims' largest eigenvals + corresponding vec */
    memset(L[0], 0, dims*dims*sizeof(double));

    gsl_matrix_view ma = gsl_matrix_view_array(B[0], n, n);
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv (&ma.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    // sort in descending order, by eigen values
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
    
    for (i = dims; i--; ) {
      double eval_i = gsl_vector_get (eval, i);
      gsl_vector_view evec_i = gsl_matrix_column (evec, i);
      
      // TODO: check: only keep the 'dims' largest eigenvals > 0
      // This will not do, I suspect?
      L[i][i] = eval_i;
      for (j = n; j--; ) {
    	  Q[j][i] = gsl_vector_get(&evec_i.vector, j);
      }
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec); 

    /* 4: Solution of classical mds: X = QL^(1/2) */
    matrix_sqrt(L, dims);
    matrix_mult(Q, L, X, n, dims, dims);
}

/**
 * Calculate the Euclidean distance between every point in A
 *
 * dij = sqrt[sum_{x=1}^{dims} (Aix - Ajx)^2]
 *
 * here dims = 2, so x takes on the values [0,1]
 *
 * @param A a m x 2 matrix
 * @param distance the matrix to be filled with distances
 * @param m the number of rows in A (number of individuals)
 */
void calc_dist(double **A, double **distance, int m) {
  int i, j;
  double ans;
  
  // distance is symmetric
  for (i = m; i--; ) {
	  distance[i][i] = 0;
	  for (j = i; j--; ) {
		  ans = sqrt((A[i][0] - A[j][0])*(A[i][0] - A[j][0]) +
				  (A[i][1] - A[j][1])*(A[i][1] - A[j][1]));
		  distance[i][j] = ans;
		  distance[j][i] = ans;
	  }
  }
}

/**
 * Calculate the Cluster Separation Score (CSS) in the window
 * We calculate
 *
 * CSS = sum_{i} sum_{j} si,j/(m*n) - (m+n)*[sum_{i=1}^{m-1} si,i+1/(m^2*(m-1)) + sum_{j=1}^{n-1} sj,j+1/(n^2*(n-1))]
 *
 * where m and n are the size of the two populations,
 * i and j are individuals from the respective populations
 * and s is the Euclidean distance between a pair of individuals in the first two axes of the MDS space.
 *
 * see Jones et al. and the master's thesis
 *
 * @param distance a matrix of distances
 * @param atracks the ids in group A
 * @param btracks the ids in group B
 * @param a the size of group A
 * @param b the size of group B
 * @return the CSS value
 */
double css(double **distance,  int *atracks,  int *btracks, int asize, int bsize) {

	/*
	 * To increase speed, this function is not divided into smaller
	 * functions
	 */
	int i, j;
    double bet_dist, a_dist, b_dist, ab_val;

    // average between-group distance
    bet_dist = 0;
    for (i = asize; i--; ) {
    	for (j = bsize; j--; ) {
    		bet_dist += distance[atracks[i]][btracks[j]];
    	}
    }
    bet_dist = bet_dist/(asize*bsize);

    // within-group distance in population a
    a_dist = 0;
    if (asize > 1) {
    	for (i = (asize-1); i--; ) {
    		a_dist += distance[atracks[i]][atracks[i+1]];
        }
    	a_dist = a_dist/(asize*asize*(asize-1));
    }

    // within-group distance in population b
    b_dist = 0;
    if (bsize > 1) {
    	for (i = (bsize-1); i--; ) {
    		b_dist += distance[btracks[i]][btracks[i+1]];
        }
    	b_dist = b_dist/(bsize*bsize*(bsize-1));
    }
  
    ab_val = (asize+bsize)*(a_dist + b_dist);

    return bet_dist - ab_val;
}

/* swap two elements.
 * inlined for speed */
inline void swap( int *a,  int *b) {
     int tmp = *a;
    *a = *b;
    *b = tmp;
}


/*
Copyright (C) 2014 Torbjorn Rognes
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
Department of Informatics, University of Oslo,
PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

inline long random_int_nrand48(long n, unsigned short state[]) {
	/* Generate a random integer in the range 0 to n-1, inclusive
	 * 	The random() function returns a random number in the range
	 * 	0 to 2147483647 (=2^31-1=RAND_MAX), inclusive.
	 * 	We should avoid some of the upper generated numbers to
	 * 	avoid modulo bias.
	 *
	 * 	inlined for speed
	 */
	long random_max = RAND_MAX;
	long limit = random_max - (random_max + 1) % n;
	long r = nrand48(state);
	while (r > limit)
		r = nrand48(state);
	return r % n;
}

/**
 * Shuffle elements in an array
 * Fisher-Yates shuffling is used to shuffle the array
 *
 * @param elms an integer array of elements to be shuffled
 * @param n the length of the array
 * @param state seeds for the prng nrand48
 */
void random_shuffle( int *elms, int n, unsigned short state[]) {
    int i, r;
    for (i = n-1; i > 0; i--) {
    	r = random_int_nrand48(i+1, state); 
    	swap(&elms[i], &elms[r]);
    }
}


/**
 * Calculate the p-value for the CSS calculation.
 * We calculate the significance by a Monte Carlo method.
 * p is estimated by
 *
 * p = (r+1)/(n+1)
 *
 * r should at least be 10. See http://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
 *
 * @param the distances between the points
 * @param tracks a list of ids for group A and B, values in [0, asize+bsize-1]
 * @param asize size of group A
 * @param bsize size of group B
 * @param score the CSS value in the window
 * @param treshold the minimum number of iterations (r >= treshold)
 * @param runs the maximum number of iterations (n <= runs)
 * @param state seeds for the prng nrand48
 */
double significance_treshold(double **distance, int *tracks, int asize, int bsize, double score, int treshold, int runs, unsigned short state[]) {
    int i, hits, nscores, ntracks;
    double newscore, p;   
    int *atracks;
    int *btracks;

    ntracks = asize + bsize;

    hits = 0;
    nscores = 0;

    while (hits < treshold && nscores < runs) {
        random_shuffle(tracks, ntracks, state);
	    atracks = &(tracks[0]);
	    btracks = &(tracks[asize]);

	    newscore = css(distance, atracks, btracks, asize, bsize);
	    if (newscore >= score) {
	    	hits++;
		}
	    nscores++;
    }

    p = (hits+1)*1.0/(nscores+1);
    return p;
}

/*
 * Calculate the stress of X.
 * We calculate
 *
 * sigma(X) = sum_{i < j}[w_ij(d_ij(X) - delta_ij)^2]
 *
 * where d_ij(X) is the distance of object i and j in X
 * and delta_ij is the dissimilarity of object i and j.
 *
 * No data is missing, so w_ij = 1 for all i,j.
 *
 * Inlined for speed.
 */
inline double stress(double **dissimilarity, double **D, int m) {
	int i, j;
	double sigma = 0;

	for (i = m; i--; ) {
		for (j = i; j--;) {
			sigma += (D[i][j] - dissimilarity[i][j])*(D[i][j] - dissimilarity[i][j]);
		}
	}
	return sigma;
}

/*
 * Copy one matrix to another matrix
 * Since dims = 2, we drop the inner loop for more speed
 * Helper method for the SMACOF alg.
 *
 * Inlined for speed.
 */
inline void copy_matrix(double **Z, double **X, int m) {
	int i;
	for (i = m; i--; ) {
		memcpy(Z[i], X[i], 2*sizeof(double));
	}
}

/*
 * Compute the Guttman Transform for X.
 * We compute
 *
 * X^k = 1/m*B(Z)*Z
 *
 * where Z is X^{k-1}
 * and the elements of B(Z) is given by
 *
 * b_ij = - w_ij*delta_ij/d_ij   for d_ij != 0
 * b_ij = 0 					 for d_ij == 0
 *
 * for i != j, d_ij = d_ij(Z), and
 *
 * b_ii = sum_{j != i} b_ij
 *
 * Inlined for speed.
 */
inline void guttman_transform(double **X, double **B, double **Z, double **D, double **dissimilarity, int m, int dims) {
	int i, j;
	double d;

	for (i = m; i--; ) {
		d = 0;
		for (j = m; j--; ){
			if (i != j) {
				if (D[i][j] < 0.00001)
					B[i][j] = 0;
				else
					B[i][j] = -1*dissimilarity[i][j]/D[i][j];
				d += B[i][j];
			}
		}
		B[i][i] = -1*d;
	}

	// X^k = 1/m*B(Z)*Z
	// B: mxm, Z: mxdims
	matrix_mult(B, Z, X, m, m, dims);
	for (i = m; i--; ) {
		X[i][0] = X[i][0]/m;
		X[i][1] = X[i][1]/m;
	}
}

/*
 * Calculate the SMACOF algorithm with random start values
 * We run the SMACOF algorithm several times with different
 * random start values
 *
 * @param dissimilarity the dissimilarity matrix
 * @param m the dimensions of the dissimilarity matrix
 * @param dims the dimension we are going to scale down to
 * @param X the matrix with the results from the MDS alg.
 * @param Z,B,D,result helper matrices for the SMACOF alg.
 * @param max_iters maximum number of iterations
 * @param n_init number of random start values
 * @param epsilon a small number, used for convergence measure
 */
void smacof_runs(double **dissimilarity, int m, int dims, double **X, double **Z, double **B, double **D, double **result, int max_iters, int n_init, double epsilon) {
	/* Vederhus computed his mds by a metric SMACOF alg, python scikit.learn, with a random init and precomputed dissimilarity
 	 dissimilarity: numerical metric of how different two data objects are. We have 0 for equal, greater score for not equal, so have dissimilarity data
 	 eps = 0.001, max_iter = 300, n_init = 4 (run 4 times with different init, get the best of all four runs out */
	int runs, i;
	double sigma, sigma_best;
	sigma_best = 99999;
	for (runs = 0; runs < n_init; runs++) {

		// initialize 'result' with random numbers between 0 and 1
		for (i = 0; i < m; i++) {
			result[i][0] = drand48();
			result[i][1] = drand48();
		}

		// 'result' stores the solution from the smacof alg.
		sigma = smacof(dissimilarity, m, dims, result, Z, B, D, max_iters, epsilon);

		/*
		 * The best according to the stress function.
		 * Stress: minimizes the so-called stress function \sigma(X)
		 * The final results will be the best output of the n_init consecutive runs in terms of stress.
		 */
		//printf("k %d sigma %g sigma_prev - sigma %g\n", k, sigma, (sigma_prev - sigma));
		if (sigma < sigma_best) {
			// result is new best solution
			// store in X
			//printf("new better solution sigma %g sigma best %g\n", sigma, sigma_best);
			copy_matrix(X, result, m);
			sigma_best = sigma;
		}
	}
	//printf("final solution k %d sigma %g sigma_prev - sigma %g\n", k, sigma, (sigma_prev - sigma));
}

/*
 * Calculating the SMACOF (Scaling by MAjorizin a COmplicated Function) algorithm for Multi-Dimensional Scaling
 * The algorithm is taken from chap 8: Modern Multidimensional Scaling: Theory and Applications
 * by Ingwer Borg and Patrick J. F. Groenen
 *
 * The goal of SMACOF is to minimize Stress(X) for all X
 *
 * Stress(X) = sum_{i < j}[w_ij(d_ij(X) - delta_ij)^2]
 *
 * where d_ij(X) is the distance of object i and j in X
 * and delta_ij is the dissimilarity of object i and j
 *
 * No data is missing, so w_ij = 1 for all i,j
 *
 * @param dissimilarity the dissimilarity matrix
 * @param dims the number of dimensions to scale down to
 * @param X the start and end point of the algorithm. Contains the result.
 * @param Z,B,D helper matrices
 * @param max_iters maximum number of iterations
 * @param epsilon a small number, used for convergence measure
 */
double smacof(double **dissimilarity, int m, int dims, double **X, double **Z, double **B, double **D, int max_iters, double epsilon) {

	int i, j, k = 0;
	double sigma, sigma_prev, d;

	copy_matrix(Z, X, m);

	// Step 2: Compute sigma_r = sigma_r(X^0)
	calc_dist(X, D, m);
	sigma = stress(dissimilarity, D, m);
	//printf("Sigma start: %g\n", sigma);

	sigma_prev = 0;

	// Step 3: while ...
	while (k == 0 || ((sigma_prev - sigma) > epsilon && k <= max_iters)) {

		sigma_prev = sigma;
		// Step 4: increase k with 1
		k++;

		// Step 5: Compute the Guttman transform X^k
		guttman_transform(X, B, Z, D, dissimilarity, m, dims);

		// Step 6: sigma_r^k =sigma_r(X^k)
		calc_dist(X,D,m);
		sigma = stress(dissimilarity, D, m);

		// Step 7: Set Z = X^k
		copy_matrix(Z, X, m);
	}

	//printf("k %d sigma %g sigma_prev - sigma %g\n", k, sigma, (sigma_prev - sigma));
	return sigma;
}


