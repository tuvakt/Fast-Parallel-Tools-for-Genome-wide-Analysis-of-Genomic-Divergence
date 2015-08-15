/*
 * cFisher.c
 *
 * Author: tuvakt
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include "cFisher.h"
#include "comparative.h"

/**
 * Compute the Fisher Exact Test (FET) described by Burke et. al. for the entire chromosome
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
void compute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen,
	     double perc, double *scores, double *stddev) {

	struct timeval before, after;
	gettimeofday(&before, NULL);
	int start, stop, npos, asize, bsize, nsamples, wcount;
	int *aidx;
	int *bidx;
	int *f;
	int *tmp;
	double *results;
	double *fetscores = NULL;
	double *samples = NULL;
	double *stdsamples;

	/* get size of a and b*/
	asize = get_population_size(apos);
	bsize = get_population_size(bpos);
	printf("asize %d\n", asize);
	printf("bsize %d\n", bsize);

	wcount = 0;
	start = 0;
	stop = wsize;
	nsamples = 100;

	// for debug
	//regend = start + wsize + 8*wstep;

	/* setup idx arrays */
	aidx = (int*)malloc(2*sizeof(int));
	bidx = (int*)malloc(2*sizeof(int));
	aidx[0] = 0; aidx[1] = 0;
	bidx[0] = 0; bidx[1] = 0;

	/* kommenter der de def. */
  	f = (int*)malloc(4*sizeof(int));
  	tmp = (int*)malloc(4*sizeof(int));
  	results = (double*)malloc(2*sizeof(double));
  	stdsamples = (double*)malloc(nsamples*sizeof(double));


  	/* prng */
  	unsigned short state[3] = {0,0,0};
  	unsigned short seed = time(NULL);
  	memcpy(state, &seed, sizeof(seed));

  	while (start + wsize <= regend + wstep) {

  		slide_right(aidx, apos, start, stop, alen);
  		slide_right(bidx, bpos, start, stop, blen);
  		npos = (aidx[1] - aidx[0])/asize;

  		if (npos > 0) {
  			// reallocing
  			// In case that ptr is a null pointer, the function behaves like malloc, assigning a new block of size bytes and returning a pointer to its beginning.
  			fetscores = (double*)realloc(fetscores, npos*sizeof(double));
  		  	samples = (double*)realloc(samples, npos*sizeof(double));
  		  fisher_exact_test(results, &(avals[aidx[0]]), &(bvals[bidx[0]]), asize, bsize, npos, f, tmp, samples, stdsamples, nsamples, fetscores, state, perc);
		  scores[wcount] = results[0];
		  stddev[wcount] = results[1];
  		}
  		start += wstep;
  		stop += wstep;
  		wcount++;
  	}

  	printf("wcount %d\n", wcount);

  	/* cleanup */
  	free(aidx);
  	free(bidx);
  	free(f);
  	free(tmp);
  	free(results);
  	free(samples);
  	free(stdsamples);
  	free(fetscores);

  	gettimeofday(&after, NULL);
  	printf("Time %g\n", time_ddiff(before, after));
}


/*
 * comparison metric of two double values
 * used by the qsort alg. in percentile
 *
 * inlined for speed
 */
inline int compare_doubles (const void *a, const void *b) {
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

/*
 * Get the percentile score
 *
 * Inlined for speed
 */
inline double percentile(double *fetscores, int n, double percentile) {
	int idx;
	double delta;
    /* sort inplace with quick sort */
    qsort(fetscores, n, sizeof(double), compare_doubles);
    idx = (n-1)*percentile; // indices from 0 to n-1
    delta = (n-1)*percentile - idx;
    return (1-delta)*fetscores[idx] + delta*fetscores[idx+1];
}

/**
 * Calculate the Fisher Exact Test (FET) score and standard deviation for one window.
 * We calculate the L10FET score for each SNP position in the window
 * and select the perc % percentile as the FET score for the entire window (L10FET perc%Q score).
 * The standard deviation of nsamples bootstrap replicates of L10FET perc%Q are also calculated.
 * The calculations are taken from Burke et. al:
 * http://www.nature.com/nature/journal/v467/n7315/abs/nature09352.html
 * Burke et. al. used perc = 0.95 and nsamples = 100
 *
 * @param results a double array of size 2 for the results of the function
 * @param avals values for group a in the window
 * @param bvals values for group b in the window
 * @param asize number of individuals in group a
 * @param bsize number of individuals in group b
 * @param npos number of SNP positions in the window
 * @param f a 2x2 frequency table, stored as an 1d array
 * @param tmp an integer array of size 4: used for the contingency table
 * @param samples a double array of size npos - used for calculating the std. dev.
 * @stdsamples a double array of size nsamples - used for calculating the std. dev
 * @nsamples number of bootstrap samples - used for calculating the std. dev
 * @fetscores npos x 1 array to hold the scores for each SNP position
 * @state a seed for the prng nrand48
 * perc the percentile score in the window. Burke et. al. used perc = 0.95
 */
void fisher_exact_test(double *results, double *avals, double *bvals, int asize, int bsize, int npos, int *f, int *tmp,
		       double *samples, double *stdsamples, int nsamples, double *fetscores, unsigned short state[], double perc) {

	int i;
	double fetscore;

	/* For each SNP position,"
	count frequencies and
   	calculate the negative log 10 fet score
	For every SNP, we calculated -log10(P) from a Fisher’s exact test (L10FET)." (Burke et al).*/
	for (i = npos; i--; ) {
		fetcount(f, avals, bvals, i, asize, bsize);
		fetscore = fet(f, tmp);
		fetscores[i] = -1.0*log10(fetscore);
	}
	/* the score at perc% percentile
	  indicates the importance of the window
	*/
	results[0] = percentile(fetscores, npos, perc);
	/* calculate std.dev. for the scores.
	 Burke et al. estimated std by calculating
     std of 100 boostrap replicate samples of L10FET5%Q
	 */

	results[1] = calc_std(fetscores, samples, stdsamples, nsamples, npos, perc, state);
}


/**
 * Create a 2x2 frequency table from the data
 *
 * @param f an empty frequency table of size 4, to be filled
 * @param avals the values for group a
 * @param bvals the values for group b
 * @param idx: SNP position in window
 * @param asize the number of individuals in group a
 * @param bsize: number of individuals in group b
 */
void fetcount(int *f, double *avals, double *bvals, int idx, int asize, int bsize) {

  int i, majsum, minsum, astart, bstart;

  astart = idx*asize;
  bstart = idx*bsize;

  majsum = 0; minsum = 0;
  for (i = astart; i < (astart+asize); i++) {
    if (avals[i] == 3) {
      majsum++;
    } else if (avals[i] == -3) {
      minsum++;
    }
  }
  
  f[0] = majsum;
  f[1] = minsum;
  
  majsum = 0; minsum = 0;
  for (i = bstart; i < (bstart+bsize); i++) {
    if (bvals[i] == 3) {
      majsum++;
    } else if (bvals[i] == -3) {
      minsum++;
    }
  }
  
  f[2] = majsum;
  f[3] = minsum;
}
 
/* Rosettacode
http://rosettacode.org/wiki/Evaluate_binomial_coefficients#C
We go to some effort to handle overflow situations

write more

*/
 
inline unsigned long gcd_ui(unsigned long x, unsigned long y) {
  unsigned long tmp;
  if (y < x) {
      tmp = x; x = y; y = tmp;   /*swap x and y*/
  }
  while (y > 0) {
      tmp = y;  y = x % y;  x = tmp;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}

unsigned long binomial(unsigned long n, unsigned long k) {
    /* formula: result = sum (n+1-i)/i */

    unsigned long i, result, g;

    if (k == 0 || k == n) return 1;
    if (k == 1) return n;
    if (k > n) return 0;
    if (k > n/2) k = n-k;

    result = 1;
    for (i = 1; i <= k; i++) {
	if (result >= ULONG_MAX/n) {  /* Max value for unsigned long -> Possible overflow */
	    /* we have to reduce the numbers */
	    unsigned long n_reduced, i_reduced;  /* reduced numerator / denominator */
	    g = gcd_ui(n, i);  n_reduced = n/g;  i_reduced = i/g;
	    g = gcd_ui(result, i_reduced);  result = result/g;  i_reduced = i_reduced/g;
	    if (result >= ULONG_MAX/n_reduced) return 0;  /* Unavoidable overflow */
	    result *= n_reduced;
	    result /= i_reduced;
	} else {
	    result *= n;
	    result /= i;
	}
		n--;
    }

    return result;
}

/**
 * Find the index corresponding to the smallest value in the array
 *
 * @param a an integer array
 * @param n the size of the array
 */
int min_idx(int *a, int n) {
    int i, idx, min;
    idx = 0;
    min = a[idx];
    for (i = 1; i < n; i++) {
    	if (a[i] < min) {
    		min = a[i];
    		idx = i;
    	}
    }
    return idx;
}

/*
 * Swap the elements.
 * Inlined for speed
 */
inline void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

/**
 * Shift table so that is becomes
 *   a0 | b0
 *   c0 | d0
 *   with a0 <= b0, c0, d0
 *
 *   This is a helper function for the metod fet.
 *
 *   @param f an integer array, with the 2x2 frequency table to be shifted
 *   @param tmp an integer array
 *   @param n the size of the array f and tmp
 */
void shift_table(int *f, int *tmp, int n) {
    int i, idx;

    /* make sure freqs stored in clockwise order: 0132 */
    swap(&f[2], &f[3]);

    /* find idx of smallest value, this will be a0 */
    idx = min_idx(f, n);

    /* fill f with values, in order
       f[0] a0, f[1] b0, f[2] c0 f[3] d0 */
    for (i = n; i--; ) {
    	tmp[i] = f[(idx+i)%n];
    }
    swap(&tmp[2], &tmp[3]);

    for (i = n; i--; ) {
      f[i] = tmp[i];
    }
}

/**
 * Create the most extreme table for the second tail,
 * see Feldman & Klinger (1963) and Zar (1987)
 *
 * This is a helper method for the method fet
 *
 * @param f an integer array for the 2x2 frequency table to be created
 * @param n the size of the array
 */
void create_table(int *f, int n) {
    int R1, R2, C1, C2, m1, midx;
   
    R1 = f[0] + f[1];
    R2 = f[2] + f[3];
    C1 = f[0] + f[2];
    C2 = f[1] + f[3];

    int A[4] = {R1, R2, C1, C2};
    midx = min_idx(A, 4);
    m1 = A[midx];

    if (R1 <= R2 && C1 <= C2) {
    	f[0] = m1 - f[0];
    	f[1] = R1 - f[0];
    	f[2] = C1 - f[0];
    	f[3] = C2 - f[1];
    } else if(R1 <= R2 && C2 <= C1) {
    	f[1] = m1 - f[1];
    	f[0] = R1 - f[1];
    	f[3] = C2 - f[1];
    	f[2] = C1 - f[0];
    } else if (R1 >= R2 && C1 <= C2) {
    	f[2] = m1 - f[2];
    	f[0] = C1 - f[2];
    	f[3] = R2 - f[2];
    	f[1] = R1 - f[0];
    } else {
    	f[3] = m1 - f[3];
    	f[1] = C2 - f[3];
    	f[2] = R2 - f[3];
    	f[0] = R1 - f[1];
    }
}

/**
 * Calulate the two-tailed Fisher Exact Test (FET) for the 2x2 contingency table f
 * Uses the short cut calculation proposed by
 * Feldman & Klinger(1963): Short cut calculation of the fisher-yates "Exact test"
 * http://link.springer.com/article/10.1007%2FBF02289576
 * and elaborated by
 * Zar (1987): A fast and efficient algorithm for the Fisher exact test
 * http://link.springer.com/article/10.3758/BF03202590
 *
 * @param f an array of size 4, the 2x2 contingency table
 * @param tmp an array of size 4
 * @return the two tailed FET
 */
double fet(int *f, int *tmp) {
    int n, R1, R2, C1, C2;
    double P, P0, P1, P2;

    R1 = f[0] + f[1]; R2 = f[2] + f[3];
    C1 = f[0] + f[2]; C2 = f[1] + f[3];
    
    n = 4;
    shift_table(f, tmp, n);

    /* calculate the one-tailed probability, P1 */
    /* calculate P0 of original table */
    P0 = fet_p(f[0], f[1], f[2], f[3]);
    //printf("P original table %g\n", P0);
    P = P0;
    P1 = P;
    /* while a0 >= 0, go on */
    while (f[0] > 0) {
    	f[1]++; f[2]++;   /* bi+1 = bi + 1, ci+1 = ci + 1; */
		P1 = (1.0*f[0]*f[3])/(f[1]*f[2])*P1;
		P += P1;
		f[0]--; f[3]--;   /* ai+1 = ai - 1; di+1 = di - 1; */
    }

    /* calculate the two tailed prob by calc the prob of the other extreme */

    /* if the row/column totals are equal
       the two-tailed prob is 2*one tailed prob. */
    if (R1 == R2 || C1 == C2) {
       	P = 2*P;
    } else {
    	create_table(f, n);
    	shift_table(f, tmp, n);
    	P2 = fet_p(f[0], f[1], f[2], f[3]);
	
    	while (P2 < P0) {
    		P += P2;
    		if (f[1] == 0 || f[2] == 0) {
    			break;
    		}
    		f[0]++; f[3]++; /* ai+1 = ai + 1, di+1 = di + 1; */
    		P2 = (1.0*f[1]*f[2])/(f[0]*f[3])*P2;
    		f[1]--; f[2]--; /*bi+1 = bi - 1; ci+1 = ci - 1; */
		}
    }
    /* TODO: What about this ? */
    if (P > 1)
    	P = 1;
    return P;
    
}

/**
 * Calculate the probability (p) of the 2x2 frequency table
 * a | b
 * c | d
 *
 * p = (a+b)!(c+d)!(a+c)!(b+d)!/(a!b!c!d!n!)
 *	 = binomial(a+b, a)*binomial(c+d,c)/binomial(n, a+c)
 *
 * where n = a+b+c+d
 *
 * @param a
 * @param b
 * @param c
 * @param d
 * @return the p-value
 */
double fet_p(int a, int b, int c, int d) {
    int ab, cd, ac, n;
    ab = a+b;
    cd = c+d;
    ac = a+c;
    n = ab + cd;
    double nom, denom;
    nom = binomial(ab, a) * binomial(cd, c);
    denom = binomial(n, ac);
    return nom/denom;
}

/**
 * Calculate the standard deviation
 *
 * @param a an double array
 * @param n the size of the array
 *
 */
double std(double *a, int n) {
    int i;
    double mu, sum;

    mu = mean(a, n);
    sum = 0;
    for (i = n; i--; ) {
    	sum += (a[i] - mu)*(a[i] - mu);
    }
    sum /= n;
    return sqrt(sum);
}

/**
 * Calculate the mean value
 *
 * @param a an double array
 * @param n the size of the array
 */
double mean(double *a, int n) {
    int i;
    double sum = 0;
    for (i = n; i--; ) {
    	sum += a[i];
    }
    return sum/n;
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
/*
 * Generate a random integer in the range 0 to n-1, inclusive.
 * The random() function returns a random number in the range
 * 0 to 2147483647 (=2^31-1=RAND_MAX), inclusive.
 * We should avoid some of the upper generated numbers to
 * avoid modulo bias.
 *
 * Inlined for speed.
 *
 */
inline int long random_int_nrand48(long n, unsigned short state[]) {
	long random_max = RAND_MAX;
	long limit = random_max - (random_max + 1) % n;
	long r = nrand48(state);
	while (r > limit)
		r = nrand48(state);
	return r % n;
}

/*
 * Random sampling with replacement
 * select n samples
 *
 * Inlined for speed
 */
inline void bootstrap_sample(double *fetscores, double *samples, int n, int npos, unsigned short state[]) {
    /* n: number of samples 
       npos: number of scores */
    int i, idx;
    for (i = n; i--; ) {
    	// draw a random index from 0 to npos-1
    	idx = random_int_nrand48(npos, state);
    	samples[i] = fetscores[idx];
    }
}


/**
 * Calculate the standard deviation of nsamples bootstrap replicate scores of L10FET perc%Q.
 * This is done by sampling 'npos' bootstrap samples from 'fetscores'. We calculate the
 * perc% percentile of these scores and store the resulting L10FET perc%Q scores in 'stdsamples'
 * Sample npos boostrap samples from 'fetscores'.
 * When we have 'nsamples' different percentile values,
 * we calculate the standard deviation from this.
 *
 * @param fetscores the L10FET scores for each SNP position in the window
 * @param samples an double array of size 'npos'
 * @param stdsamples an double array of size 'nsamples'
 * @param nsamples number of samples for std. dev calculation
 * @param npos number of SNP positions in the window
 * @param perc the percentile of the L10FET scores
 * @param state a seed for the prng nrand48
 */
double calc_std(double *fetscores, double *samples, double *stdsamples, int nsamples, int npos, double perc, unsigned short state[]) {
	int i, j;
	for (i = 0; i < nsamples; i++) {
		bootstrap_sample(fetscores, samples, npos, npos, state);
		stdsamples[i] = percentile(samples, npos, perc);
	}
	return std(stdsamples, nsamples);
}


