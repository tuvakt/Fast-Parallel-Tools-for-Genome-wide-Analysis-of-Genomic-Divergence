cimport numpy as np
import numpy as np

cdef extern from "threadcss.h":
     void threadcompute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, int treshold, int runs, int drosophila, int mds, double *scores, double *p)
     pass

cdef extern from "css.h":
     pass

cdef extern from "comparative.h":
     pass

def cluster_separation_scorer(np.ndarray[np.float64_t, ndim=1] avals, np.ndarray[np.float64_t, ndim=1] bvals, np.ndarray[np.int32_t, ndim=1] apos, np.ndarray[np.int32_t, ndim=1] bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, int treshold, int runs, int drosophila, int mds, np.ndarray[np.float64_t, ndim=1] scores, np.ndarray[np.float64_t, ndim=1] p):
    threadcompute(<double*> avals.data, <double*> bvals.data, <int*> apos.data, <int*> bpos.data, regstart, regend, wsize, wstep, alen, blen, treshold, runs, drosophila, mds, <double*> scores.data, <double*> p.data)