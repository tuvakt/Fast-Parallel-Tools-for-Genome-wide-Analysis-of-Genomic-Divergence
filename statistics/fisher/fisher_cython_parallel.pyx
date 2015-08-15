cimport numpy as np
import numpy as np

cdef extern from "threadfisher.h":
     void threadcompute(double *avals, double *bvals, int *apos, int *bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, double perc, double *scores, double *stddev)
     pass

cdef extern from "cFisher.h":
     pass

cdef extern from "comparative.h":
     pass

def fisher_exact_tester(np.ndarray[np.float64_t, ndim=1] avals, np.ndarray[np.float64_t, ndim=1] bvals, np.ndarray[np.int32_t, ndim=1] apos, np.ndarray[np.int32_t, ndim=1] bpos, int regstart, int regend, int wsize, int wstep, int alen, int blen, double perc, np.ndarray[np.float64_t, ndim=1] scores, np.ndarray[np.float64_t, ndim=1] stddev):
    threadcompute(<double*> avals.data, <double*> bvals.data, <int*> apos.data, <int*> bpos.data, regstart, regend, wsize, wstep, alen, blen, perc, <double*> scores.data, <double*> stddev.data)