import os,sys,inspect
import scipy.stats as ss
import numpy as np
import math
import mpmath as mp
import cmath

from chosen_parameter_sets import *

mp.mp.dps = 100

def chi2uppertail(x,k):
    x,k = mp.mpf(x), mp.mpf(k)
    #http://mpmath.org/doc/current/functions/expintegrals.html
    # upper_tail = gamma(k/2)/gamma(k/2) - lower_tail
    return 1 - mp.gammainc(k/2, 0, x/2, regularized=True)


def A_chi2_compute_p(target, numbins, verbose=False):
    cv = 0
    while mp.mpf(chi2uppertail(cv,numbins - 1 )) >  (mp.mpf(2)**(-target)):
        cv = cv + 1

    if verbose:
        print "Chi2: Critical value Level 2^(-" + repr(target) + ") = " + repr(cv)

    return cv


def binomial_test_compute_th(kappa, d , num_bins, verbose=True):
    pr = 1/float(num_bins)        # probability of a bin
    total_pr_threshold = mp.mpf(2)**(-kappa)
    total_pr = 0
    c = d
    while total_pr < total_pr_threshold:
        pr_c = mp.mpf(ss.binom.pmf(c,d,pr))
        total_pr += pr_c
        c -= 1
    if verbose:
        print "Binomial test: Threshold value Level 2^(-" + repr(kappa) + ") = " + repr(c+1)
    return c+1


def get_thresholds():
    for paramSet in getParameterSets():
        
        num_bins = paramSet.p
        thBinTest = binomial_test_compute_th(paramSet.kappa, paramSet.d, num_bins, False)
        print paramSet.name
        print "   binomial test: ", paramSet.p, thBinTest
        
        ## code computing the binomial threshold when a histogram is used
        num_bins = 2**int(math.log(paramSet.d/10, 2))
        #thBinTest2 = binomial_test_compute_th(paramSet.kappa, paramSet.d, num_bins, False)
        #print "   binomial test (2):", num_bins, paramSet.p/num_bins, thBinTest2
        
        cv = A_chi2_compute_p(paramSet.kappa, num_bins, False)
        print "   chi2 test        :", num_bins, cv


get_thresholds()



