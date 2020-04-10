#!/usr/bin/env python

#	ct_compute_hi.py
#	2018-09-21	Markku-Juhani O. Saarinen <mjos@pqshield.com>
#	Copyright (C) 2019, PQShield Ltd. Please see LICENSE.

#	Compute the number of iterations that the constant time sampler requires
#	so that the weight is "h" with probability 1-2^-kappa

import numpy as np
import mpmath as mp

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,currentdir + "/parameter_definitions")

from chosen_parameter_sets import getParameterSets, R5_paramSet

# precision (in decimal digits) -- 100 is more than enough for fp > 2^-256
mp.mp.dps = 100

# base-2 logarithm string
def str_log2(x):
	x = mp.fabs(x)
	return ("2^" + mp.nstr(mp.log(x, 2), 6, strip_zeros=False)).ljust(12)

# given dimension d and weight h, determine number of "passes" required

def compute_hi(d, h,  alg, kappa=128, verbose=False):

    # failure probability threshold 2^-kappa;
    fp = mp.mpf(2)**(-kappa)

    # rejection rate of the uniform sampler part
    div = 65536 / d								# PARAMS_RS_DIV
    lim = div * d								# PARAMS_RS_LIM
    r = mp.mpf(lim)/mp.mpf(65536)

    # bp[w] is the "Bernoulli" transition probability from weight w to w+1.
    # computed as
    #       probability of sampling a suitable element (r)
    #       times
    #       probability of the element non being occupied (d-i/d)
    bp = np.array([r * mp.mpf(d-i)/mp.mpf(d) for i in range(h)])

    # weight probability, wp[i] = Pr(w(eight)=i).
    # The initial distribution is Pr(w=0) = 1, Pr(w>0) = 0
    wp = np.array([mp.mpf(0)] * (h+1))
    wp[0] = mp.mpf(1);

    # i is the number of passes (calls to drbg).
    # after one iteration (i=1),
    #       p[weight = 0] = p[w=0]* p[staying at w=0] = p[w=0]*(1-bp(0))
    #       p[weight = 1] = p[w=0]* p[moving to w=1] = p[w=0]*bp(0)
    #       p[weight > 1] = 0
    # after i iterations (i=1),
    #       p[weight = w] = p[w-1]*bp[w-1] +  p[w]*(1-bp(w))
    for i in range(9999):

        # end condition; Pr(w=h) > 1-fp or equivalently Pr(w<h)+Pr(w>h) < fp
        if (wp[h] > mp.mpf(1) - fp):
            break;

        # redistribute probability mass; can be done in place
        a = wp[0];
        wp[0] = mp.mpf(0);
        for j in range(h):
            b = wp[j + 1];
            wp[j] += (mp.mpf(1) - bp[j]) * a	# w     unchanged probability
            wp[j + 1] = bp[j] * a				# w++   probability
            a = b
        wp[h] += b								# wp[h] mass approaches 1

    # verbose information
    if verbose == True:
        print alg.ljust(20) + "d=" + str(d).ljust(6) + "h=" + str(h).ljust(6),
        print "HMAX=" + str(i).ljust(6) + "fp=" + str_log2(wp[0:h].sum())
    return i

#	parameters
for paramSet in getParameterSets():
    compute_hi(paramSet.d,     paramSet.h, paramSet.name, paramSet.kappa, True)
