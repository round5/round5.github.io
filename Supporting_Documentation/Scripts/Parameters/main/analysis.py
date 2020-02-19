# -*- coding: utf-8 -*-
#!/usr/bin/env python
from __future__ import division
from math import log, ceil, floor, sqrt, pi, exp, factorial as fac
import os,subprocess,time,datetime,platform
import operator as op
# from sympy import isprime,nextprime
# from collections import defaultdict
from ast import literal_eval
from decimal import Decimal
import sys, os
import numpy as np
import mpmath as mp

log_infinity=999999999



#----------------  Probability distribution convolution routines  -----------------
#
#    Copyright (c) 2019, PQShield and Koninklijke Philips N.V.
#    Markku-Juhani O. Saarinen, Koninklijke Philips N.V.
#
#    Failure probability computation for Round5. Custom-made, super fast.

# epsilon in probability computations
prq_eps = mp.mpf(2)**(-300)

# precision (number of decimal digits)
mp.mp.dps = 100

# just a (very precise) zero vector

def mp_zeros(l):
    return np.array([mp.mpf(0)] * l)

# Computation of v = v*f where v has lenght q and f is a uniform probability generation function f
# with n non-zero components each of size 1/n. d is "scaling". The support is {0,d,2d,..,(n-1)d}.
# The scaled function f'(t) and the original function f(t) are related as f'(t/d) = f.
# It requires n to be a power of 2 and (n-1)d < q.
#
#  Function f with scaling 1
#   1/n ***********
#       |
#       |
#       |
#       |__________***************
#                  n              q

def conv_comb(v, n, d):
    l = len(v)
    i = 1
    while i < n:
        v += np.roll(v, (d * i) % l)
        i <<= 1
    return v / n

# Convolution of input function v with two probability generating functions with support
# {-q/2p,-q/2p+1,..,q/2p-1} and {-q/2p+1-q/2p+1,..,q/2p}.
# Used to find the distribution of the sum of h/2 variables distributed as U{-q/2p,-q/2p+1,..q/2p -1}
# minus the sum of h/2 such variables
# This is computed as
#
#  < ---------        h/2     --------- >
# (v+ - v-) * (v+ - v-) * ... * (v+ - v-)
#
def spec_conv_comb_qp(v,q,p):
    w=conv_comb(v,q//p,1)   # first convolution  v = v * f where f = pgf of U({0,1,..,n=q/p-1})
    w=np.roll(w, -q//(2*p)) # shift to the left to recenter it (less probable values in the center around q/2)  so that the convolution is  with pfg of U({-q/2p,-q/2p+1,..,q/2p-1})
    w=conv_comb(w,q//p,1)   # second convolution v = v * f where f = pgf of U({0,1,..,n=q/p-1})
    w=np.roll(w,1-q//(2*p)) # shift to the left to recenter it (less probable values in the center around q/2)  so that the convolution is  with pfg of U({-q/2p+1,-q/2p+2,…,q/2p})
    return w


# convolution of input function v with uniform varibles with rounding noise t->p scaled q/t in support q.
def spec_conv_comb_pt(v,q,p,t):
    w = conv_comb(v,p//t,q//p)      # convolution v = v * f where f= pgf of U(({j *(q/p) | 0 <= j <= (p/t)-1})
    w = np.roll(w, q//p-(q//(2*t))) # shift to the left to recenter it (less probable values in the center around q/2) so now w = v* pgf of U({j*(q/p) + q/p –(q/2t) | 0<=j<=(p/t)-1} = v* pgf of U({ (p/q)*k | 1-(p/2t) <= k <= p/2t})
    return(w)

def conv_repetition(v,q):
    w = np.zeros(q)
    w[0] = v[0]
    w[q//2] = v[q//2]
    for i in range(1,q//2-1,1):
        w[i] = v[i] + v[q-i]
    w = np.convolve(w,w)
    return w

# add_coeff is true if and only if we consider RLWR without special reduction

def round5_fp_fast(n, h, q, p, t, f, mu, B, repetition_code, add_coeff=False, verbose=False ):
    ## use this for full precision
    ###   dif = mp_zeros(q)
    ###   dif[0] = mp.mpf(1.0)
    ## use this for double precision
    dif = np.zeros(q)
    dif[0] = 1.0
    # determine number of times to convolve the term differences in term s*e.
    #If prime cyclotomic, it is h, if ring swithing or scaler, just h/2.
    # Only h(2h) convolutions as we start with probability gen function
    # for difference between two random variables with q-p rounding noise
    # This can be done as we consider balanced secrets
    w=h;
    if add_coeff:
        w=2*w;
    # computation of convolution of w variables (we have the sum of w variables)
    for i in range(w):
        dif = spec_conv_comb_qp(dif, q, p)
    # convert to mpf if we were double
    dif = np.array([mp.mpf(dif[i]) for i in range(len(dif))])
    # normalize, just in case
    dif /= dif.sum()
    # second convolution with the p to t noise. Single time to account for coefficient compression from p to t.
    dif = spec_conv_comb_pt(dif, q, p, t)
    # if information is encoded twice (as in Newhope), then, convolve again. normalize, just in case.

    if repetition_code == 1:
        # Decryption failure ranges for case when parameter B>1
        up_fail = (q + (1<<B)) // (1<<(B+1))
        down_fail = (q*( (1<<(B+1)) - 1) + (1<<B)) // (1<<(B+1))
        # Bit failure probability (bfp). Distribution goes from 0 to q-1, highly provable values are centered around 0, thus, we sum the central part.
        bfp = dif[up_fail:down_fail].sum();
    else:
        dif = conv_repetition(dif, q)
        bfp = dif[q//2:dif.size].sum();

    # Failure probability after error correction
    ffp = mp.mpf(0);
    for j in range(f + 1, mu + 1):
        ffp += (mp.binomial(mu, j) * mp.power(bfp, j) *
                mp.power(mp.mpf(1) - bfp, (mu - j)))
    return float(mp.log(ffp, 2))

#----------------  utilities.py (Reference:https://github.com/pq-crystals/kyber/tree/master/scripts)  -----------------


def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


# Create probability distribution for a sparse, ternary distribution
# with Hamming weight h. +1, -1 each occur with probability h/2n,
# 0 occurs with probability 1-(h/n).
# The parameter "sparseness" denotes the ratio hamming_weight/n,
# where n is the length of the vector and hamming_weight the number of non-zeroes in it.
def build_sparse_ternary_law(sparseness):
    # Sanity check. sparseness should be a fraction
    if sparseness>=1.:
        print "build_sparse_ternary_law() received a sparseness parameter >= 1. Abort."
        exit()
    D = {}
    domain=[-1,0,1]
    D[domain[0]] = float(sparseness/2.)
    D[domain[1]] = float(1.0 - sparseness)
    D[domain[2]] = float(sparseness/2.)
    return D


# Binary search on a non-decreasing function. Adapted from John Schank's PARI GP code.
# func(x) non-decreasing; maxLT() finds maximum x such that f(x) <= target
# Parameter func is a lambda function, defined in and passed as parameter from the calling function.
def maxLT(func,target,low,high):
    #print "\nmaxLT(), func is",func,"target is",target,"low is",low,"high is",high
    if(low >= (high-1)):
        return low
    midpt = ceil((low+high)/2)
    if(func(midpt)<=target):
        return maxLT(func,target,midpt,high)
    else:
        return maxLT(func,target,low,midpt)
    # endif


# Get the largest in the domain of a given PDF, that has non-zero probability
# If the PDF is not zero-centered, we first make it so.
def max_coeff(coeffDist):
    nonzero_keys = []
    for key in coeffDist.keys():
        if coeffDist[key]!=0: nonzero_keys.append(key)
    abs_keys = [abs(nonzero_keys[i]) for i in range(len(nonzero_keys))]
    return max(abs_keys)


#-------------------------------------------------  secanalysis_functions.py  ------------------------------------------------------------


#Calculate log base 2 of an input
def log2(n):
    res=float(log(n)/log(2))
    return res


#Calculate the root Hermite factor for a given BKZ blocksize
def get_rhf(b):
    if((b>=2) and (b<50)):    # BKZ with block-size=2 is
                              # exactly the LLL algorithm
        return 1.0219         # best known practical value
                              # of LLL's root Hermite factor.
                              # Theoretical upper bound of LLL's RHF is 1.0754.
                              # See Gama's "Predicting Lattice Reduction".
    else:                     # BKZ is the generalization of LLL to block-sizes>2.
                              # Relation to obtain RHF from b only makes sense for b>50
        return ( (pi*b)**(1./b) * b / (2*pi*exp(1)))**(1./(2.*b-2.))

# log_2 of best plausible Generic Cost of SVP in dimension b. For random lattices,
# hence has an extra 'b' multiplicative factor as compared to Ideal lattices.
def svp_core_sieving(b, ring, exponent):
    if(ring==0):
        return (float(log2(b)) + float(b*exponent)) # Cost for random lattices,
                                                    # note the additional log b
    else:
        return float(b*exponent)

# log_2 of best plausible Quantum Cost of SVP in dimension b. For random lattices,
# hence has an extra 'b' multiplicative factor as compared to Ideal lattices.
def svp_quantum_core_sieving(b, ring):
    return svp_core_sieving(b, ring, 0.265)

# log_2 of best plausible Classical Cost of SVP in dimension b. For random lattices,
# hence has an extra 'b' multiplicative factor as compared to Ideal lattices.
def svp_classical_core_sieving(b, ring):
    return svp_core_sieving(b, ring, 0.292)

# Analogous to Q-core Enum + O(1) of "Estimate all the {LWE,NTRU} Schemes" paper.
# Used to estimate Enum attack-cost in NTRU-HRSS-KEM.
def svp_quantum_enumeration(b, ring=0):
    return float((0.18728*b*log(b, 2) - 1.0192*b + 16.1)/2.)

# Classical version of the above Enumeration attack.
def svp_classical_enumeration(b, ring=0, exponent=None, dim=None, rounds=1):
    return float(0.18728*b*log(b, 2) - 1.0192*b + 16.1)


# Calculate the standard deviation of the LWR rounding-down + scaling-down error
def lwr_error_sd(q,p):
    return sqrt( ((((q*1.)/(p*1.))**2) - 1.)/12.)


# Calculate the norm of the i-th Gram-Schmidt vector in a reduced q-ary
# lattice basis of dimension d, using the Geometric Series Assumption.
# Lattice may be scaled by a scaling factor to exploit small secrets
def gram_schmidt_norm(d,m,q,rhf,i,omega):
    return Decimal(q**(m/d)) * Decimal(omega**((d-m-1)/d)) * (rhf**Decimal(d-(2*(i-1))))


# Calculate the log of the scaling/weighting factor used to build
# weighted lattices for attacks on small-secret LWE
def log_scaling_factor(sigma,theta,is_hnf=False,is_binary=False,is_ternary=False,is_sparsehnf=False,is_spbinary=False,is_spternary=True):
    if(is_hnf):
        log_w = 0.0
    if(is_binary):
        log_w = float(1.0+log2(sigma))
    if(is_ternary):
        log_w = float(log2(sigma)+(0.5*(log2(3.0) - 1)))
    if(is_sparsehnf):
        log_w = 0.0             #TODO: Confirm if this is correct!!
    if(is_spbinary):
        log_w = float(log2(sigma)-(0.5*(log2(theta) + log2(1.0-theta))))
    if(is_spternary):
        log_w = float(log2(sigma)-(0.5*log2(theta)))
    return log_w


# Calculate the log of the scaling factor specifically for the Dual attack,
# computed by balancing norms of secret and error per-vector, and not per-component.
def log_dual_scaling_factor(sigma,n,m,hamm_wt,is_hnf=False,is_binary=False,is_ternary=False,is_sparsehnf=False,is_spbinary=False,is_spternary=True):
    if(is_spternary):
        log_w = float( log2(sigma)+(0.5*(log2(m)-log2(hamm_wt))) )
    else:
        # WARNING!!! See assumption^.
        log_w = 0.0
    return log_w


#Calculate the norm of the unique shortest vector (s,e,1) in the Primal Embedding lattice
def minNorm(n,m,sigma,theta,omega,is_hnf=False,is_binary=False,is_ternary=False,is_sparsehnf=False,is_spbinary=False,is_spternary=True):
    h=floor(n*theta)
    norm=0.0
    if(is_hnf):
        norm=sigma*sqrt(n+m)
    if(is_binary):
        norm=sqrt(float((omega**2)*(n/2)) + float(m*(sigma**2)))
    if(is_ternary):
        norm=sqrt(float(((omega**2)*2*n)/3) + float(m*(sigma**2)))
    if(is_sparsehnf):
        norm=sigma*sqrt(h+m)
    if(is_spbinary):
        norm=sqrt(float((omega**2)*h) + float(m*(sigma**2)))
    if(is_spternary):
        norm=sqrt(float((omega**2)*h) + float(m*(sigma**2)))
    return float(norm)


#Calculate the norm of the projection \widetilde(v) of the vector v=(s,e,1) on the vector space spanned by the last b-Gram-Schmidt vectors of the embedding lattice.
#v_norm: The norm of v=(s,e,1); v_dim: its dimension.
def projectedNorm(v_norm,b,v_dim):
    return float(v_norm*sqrt(float(b/v_dim)))


#Calculate the log of the volume of the lattice into which (primal) embedding is done, scaled or otherwise by a scaling/weighting factor omega.
#Calculating in log and not absolute values because q^d is huge and might quickly go out of hand.
def log_embedLatticeVol(q,d,n,omega):
    log_vol=float((d-n-1)*(log2(q))) + float(n*log2(omega))
    return log_vol


#Success condition of Primal attack against (R)LWE. Returns true if attack succeeds for given parameters and false otherwise.
def primal_attack_success(sd,b,q,n,m,theta,is_hnf=False,is_binary=False,is_ternary=False,is_sparsehnf=False,is_spbinary=False,is_spternary=True):
    #Dimension of the lattice considered in the Primal attack.
    d=m+n+1
    #Root-Hermite factor for a BKZ lattice reduction algorithm running with given block-size b.
    log_rhf=log2(get_rhf(b))
    #Log2 of scaling factor for rescaling the lattice in case of small-secrets.
    log_omega=log_scaling_factor(sd,theta,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary)
    omega=float(2**log_omega)
    #Norm of the unique shortest vector (s,e,1) that the primal attack searches for.
    norm_v=minNorm(n,m,sd,theta,omega,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary)
    dim_v=m+n+1
    #Norm of the projection of (s,e,1) on the span of the last b Gram-Schmidt vectors
    proj_norm_v=projectedNorm(norm_v,b,dim_v)
    #Log2 of volume of the lattice considered in the Primal attack, possibly scaled due to small-secrets.
    log_latt_vol=log_embedLatticeVol(q,d,n,omega)
    #Check embedding success condition
    lhs=log2(proj_norm_v)
    rhs=float((b+b-d)*log_rhf) + float(log_latt_vol/d)
    if(lhs<=rhs):
        #Primal attack succeeded
        return True
    else:
        #Primal attack failed
        return False


#for a given q,p,n calculate the cost of running the primal attack (if it succeeds), otherwise return infinity
def primalcost(b,q,p,n,m,theta,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,ring,svp_cost_model):
    #Calcuate the root Hermite factor
    rhf=get_rhf(b)
    #Calculate the standard deviation of the LWR rounding error
    sigma=lwr_error_sd(q,p)
    #Calculate primal attack success condition
    if(primal_attack_success(sigma,b,q,n,m,theta,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary)):
        # Primal attack succeeds. Calculate the actual cost of running it, assuming one SVP call
        return svp_cost_model(b, ring)
    else:
        return log_infinity     #Primal attack fails. Cost of running it is undefinably large.


# Primal attack analysis against a LWR noisy El-Gamal encryption scheme,
# with ciphertext of the form (u,v) \in ZZ_p^m x ZZ_t^\mu, where \mu<=m.
# Further, the component v has the message added to it, as (t/2^B)m.
def primalcost_ct(b,q,p,t,n,h,m,mu,B,ring,svp_cost_model):

    # Assuming sparse-ternary secrets.


    # Scaling factors for the different parts of the short(est)
    # vector v_w = ( w1.e1,w2.(q.2^B/p)e2,w3.1,w4.r )

    qbyp = 1<<int(log(q,2.)-log(p,2.))                                     # q/p as a power of 2
    pbyt = 1<<int(log(p,2.)-log(t,2.))                                     # p/t as a power of 2
    sigma_qp = sqrt( (qbyp*qbyp - 1.)/12. )
    sigma_pt = sqrt( (pbyt*pbyt - 1.)/12. )
    w = max( sigma_qp, 1.*qbyp*(2**B)*sigma_pt, 1, sqrt(h*1./n)  )
    w1 = w*1./sigma_qp                                                  # scales e1
    w2 = w*1. / (qbyp*(2**B)*sigma_pt)                                  # scales q.2^B/p e2
    w3 = w                                                              # scales c
    w4 = w*1./sqrt(h*1./n)                                              # scales the secret r


    # Norm of the short(est) vector v_w

    norm_vw = sqrt( ( ((w1*sigma_qp)**2)*m ) +\
                    ( (( (2**B)*qbyp*w2*sigma_pt )**2)*mu ) +\
                    ( w3**2 ) +\
                    ( (w4**2)*h ) )


    # Norm of the projection of v_w on the last b Gram-Schmidt vectors

    norm_vb = sqrt( (1.*norm_vw*norm_vw*b)/(m+mu+1+d) )


    # Volume of the rescaled lattice, in log2

    log_latt_vol = (1.*m*(log(q,2.)+log(w1,2.))) + (1.*mu*log(w2,2.)) +\
                   log(w3,2.) + (1.*d*log(w4,2.))


    # Norm of the (d-b)th Gram-Schmidt vector, in log2

    log_rhf = log(get_rhf(b),2.)
    log_gs_norm = (1.*(b+b-m-mu-1-d-1)*log_rhf) +\
                  (1.*log_latt_vol/(m+mu+1+d))

    if log(norm_vb,2.)<log_gs_norm:
        # Primal attack on ciphertext succeeds
        return svp_cost_model(b, ring)
    else:
        # Primal attack fails
        return log_infinity


def optimize_primalattack_ct(q,p,t,d,h,mu,
                             B,ring,
                             svp_cost_model,
                             max_m,           # either for u or v, depending on flags below
                             attack_u,        # Optimal samples are computed from u
                             attack_v):       # Optimal samples are computed from v
    best_cost=log_infinity                    # Init
    best_b=2                                  # Optimal BKZ blocksize
    best_m=1                                  # Optimal samples

    # Use tiny steps if not doing a fast search
    b_step = 1
    m_step = 1

    b_min = 60
    b_max = d+d+mu+1                          # conservative
    if(b_max<b_min):
        return (0,0,0)

    if fastsearch:
        b_step = (b_max - b_min) >> 4
    if b_step < 1:
        b_step = 1
    print "\n\nAttacking u:",attack_u,"\nattack v:",attack_v,"\nb_min:",b_min,"\tb_max:",b_max,"\tb_step:",b_step,"\nmax_m:",max_m,"\n"

    while b_step >= 1:
        for b in range(b_min, b_max + b_step, b_step):
            current_cost = svp_cost_model(b,ring)
            if(current_cost>best_cost):
                break

            ms_min = max(1,(b-d))
            ms_max = max_m

            if fastsearch:
                m_step = (ms_min - ms_max) >> 4
            if m_step < 1:
                m_step = 1

            while m_step >= 1:

                local_best = log_infinity
                local_opt = -1

                for m in range(ms_min, ms_max + m_step, m_step):

                    if attack_u:
		        attack_cost = primalcost_ct(b,q,p,t,d,h,m,mu,B,
                                                    ring,
                                                    svp_cost_model)
                    else:
                        attack_cost = primalcost_ct(b,q,p,t,d,h,d,m,B,
                                                    ring,
                                                    svp_cost_model)

                    if (attack_cost < local_best):
                        local_best = attack_cost
                        local_opt = m

                if (local_best < best_cost):
                    # Have found a b,m pair for which the attack cost
                    # is even cheaper than what it was so far.
                    best_cost = local_best
                    best_b = b
                    best_m = local_opt

                # adjust m range based on local best
                if (local_opt >= 0):
                    ms_dist = (ms_max - ms_min + 4) >> 2
                    ms_min = max(ms_min, local_opt - ms_dist)
                    ms_max = min(ms_max, local_opt + ms_dist)
                else:
                    break       # nothing found in this iteratpn

                m_step >>= 1

        # adjust the ranges based of best found..
        b_dist = (b_max - b_min + 4) >> 2
        b_min = max(b_min, best_b - b_dist)
        b_max = min(b_max, best_b + b_dist)
        b_step >>= 1

    # Recompute
    if attack_u:    correct_cost = primalcost_ct(best_b,q,p,t,d,h,best_m,mu,B,ring,svp_cost_model)
    else:           correct_cost = primalcost_ct(best_b,q,p,t,d,h,d,best_m,B,ring,svp_cost_model)

    if((best_b != 2) and (best_cost != correct_cost)):
        return (-1,-1,-1)
    else:
        return (best_b,best_m,best_cost)



# Optimize the primal attack on a Round5 ciphertext
def analyze_primalattack_ct(d,n,h,
                            q,p,t,
                            B,nbar,mbar,mu,
                            bkzexp,
                            fastsearch,
                            svp_cost_model):
    if d==n:
        ring=1
    else:
        ring=0

    # Check which of u,v has the greater noise
    qbyp = 1<<int(log(q,2.)-log(p,2.))           # q/p as a power of 2
    pbyt = 1<<int(log(p,2.)-log(t,2.))           # p/t as a power of 2
    sigma_qp = sqrt( (qbyp*qbyp - 1.)/12. )
    sigma_pt = sqrt( (pbyt*pbyt - 1.)/12. )
    attack_u = False
    attack_v = False
    if sigma_qp < sigma_pt:
        # Attack the weaker component
        attack_u = True
    else:
        attack_v = True

    if attack_u:       max_m = d              # For each secret r, U contains max. d samples in both ring and nonring cases
    else:
        if ring==0:    max_m = nbar           # For each secret r, v contains nbar samples in the nonring case
        else:          max_m = mu             # For each secret r, v contains mu samples in the ring case

    return optimize_primalattack_ct(q,p,t,d,h,mu,B,ring,svp_cost_model,max_m,attack_u,attack_v)



#Calculation of the distinguishing advantage of the Dual Attack and the corresponding attack cost, in log_2
#Note: the parameter "use_lwr_reduction" determines whether the (practical) dual attack cost (of solving Decision LWR) is calculated,
#or the (theoretical) dual attack cost (of solving Search LWE by reducing it to Decision LWR,
#via the Bogdanov reduction, Theorem 3 - https://eprint.iacr.org/2015/769.pdf) is calculated.
def dualcost(b,q,p,n,m,theta,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,ring,use_lwr_reduction,svp_cost_model):
    #Dimension of the dual lattice
    d=m+n
    #Root-Hermite factor for a BKZ lattice reduction algorithm running with given block-size b.
    rhf=get_rhf(b)
    #Calculate the standard deviation of the LWR rounding error
    sd=lwr_error_sd(q,p)
    #Log2 of scaling factor for rescaling the lattice in case of small-secrets.
    log_omega=log_dual_scaling_factor(sd,n,m,int(floor(theta*n)),is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary)
    omega=float(2**log_omega)
    #Norm of the shortest vector in the scaled dual lattice. omega is the scaling factor calculated above.
    norm_v=(rhf**(d-1)) * ((q/omega)**(1. * n/d))
    tau= norm_v * sd / q
    log2_eps = -0.5 + (- 2 * pi * pi * tau * tau) / log(2)
    cost_of_svp=svp_cost_model(b,ring)
    total_cost = max(0, - 2 * log2_eps - svp_core_sieving(b,ring,exponent=0.2075)) + cost_of_svp

    security_gap=0.0    #Default value
    #Sanity check. Currently, no efficient reductions to DECISION Ring-LWR are known.
    if((use_lwr_reduction==1) and (ring!=0)):
        print "Error.Cannot calculate Ring-LWE advantage based on Ring-LWR advantage. No efficient reductions for Decision-RLWR known. Abort!"
        exit()
    elif((use_lwr_reduction==1) and (ring==0)):
        #print "Calculating cost of solving Search-LWE by reducing to Decision-LWR."
        E=abs(float((2.**(log2_eps))/(4*q*m)) - float(2.**(n-(m*log2(p)))))    #Since E is getting squared, can disregard the negative sign to prevent log() from throwing UNDEFINED errors.
        if(E==0):
            log2_eps_lwe=0
        else:
            log_E_sq=2*log2(E)
            log2_eps_lwe=log_E_sq - m
        total_cost_lwe = max(0, - 2 * log2_eps_lwe - svp_core_sieving(b,ring,exponent=0.2075)) + cost_of_svp
        security_gap=total_cost_lwe-total_cost
    return (total_cost,security_gap)


## Hybrid meet-in-the-middle and lattice reduction attack analysis ##


# Calculate the cost of guessing (either via Grover's search or via meet-in-the-middle) the non-zero components of a vector of length r, whose non-zero components are distributed according to the probability distribution coeffDist.
def mitm_cost(coeffDist,r):
    # Can use either Grover's search (quantum approach) or MITM (classical approach) to search r coefficients
    # using 2^(0.5 * entropy(coeffDist)) queries
    plogp = 0
    for symbol in coeffDist:
        probab=coeffDist.get(symbol,0)
        plogp+=(probab * log2(probab))
    cost_mitm=floor(0.5*r*(-plogp))
    return cost_mitm


# Set the root Hermite factor requirement such that the last Gram Schmidt vector in the reduced lattice basis (of dimension d)
# has length greater than twice the maximum coefficient size in coeffDist, the distribution of the target vector's components.
def hybrid_hermite(d,m,q,coeffDist,omega):
    #WARNING!! Assumning a sparse, ternary distribution; maximum coefficient size is 1.
    max_coeff_sz = max_coeff(coeffDist)    # Get the largest in the domain of the given PDF, that has non-zero probability.
    min_gs_vec=2*omega*max_coeff_sz
    rhf_min=float(float((m/d)*log2(q)) + float(((d-m-1)/d)*log2(omega)) - (log2(min_gs_vec)))/d
    return Decimal(2**rhf_min)       #using Decimal data type for higher accuracy


# Back-calculate to estimate the largest BKZ block-size that results in a given root Hermite factor
def hybrid_blocksz(dim,rhf):
    bs = maxLT(lambda x: (-1*get_rhf(x)),
        -rhf, 2, dim)
    return (bs+1)


def hybrid_cost_estimate_nist(q,p,n,max_m,secretkey_law_param,coeffDist,r, ring, is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary, svp_cost_model):
    sigma=lwr_error_sd(q,p)    #updated non-0-centered LWR rounding error per P. Nguyen's email
    omega=float(2**(log_scaling_factor(sigma,secretkey_law_param,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary)))
    ## Find an optimal (m,RHF) pair satisfying the following two conditions:
    ## 1. The first G.S. vector has norm at most q, i.e., this is a necessary condition that lattice reduction succeeded in reducing the lattice basis.
    ## 2. The last G.S. vector has norm greater than twice the largest coefficient in coeffDist, the distribution of the target vector's components.

    ## max (m, h). m s.t. b*_1 < q; h s.t. b*_{m+(n-k)} > s
    ## Notice that with different m's tried out, omega remains the same since it is a function of constant n,q,p,theta
    m_opt = maxLT(lambda m: gram_schmidt_norm(m+(n-r),m,q,hybrid_hermite(m+(n-r),m,q,coeffDist,omega),1,omega),
        q, 2, max_m)
    d_opt=m_opt+(n-r)
    rhf_opt=hybrid_hermite(d_opt,m_opt,q,coeffDist,omega)
    blocksz_opt=hybrid_blocksz(d_opt,rhf_opt)
    cost_lr=svp_cost_model(blocksz_opt, ring)
    cost_MITM=mitm_cost(coeffDist,r)
    return (cost_lr,int(floor(blocksz_opt)))


# Calculate the optimal meet-in-the-middle dimension r that results in (almost) eqaluuy balanced meet-in-the-middle and lattice reduction costs.
# Parameter coeffDist is a dictionary representing the distribution of the secret distribution, containing a (symbol)->(probability) mapping.
def hybrid_tradeoff_nist(q,p,n,max_m,secretkey_law_param,coeffDist,ring, is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary, svp_cost_model):
    r_opt = maxLT(lambda r: (mitm_cost(coeffDist,r) - hybrid_cost_estimate_nist(q,p,n,max_m,secretkey_law_param,coeffDist,r,ring,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary, svp_cost_model)[0]),
        0, 0, n)
    return r_opt


# Allowed secret-key distributions (represented by the last parameter): Sparse-tri(/er)nary secrets, Centered binomial.
# Default SK distribution is the first.
# For the 1st secret_key distribution, the parameter "secretkey_law_param" denotes the ratio hamming_weight/vector_length.
# For the 2nd secret_key distribution, secretkey_law_param denotes the parameter of the Centered binomial distribution.
def hybrid_attack_cost_r_nist(q,p,n,max_m,secretkey_law_param,ring,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,svp_cost_model,secretkey_law=build_sparse_ternary_law):
    #probability distribution for the coefficients of the secret-key in Round2, implemented as a Python dictionary
    coeffDist=secretkey_law(secretkey_law_param)
    #calculate optimal meet-in-the-middle dimension by trying to balance MITM cost and lattice reduction cost in the hybrid attack.
    r=hybrid_tradeoff_nist(q,p,n,max_m,secretkey_law_param,coeffDist,ring,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,svp_cost_model)
    cost_mitm=mitm_cost(coeffDist,r)
    (cost_lr,blocksz_opt)=hybrid_cost_estimate_nist(q,p,n,max_m,secretkey_law_param,coeffDist,r,ring,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,svp_cost_model)
    return (cost_mitm,cost_lr,r,blocksz_opt)


## END Hybrid meet-in-the-middle and lattice reduction attack analysis ##


# For a given lattice instance, optimize over BKZ blocksize b and number of (R)LWR samples m to get the minimal primal or dual attack cost.
def optimizeattack_b_m(q,p,n,max_m,theta,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,attack_type,ring,use_lwr_reduction,fastsearch,svp_cost_model):
    if (attack_type != 1) and (attack_type != 2):
        print "Unknown optimizeattack_b_m() attack_type = ", attack_type
        exit()

    best_cost=log_infinity    #Initialize optimal attack cost to infinity.
    best_b=2                  #Least possible value     #Optimal BKZ blocksize that results in minimum attack cost.
    best_m=1                  #Least possible value     #Optimal number of LWR samples that results in minimum attack cost.
    b_min=60                  #Smallest value of BKZ blocksize that we consider.
    #Below value for the maximum value of b that can be considered as available to an attacker has been taken from NewHope.
    dual_attack_gap=0.0       #Security gap in the cost of solving Decision LWR via dual attack, versus the cost of
                              #solving Search LWE by reducing it to Decision LWR.
    b_max=(2*n)+max_m+1
    if(b_max<b_min):
        return (0,0,0)

    # Use tiny steps if not doing a fast search
    b_step = 1
    m_step = 1

    if fastsearch:
        b_step = (b_max - b_min) >> 4
    if b_step < 1:
        b_step = 1

    while b_step >= 1:

        for b in range(b_min, b_max + b_step, b_step):
            current_cost = svp_cost_model(b,ring)
            if(current_cost>best_cost):
                break        #Cost of running the attack is higher than the (previously) best cost. No point checking any more values of b,m since minimal cost has already been found

            ms_min = max(1,(b-n))
            ms_max = max_m

            if fastsearch:
                m_step = (ms_min - ms_max) >> 4
            if m_step < 1:
                m_step = 1

            while m_step >= 1:

                local_best = log_infinity
                local_opt = -1

#               print "m_step = ", m_step

                for m in range(ms_min, ms_max + m_step, m_step):

                    if(attack_type==1):
                        #Primal attack
                        attack_cost = primalcost(b, q, p, n, m, theta, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary, ring, svp_cost_model)
                    elif(attack_type==2):
                        #Dual attack (XXX this was commented out)
                        (attack_cost, dual_attack_gap) = dualcost(b, q, p, n, m, theta, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary, ring, use_lwr_reduction, svp_cost_model)

                    if (attack_cost < local_best):
                        local_best = attack_cost
                        local_opt = m


                if (local_best < best_cost):
                    #Have found a b,m pair for which the attack cost is even cheaper than what it was so far.
                    best_cost = local_best
                    best_b = b
                    best_m = local_opt

                # adjust m range based on local best
                if (local_opt >= 0):
                    ms_dist = (ms_max - ms_min + 4) >> 2
                    ms_min = max(ms_min, local_opt - ms_dist)
                    ms_max = min(ms_max, local_opt + ms_dist)
                else:
                    break       # nothing found in this iteratpn

                m_step >>= 1

        # adjust the ranges based of best found..
        b_dist = (b_max - b_min + 4) >> 2
        b_min = max(b_min, best_b - b_dist)
        b_max = min(b_max, best_b + b_dist)
        b_step >>= 1


    if  (attack_type==1):
        correct_cost = primalcost(best_b, q, p, n, best_m, theta, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary, ring, svp_cost_model)
    else:
        (correct_cost, dual_attack_gap) = dualcost(best_b, q, p, n, best_m, theta, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary, ring, use_lwr_reduction, svp_cost_model)

    if((best_b != 2) and (best_cost != correct_cost)):
        return (-1,-1,-1,-1)
    else:
        return (best_b,best_m,best_cost,dual_attack_gap)

## Updated analysis of SILKE


# Calculation of the distinguishing advantage of the Extended Dual Attack and the corresponding attack cost, in log_2
# Additional attack parameters:  k, the number of guessed components, and ell, the maximum Hamming weight of the guess
# Implicit assumptions: 0<=k<=n, 0<= ell <= min(h,k) and m<=n. Here h is the Hamming weight
# Note: the parameter "use_lwr_reduction" determines whether the (practical) dual attack cost (of solving Decision LWR) is calculated,
# or the (theoretical) dual attack cost (of solving Search LWE by reducing it to Decision LWR,
# via the Bogdanov reduction, Theorem 3 - https://eprint.iacr.org/2015/769.pdf) is calculated.
def dualcost_ell(b, q, p, n, m, theta, k, ell, log_som, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary,
                 ring, use_lwr_reduction, svp_cost_model, log_guessing_cost, log_parallelization_power, EDA_BKZ_LLL, LLL_cost):
    
   
    # Hamming weight of the secrets
    h = int(floor(theta * n))
    # Dimension of the dual lattice
    d = m + n - k
    # Root-Hermite factor for a BKZ lattice reduction algorithm running with given block-size b.
    rhf = get_rhf(b)
    # Calculate the standard deviation of the LWR rounding error
    sd = lwr_error_sd(q, p)
    # Log2 of scaling factor for rescaling the lattice in case of small-secrets.
    log_omega = log_dual_scaling_factor(sd, n, m, h - ell, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary)
    #log_omega = log_dual_scaling_factor(sd, n, m, h, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary)
    # This is to verify that this leads to a higher cost.
    
    omega = float(2 ** log_omega)
    # Norm of the shortest vector in the scaled dual lattice. omega is the scaling factor calculated above.
    norm_v = (rhf ** (d - 1)) * ((q / omega) ** (1. * (n - k) / d))
    
    if EDA_BKZ_LLL == True:
        norm_v = 1.1*norm_v
    
    tau = norm_v * sd / q
    log2_eps = -0.5 + (- 2 * pi * pi * tau * tau) / log(2)

    # Lattice contribution
    cost_of_svp = svp_cost_model(b, ring)
    
    if EDA_BKZ_LLL == False:
        total_cost = max(0, - 2 * log2_eps - svp_core_sieving(b, ring, exponent=0.2075)) + cost_of_svp
    else:
        total_cost = cost_of_svp


    # Guessing contribution
    # direct computation results in float division by zero, hence condition on som
    #som = sum(binomial(k, w) * 2 ** w for w in range(ell + 1))
    #log_som = log(som, 2)

    if EDA_BKZ_LLL == True:
        log_som = log_som + LLL_cost

    log_guess_cost = log_som - 2 * log2_eps
    #log_guessing_cost = 0 # this is the cost of guessing.
    log_guess_cost += log_guessing_cost
    #log_parallelization_power = 8 # this is cost reduction due to potential parallelization.
    log_guess_cost -= log_parallelization_power
    
    total_cost = max(log_guess_cost, total_cost)
    # approximate log(guess_cost+ 2**(total cost)) by max(log_guess_cost,total cost)
    # Taking into account that the secret vector has weight larger than l in k positions
    som = sum(binomial(h, w) * binomial(n - h, k - w) for w in range(ell + 1))
    r = ceil(binomial(n, k) / som)
    total_cost = total_cost + log(r, 2)

    return (total_cost)

#function to find optimal parameters for guessing_dual attack using binary search
def optimize_guessing_dual(b_minn, q, p, n, max_samples, theta, is_hnf, is_binary, is_ternary,
                                           is_sparsehnf, is_spbinary, is_spternary, ring, use_lwr_reduction,
                                           svp_cost_model, EDA_BKZ_LLL):
   
    mincost = 99999;
#    return (mincost, 99, 99, 99, 99,99,99,99,99)

    
    h = int(floor(theta * n))

    b = b_minn;

    log_guessing_cost = 0           # guessing  has a cost of 2^{log_guessing_cost} units
    log_parallelization_power = 0   # there are 2^{log_parallelization_power} performing guesses in parallel
    LLL_cost = 0                    # cost of LLL assuming BKZ + LLL
    
    # Initial step sizes
    step = 256
    lstep = 32
    
    speed = 2
    
    bmin = 20 #b_minn
    bmax = n
    kmin = 1
    kmax = n - h
    mmin = 1
    mmax = max_samples + 1
    
    ## using binary search
    (bopt, kopt, lopt, mopt) = (b_minn, 0, 0, 0)

    while (step > 1):
        for b in range(bmin,bmax,step):
            target = svp_cost_model(b, ring) -  0.2075 * b
            for k in range(kmin, kmax, step):
                u = min(k, h - 1)
                for ell in range (0, u, lstep):
                    som = sum(binomial(k, w) * 2 ** w for w in range(ell+1))
                    log_som = log(som, 2)
                    if int(log_som) < mincost:
                        for m in range(mmin, mmax, step):
                            newcost = dualcost_ell(b, q, p, n, m, theta, k, ell, log_som, is_hnf, is_binary, is_ternary,
                                                   is_sparsehnf, is_spbinary, is_spternary, ring, use_lwr_reduction, svp_cost_model,
                                                   log_guessing_cost, log_parallelization_power, EDA_BKZ_LLL, LLL_cost)
                            if (newcost < mincost):
                                mincost = newcost
                                (bopt, kopt, lopt, mopt) = (b, k, ell, m)

        step = max(1, int(step/speed))
        lstep = max(1, int(lstep/speed))
        
        # update values
        bmin = max (b_minn, bopt - 2*step)
        bmax = min (n, bopt + 2*step)
        kmin = max (1, kopt - 2*step)
        kmax = min (n - h, kopt + 2*step)
        mmin = max (1, mopt - 2*step)
        mmax = min (max_samples + 1, mopt + 2*step)
        
    ropt = log(binomial(n, kopt), 2) - log(sum(binomial(h, w) * binomial(n - h, kopt - w) for w in range(lopt + 1)), 2)

    d = mopt + n - kopt
    # Root-Hermite factor for a BKZ lattice reduction algorithm running with given block-size b.
    rhf = get_rhf(bopt)
    # Calculate the standard deviation of the LWR rounding error
    sd = lwr_error_sd(q, p)
    # Log2 of scaling factor for rescaling the lattice in case of small-secrets.
    log_omega = log_dual_scaling_factor(sd, n, mopt, h - lopt, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary)
    omega = float(2 ** log_omega)
    # Norm of the shortest vector in the scaled dual lattice. omega is the scaling factor calculated above.
    norm_v = (rhf ** (d - 1)) * ((q / omega) ** (1. * (n - kopt) / d))
    tau = norm_v * sd / q
    log2_eps = -0.5 + (- 2 * pi * pi * tau * tau) / log(2)
    cost_of_svp = svp_cost_model(bopt, ring)
    bkzcost = max(0, - 2 * log2_eps - svp_core_sieving(bopt, ring, exponent=0.2075)) + cost_of_svp
    # add the costs for guessing
    # direct computation results in float division by zero, hence condition on som
    som = sum(binomial(kopt, w) * 2 ** w for w in range(lopt + 1))
    log_guess_cost = log(som, 2) - 2* log2_eps
    
    return (mincost, bopt, kopt, lopt, mopt,log_guess_cost,bkzcost,ropt,cost_of_svp)




