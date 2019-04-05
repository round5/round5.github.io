#!/usr/bin/env python
from __future__ import division
from math import log, ceil, floor, sqrt, pi, exp, factorial as fac
import os,subprocess,time,datetime,platform
import operator as op
# from sympy import isprime,nextprime
# from collections import defaultdict
from ast import literal_eval
from decimal import Decimal
import sys
import numpy as np
import mpmath as mp

log_infinity=999999999


#----------------  Probability distribution convolution routines, author: Markku-Juhani O. Saarinen <mjos@iki.fi>  -----------------

#   r5fp.py
#   2018-08-25  Markku-Juhani O. Saarinen <mjos@iki.fi>
#   round5_parameters.py
#   2018-07-18  Sauvik Bhattacharya <sauvik.bhattacharya@philips.com>

#   Error estimates for Round5 via convolutions.

# epsilon in probability computations
prq_eps = mp.mpf(2)**(-300)

# precision (number of decimal digits)
mp.mp.dps = 100

# just a (very precise) zero vector

def mp_zeros(l):
    return np.array([mp.mpf(0)] * l)

# base-2 logarithm string

def str_log2(x):
    x = mp.fabs(x)
    if x <= prq_eps:
        return "0".ljust(12)
    return ("2^" + mp.nstr(mp.log(x, 2), 6, strip_zeros=False)).ljust(12)

# gives variance, 1-sum = eps, and sum over "failure" range, as a string

def str_var(v):
    l = len(v)
    var = mp.mpf(0)
    eps = 1.0
    for i in range (-l/2,l/2):
        var += v[i] * (i * i)
        eps -= v[i]
    if mp.fabs(eps) > prq_eps:
        return "!!! eps=" + str_log2(eps)
        exit()

    return mp.nstr(var).ljust(12)

# gives variance, 1-sum = eps, and sum over "failure" range, as a string

def str_stats(v):
    l = len(v)
    return ("var=" + str_var(v) +
            "bfp=" + str_log2(v[(l/4):(3*l/4)].sum()))


# failure probability computation. Custom-made, super fast.

'''
    Idea for "comb" convolutions is similar to polynomial multiplications:

         2**n-1           n-1
    X * ( SUM x^i ) = X * PROD  (x^(2**i) + 1)
          i=0             i=0

    Followed by a centering and mirroring etc.

    Algorithm is entirely additive so smaller mantissa precision is required.
'''

# comb length n, "scaling" d

def conv_comb(v, n, d):
    l = len(v)
    i = 1
    while i < n:
        v += np.roll(v, (d * i) % l)
        i <<= 1
    return v / n

# q->p rounding error multiplied with with +-d distribution

# Warning. Below, it is needed that 2|t|p|q because of
# Python's stupid division

def conv_comb_qp(v, q, p, d):
    v = conv_comb(v, q / p, d);

    # center, mirror it
    v += np.roll(v, d)
    v = np.roll(v, q - (((d * (q // p)) // 2) % q))
    return v / 2

# p->t rounding

def conv_comb_pt(v, q, p, t):
    v = conv_comb(v, q / t, 1)
    # p->t via q->p just puts it at an angle
    v = np.roll(v, q - (((q // t) + (q // p)) // 2))
    return v

def round5_fp_fast(dim, n, h, q, p, t, f, mu, B, add_coeff=False, verbose=False):

#   dif = mp_zeros(q)                   # full precision
#   dif[0] = mp.mpf(1.0)

    dif = np.zeros(q)                   # double precision approximation
    dif[0] = 1.0

    d = 0
    for w in h:                         # w = weight of +-d
        d += 1
        if add_coeff:
            w = 2 * w                   # closer to 2 * w - 1 if h ~= w/2
        w *= 2                          # double it here; can use floats
        if verbose:
            print "\nw" + str(d), "=", w

        for i in range(w):
            dif = conv_comb_qp(dif, q, p, d)
            if verbose:
                sys.stdout.write('.')
                sys.stdout.flush()
    if verbose:
        print "done"

    # convert to mpf if we were double
    dif = np.array([mp.mpf(dif[i]) for i in range(len(dif))])
    dif /= dif.sum()                    # normalize, just in case

    if verbose:
        print "eqp*R + eqp'*S".ljust(16) + str_stats(dif)

    dif = conv_comb_pt(dif, q, p, t)
    if verbose:
        print "rounded p->t".ljust(16) + str_stats(dif)

    # Decryption failure ranges for case when parameter B>1
    up_fail = (q + (1<<B)) // (1<<(B+1))
    down_fail = (q*( (1<<(B+1)) - 1) + (1<<B)) // (1<<(B+1))

    bfp = dif[up_fail:down_fail].sum();

    if verbose:
        print "per bit failure probability:", str_log2(bfp)
        print bfp

    ffp = mp.mpf(0);
    for j in range(f + 1, mu + 1):
        ffp += (mp.binomial(mu, j) * mp.power(bfp, j) *
                mp.power(mp.mpf(1) - bfp, (mu - j)))

    if verbose:
        print "more than", f, "errors in", mu, "bits:", str_log2(ffp)
        print ffp

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

# log_2(Cost) of __quantum__-enabled Enumeration attacks on Ideal/Euclidean lattice
# of dimension dim, using BKZ-2.0 with blocksize b, assuming only 1 call to the SVP oracle.
# Analogous to Q-core Enum + O(1) of "Estimate all the {LWE,NTRU} Schemes" paper.
# Used to estimate Enum attack-cost in NTRU-HRSS-KEM.
def svp_quantum_enumeration(b, ring=0):
    return float((0.18728*b*log(b, 2) - 1.0192*b + 16.1)/2.)


# Classical version of the above Enumeration attack. According to some opinions,
# integrating Grover into Enumeration is impractical,
# making the one below the only (as of now) practical Enumeration based attack.
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


#Albrecht sparse attack analysis.


#total running time complexity in bits, for Albrecht's "Drop A vectors corresponding to zero-components in the secret vector" attack (naive embodiment)
#time_bkz: Cost of running a lattice reduction with BKZ-b on the given lattice of dimension (n-k).
def albrecht_attack_cost_naive(n,k,h,time_bkz):
    repetitions = log2((binomial(n,k)) // (binomial((n-h),k)))    #Must be an integer. We use // (equivalent to int(floor())) to get a conservative (defender viewpoint) cost estimate of the attack.
    repetitions = max(0,repetitions)    #Ignore negative values
    #print "albrecht_attack_cost_naive(). repetitions",repetitions,"time_bkz",time_bkz
    return repetitions + int(floor(time_bkz))


#total running time complexity in bits, for Albrecht's "Drop sparse columns" attack (adaptive embodiment: tries to recover from the error of mistakenly dropping a non-sparse column)
def albrecht_attack_cost_adaptv(n,k,h,time_bkz):
    j = 0
    repetitions = 0
    l = k - 1            #Maximum number of erroneoulsy dropped non-sparse columns
    checks = 0
    while(j < l):
        drop_probab_j = (binomial((n-h),(k-j)) * binomial(h,j)) // (binomial(n,k))
        repetitions += drop_probab_j
        checks += binomial(k,j) * (2**j)
        j += 1
    if(repetitions>=log_infinity or checks>=log_infinity):
        return log_infinity + int(floor(time_bkz))
    else:
        #print "repetitions",repetitions,"checks",checks
        if(repetitions!=0):
            repetitions = checks // repetitions
        #else:    repetitions = 0
        repetitions = max(0,repetitions)    #ignore negative values
        return repetitions + int(floor(time_bkz))


#Optimize Albrecht's sparse column dropping attack for a given n and hammwt.
def optimize_albrecht_spattack(q,p,n,max_samples,theta,ring, use_lwr_reduction,is_hnf,is_binary,is_ternary,is_sparsehnf,is_spbinary,is_spternary,fastsearch,svp_cost_model):
    min_alb_cost=log_infinity
    hammwt=int(floor(theta*n))
    k_opt = 0
    attack_types=[1,2]    #we check for both the Primal and dual flavors of Albrecht's attack.
    for attack_type in attack_types:

        (b, m, cost_lr, dual_gap) = optimizeattack_b_m(q, p, n, max_samples, theta, is_hnf, is_binary, is_ternary, is_sparsehnf, is_spbinary, is_spternary, attack_type, ring, use_lwr_reduction, fastsearch, svp_cost_model)

        k_min = 0
        k_max = n-hammwt

        # Use tiny steps if not doing a fast search
        k_step = 1

        if fastsearch:
            k_step = (k_max - k_min) >> 4
        if (k_step < 1):
            k_step = 1


        while (k_step >= 1):


            for k in (k_min, k_max, k_step):    #The number of vectors in A that the attacker can choose to ignore.
                #First try using the Primal attack flavor. Usually stronger than the Dual attack flavor.

                cost_naive = albrecht_attack_cost_naive(n, k, hammwt, cost_lr)
                if (cost_naive < min_alb_cost):
                    min_alb_cost = cost_naive
                    k_opt = k

                cost_adaptv = albrecht_attack_cost_adaptv(n, k, hammwt, cost_lr)
                if (cost_adaptv < min_alb_cost):
                    min_alb_cost = cost_adaptv
                    k_opt = k

            # adjust ranges
            k_dist = (k_max - k_min + 4) >> 2
            k_min = max(k_min, k_opt - k_dist)
            k_max = min(k_max, k_opt + k_dist)
            k_step >>= 1

    return (min_alb_cost, k_opt)


#Load values of k (ignored components) and total attack cost from given files
def visualize_albrecht_attack_k_cost(k_filename,albcost_filename):
    k_file=open(k_filename,'r')
    albcost_file=open(albcost_filename,'r')
    k_list=[]
    albcost_list=[]
    for line in k_file:
        line=line.replace(","," ")
        line=line.split()
        if line:
            k_list=[int(i) for i in line]
    for line in albcost_file:
        line=line.replace(","," ")
        line=line.split()
        if line:
            albcost_list=[float(i) for i in line]
    k_file.close()
    albcost_file.close()
    return (k_list,albcost_list)


## End analysis of Albrecht's A-vector dropping attack against sparse-LWE, sparse-LWR ##


#Return the next largest prime from a list of prime nplus1's such that the nplus1-th cyclotomic polynomial \phi_nplus1(x) is irreducible modulo 2, and modulo 2**k for any k.
def next_cyclotomic_prime_nplus1_qpow2(nplus1_input):
    nplus1_list_pow2q = [2, 3, 5, 11, 13, 19, 29, 37, 53, 59, 61, 67, 83, 101, 107, 131, 139, 149, 163, 173, 179, 181, 197, 211, 227, 269, 293, 317, 347, 349, 373, 379, 389, 419, 421, 443, 461, 467, 491, 509, 523, 541, 547, 557, 563, 587, 613, 619, 653, 659, 661, 677, 701, 709, 757, 773, 787, 797, 821, 827, 829, 853, 859, 877, 883, 907, 941, 947, 1019, 1061, 1091, 1109, 1117, 1123, 1171, 1187, 1213, 1229, 1237, 1259, 1277, 1283, 1291, 1301, 1307, 1373, 1381, 1427, 1451, 1453, 1483, 1493, 1499, 1523, 1531, 1549, 1571, 1619, 1621, 1637, 1667, 1669, 1693, 1733, 1741, 1747, 1787, 1861, 1867, 1877, 1901, 1907, 1931, 1949, 1973, 1979, 1987, 1997, 2027, 2029, 2053, 2069, 2083, 2099, 2131, 2141, 2213, 2221, 2237, 2243, 2267, 2269, 2293, 2309, 2333, 2339, 2357, 2371, 2389, 2437, 2459, 2467, 2477, 2531, 2539, 2549, 2557, 2579, 2621, 2659, 2677, 2683, 2693, 2699, 2707, 2741, 2789, 2797, 2803, 2819, 2837, 2843, 2851, 2861, 2909, 2939, 2957, 2963]
    i = 0
    while (nplus1_list_pow2q[i]<=nplus1_input and (i+1)!=len(nplus1_list_pow2q)):
        i = i + 1
    return nplus1_list_pow2q[i]
