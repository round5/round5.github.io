############## Use Example ##############
#
# WARNING. This code is terribly slow for large parameters like n > 2**15 :
#           (It may takes few hours, or even one day(!))
#           To reproduce the result in the paper,
#           we recommend to try rather small parameters.
#
# > sage: from estimator_primal_hybrid import *
# > sage: n = 2048; q = 2**45; alpha = 8/q; h = 64;
# > sage: cost_model = BKZ.enum # enumeration cost model for SVP oracle
# > sage: # cost_model = BKZ.sieve # sieving cost model for SVP oracle
# > sage: primal_hyb(n, alpha, q, secret_distribution=((-1,1),h), reduction_cost_model=cost_model)
#
#    ***** Final Best rop : 123.928967 with (delta_0, r, m) = (1.008587, 40507, 24387)
#                        log(p_M, p_NP, p_s, p_c) = (-51.312136, 0.000000, -5.091901, -29.849079)
#
#
# Remark. 'mitm = 1' flag represents that the MitM hybrid attack cost estimation is applied.
#           Otherwise if 'mitm = 0,' that cost is estimated from the exhaustive hybrid attack.
#
# Remark. The output will be printed within command-line (in default setting) 
#           and a text file will be generated in name "n_logq_h.txt" in the same folder.
#
#       


from estimator import *
from sage.functions.error import Function_erf
import numpy as np
import logging
import sage.crypto.lwe

oo = PlusInfinity()

def my_binary_search(f, start, stop, param, predicate=lambda x, best: x<=best, *arg, **kwds):
    kwds[param] = stop
    D = {}
    D[stop] = f(*arg, **kwds)
    best = D[stop]
    b = ceil((start+stop)/2)
    direction = 0
    while True:
        if b not in D:
            kwds[param] = b
            D[b] = f(*arg, **kwds)
        if b == start:
            best = D[start]
            break
        if not predicate(D[b], best):
            if direction == 0:
                start = b
                b = ceil((stop+b)/2)
            else:
                stop = b
                b = floor((start+b)/2)
        else:
            best = D[b]
            if b-1 not in D:
                kwds[param] = b-1
                D[b-1] = f(*arg, **kwds)
            if predicate(D[b-1], best):
                stop = b
                b = floor((b+start)/2)
                direction = 0
            else:
                if b == stop:
                    break
                elif b+1 not in D:
                    kwds[param] = b+1
                    D[b+1] = f(*arg, **kwds)
                if not predicate(D[b+1], best):
                    break
                else:
                    start = b
                    b = ceil((stop+b)/2)
                    direction = 1
    return best

def prob_M(h_M, h, n, r):
    result = 0
    for i in range(h_M+1):
        result += binomial(h, 2*i) * binomial(n - h, r - 2*i)
    return result / binomial(n,r)
    
def sizeS(h_M, r):
    result = 0
    for i in range(h_M+1):
        result += 2**i * binomial(r, 2*i)
    return result

def prob_NP_and_GSA(m, n, alpha, q, r, delta_0, opt_GSA, h, h_M = False):
    dim = m + n-r + 1    
    scaling_factor = RR(sqrt(pi) / (2*alpha*q)) 
    if h_M is False:
        h_ = h * (n-r) / n + 1
        scale = RR(sqrt((n+1-r)/h_)*alpha*q/sqrt(2*pi))
    else:
        scale = RR(sqrt((n+1-r)/(h-2*h_M+1))*alpha*q/sqrt(2*pi))

    if opt_GSA is True:
        R = Mod_GSA(m, q, dim, delta_0, scale)
    else:
        R = GSA(m, q, dim, delta_0, scale)        
     
    probs = [RR(R[i] * scaling_factor).erf() for i in range(len(R))]
    return prod(probs), m, dim

def time_NP(dim):
    return dim / 2**(1.06)

def GSA(m, q, dim, delta_0, scale):
    det = RR((q**m * scale**(dim-m))**(1/dim))
    b = [delta_0**(dim-2*i) * det for i in range(dim)]
    return b

def Mod_GSA(m, q, dim, delta_0, scale):
    
    # Modified GSA for q-ary lattices 
    #       proposed by [Wun16]

    k = min(int(sqrt((dim - m)/log(delta_0, q/scale)).n()), dim)
    b_1 = [q for i in range(dim-k)]
    b_2 = [RR(delta_0**(k - 2*(i-(dim-k)-1)) * (q/scale)**((k-dim+m)/k) * scale) for i in range(dim - k, dim)]
    return b_1 + b_2

def cost_function_delta(n, q, alpha, h, delta_P, m_max,
                        reduction_cost_model, opt_GSA, verbose):
    
    line = 'Current delta_P = %f \n' % delta_P
    if verbose is True:
        print line
    
    kwds = {"n":n, "q":q, "alpha":alpha, "h":h, "delta_P":delta_P, "m_max":m_max,
            "reduction_cost_model": reduction_cost_model,
            "opt_GSA":opt_GSA, "verbose":verbose}

    best_delta = my_binary_search(cost_function_r, start = 0, stop = n, param="r", 
                                predicate=lambda x, best: RR(x["rop"])<=RR(best["rop"]), **kwds)
             
    line = '    * So far Best rop : %f with (delta_0, r, m) = (%f, %d, %d)\n' % (float(log(best_delta["rop"], 2)), best_delta["delta_0"], best_delta["r"], best_delta["m"]) 
    if verbose is True:   
        print line
    line = '                            log(p_M, p_NP) = (%f, %f)\n' % (RR(log(best_delta["p_M"],2)), RR(log(best_delta["p_NP"],2)))
    if verbose is True:   
        print line       

    return best_delta

def cost_function_r(n, q, alpha, h, r, delta_P, m_max,
                 reduction_cost_model, opt_GSA, verbose):
    
    kwds = {"n":n, "alpha":alpha, "q":q, "r":r, "delta_0":delta_P, "opt_GSA":opt_GSA, "h":h}
    
    p_NP, m, dim = my_binary_search(prob_NP_and_GSA, start = 50, stop = 2*n, param="m", 
                                predicate=lambda x, best: RR(x[0])>=RR(best[0]), **kwds)

    current = lattice_reduction_cost(reduction_cost_model, delta_P, dim, B=log(q, 2))
    
    p_M = 1.0
    

    total = 0
    for i in range(h+1):
        total += binomial(r, i) * 2**i

    ratio = RR(log(total, 2))

    if RR(p_NP).is_NaN() or n == r:
        current["rop"] = oo

    else:
        time_lat = current["rop"]
        time_post = 0
        S_i_sizes = [binomial(r, i) * 2**i for i in range(max(0, h-n+r), min(r, h)+1)]
        q_i = [binomial(n-r, h-i) * 2**(-i) / binomial(n, h) for i in range(max(0, h-n+r), min(r, h))]
        p_M = RR(0)
        sum_temp = RR(0)
        sizeS = 0
        while time_post < time_lat and len(q_i) > 0:
            index_i_current = np.argmax(q_i)
            size_S_i_current = S_i_sizes[index_i_current]
            del S_i_sizes[index_i_current]
            q_i_current = q_i[index_i_current]
            del q_i[index_i_current]
            sizeS += size_S_i_current
            p_M += RR(q_i_current * size_S_i_current)
            sum_temp += RR(pow(2, (2 * log(q_i_current, 2) / 3)) * size_S_i_current)
            time_post = RR(sqrt(sum_temp**3) / p_M ) * time_NP(dim)
            ratio = RR(log(total,2)) - RR(log(sizeS,2))
        
        probability = p_M * p_NP
        current["rop"] = time_lat + time_post
        current = current.repeat(1/probability, select={"m": False, "red": False})  
    
    current["ratio"] = ratio
    current["p_M"] = p_M
    current["p_NP"] = p_NP
    current["r"] = r
    current["m"] = m
    current = current.reorder(["rop"])
    return current

def primal_hyb(n, alpha, q, secret_distribution, 
               m = oo, h = None, success_probability=0.99, 
               reduction_cost_model=reduction_default_cost, 
               init_delta_P = None, opt_GSA = False, verbose=True):

    primald = partial(drop_and_solve, primal_usvp)
    kwds = {"postprocess":False}
    primalcost = primald(n, alpha, q, secret_distribution=secret_distribution, reduction_cost_model=reduction_cost_model, **kwds)
    
    duald = partial(drop_and_solve, dual_scale)
    dualcost = duald(n, alpha, q, secret_distribution=secret_distribution, reduction_cost_model=reduction_cost_model)

    if h is None:
        h = SDis.nonzero(secret_distribution, n)

    line = 'n = %d, logq = %d, stddev = %f, HW = %d\n' % (n, int(log(q,2)), float(stddevf(alpha*q)), h)
    print line    
    line = 'c.f. Primal cost = %5.1f, Dual cost = %5.1f\n' % (float(log(primalcost["rop"],2)), float(log(dualcost["rop"],2)))
    print line
    

    n, alpha, q, success_probability = Param.preprocess(n, alpha, q, success_probability)
    RR = parent(alpha)

    if not SDis.is_bounded_uniform(secret_distribution):
        raise NotImplementedError("Only bounded uniform secrets are currently supported.")

    m_max = m        
    
    if init_delta_P is None:
        init_delta_P = min(dualcost["delta_0"], primalcost["delta_0"])
        
    best = cost_function_delta(n, q, alpha, h, init_delta_P, m_max,
                               reduction_cost_model=reduction_cost_model, opt_GSA=opt_GSA, verbose=verbose)
    
    delta_step = min(0.005, (init_delta_P - 1)/2)

    while delta_step > 0.000005:
        
        current_pre = None
        if best["delta_0"] != init_delta_P:
            current_pre = cost_function_delta(n, q, alpha, h, best["delta_0"] - delta_step, m_max,
                                              reduction_cost_model=reduction_cost_model, 
                                              opt_GSA=opt_GSA, verbose=verbose)
            
        current_post = cost_function_delta(n, q, alpha, h, best["delta_0"] + delta_step, m_max,
                                          reduction_cost_model=reduction_cost_model, 
                                          opt_GSA=opt_GSA, verbose=verbose)

        delta_step /= 2
        min_cost = 0
        
        cost_pre = oo
        cost_old = RR(log(best["rop"],2))
        cost_post = RR(log(current_post["rop"],2))
        
        min_cost = min(cost_old, cost_post)
        if current_pre is not None:
            cost_pre = RR(log(current_pre["rop"],2))
            min_cost = min(min_cost, cost_pre)

        if min_cost == cost_post:
            best = current_post
        elif min_cost == cost_pre:
            best = current_pre
        
        if 0 < cost_old - min_cost < 0.5:
            break

        line = '** So far Best rop : %f with (d elta_0, r, m) = (%f, %d, %d)\n' % (float(log(best["rop"], 2)), best["delta_0"], best["r"], best["m"]) 
        if verbose is True:
            print line
        line = '                        log(p_M, p_NP) = (%f, %f)\n' % (RR(log(best["p_M"],2)), RR(log(best["p_NP"],2))) 
        if verbose is True:
            print line
        

    line = '***** Final Best rop : %f with (delta_0, r, m, log ratio) = (%f, %d, %d, %f)\n' % (float(log(best["rop"], 2)), best["delta_0"], best["r"], best["m"], best["ratio"]) 
    print line
    line = '                        log(p_M, p_NP) = (%f, %f)\n' % (RR(log(best["p_M"],2)), RR(log(best["p_NP"],2))) 
    print line

    return best


def getParams():
    from chosen_parameter_sets import getParameterSets
    r5paramsets = getParameterSets()
    params = []
    for i in range(len(r5paramsets)):
        ring_flag = True if r5paramsets[i].d==r5paramsets[i].n else False
        next = [r5paramsets[i].d, r5paramsets[i].q, r5paramsets[i].p, r5paramsets[i].h, ring_flag, r5paramsets[i].name]
        params.append(next)
    return params


f = open("Quantum_Sieve_Round5.txt", "w")
f.write("Security Estimations based on Quantum Sieving\n")
f.close()
#f = open("Quantum_Enum_Round5.txt", "w")
#f.write("Security Estimations based on Quantum Enumeration\n")
#f.close()

init_delta_P = None

Params = getParams()

for i in range(len(Params)):
    n = Params[i][0]; q = RR(Params[i][1]); p = RR(Params[i][2]); h = Params[i][3]; ring_flag = Params[i][4]; param_name = Params[i][5]
    stddev = RR(sqrt(((q/p)**2 - 1) / 12))
    alpha = RR(sqrt(2*pi) * stddev / q)

    cost_model = partial(BKZ.R5ringSieve, mode = "quantum") if ring_flag else partial(BKZ.R5nonringSieve, mode = "quantum")
    result = primal_hyb(n, alpha, q, secret_distribution=((-1,1),h), reduction_cost_model=cost_model, init_delta_P=init_delta_P, verbose=True)
    f = open("Quantum_Sieve_Round5.txt", "a")
    f.write(param_name + " / (n, q, p, h, ring) = (%d, %d, %d, %d, %s) / Bit - security = %f\n" % (n, q, p, h, ring_flag, float(log(result["rop"], 2))))
    f.close()

    #cost_model = partial(BKZ.R5Enum, mode = "quantum")
    #result = primal_hyb(n, alpha, q, secret_distribution=((-1,1),h), reduction_cost_model=cost_model, init_delta_P=init_delta_P, verbose=True)
    #f = open("Quantum_Enum_Round5.txt", "a")
    #f.write(param_name + " / (n, q, p, h, ring) = (%d, %d, %d, %d, %s) / Bit - security = %f\n" % (n, q, p, h, ring_flag, float(log(result["rop"], 2))))
    #f.close()