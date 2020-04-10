load quantum_hybrid.sage

sage.repl.load.load('r5_parameter_set.py', globals())
sage.repl.load.load('main_parameter_definitions.py', globals())
sage.repl.load.load('main_parameter_sets_ee.py', globals())
sage.repl.load.load('chosen_parameter_sets.py', globals())


def find_r_quantum_improved_optimizer(n, q, h, sigma, secret_structure, bkzmodel, nr_rotations, scaling, print_results, preci, verbose):

    min_cost = 999999999
    step = 256
    
    k_lowest = 0
    k_highest = 3*n/4 #n 
    k_low = k_lowest
    k_high = k_highest
    k_opt = k_lowest

    # ratio between the keys that the attacker tries and the total number of keys
    s_lowest = 0
    s_highest = 512
    s_low = s_lowest
    s_high = s_highest
    s_opt = s_lowest
    
    m_lowest = 0
    m_highest = n
    m_low = m_lowest
    m_high = m_highest
    m_opt = m_lowest
    
    while (step > 1) :
        for m in range (m_low, m_high, step):
            for k in range (k_low, k_high, step):
                for s in range (s_low, s_high, step):
		    #print n, q, m, h, sigma, secret_structure, bkzmodel, nr_rotations, scaling, 2**(-s), [k], print_results, preci
                    out = find_r_quantum_improved(n, q, m, h, sigma, secret_structure, bkzmodel, nr_rotations, scaling, 2**(-s), [k], print_results, preci, verbose)
                    if out[0] < min_cost:
                        k_opt = k
                        s_opt = s
                        m_opt = m
                        min_cost = out[0]

        k_low = max(k_lowest, k_opt - step)
        k_high = min(k_highest, k_opt + step)
        s_low = max(s_lowest, s_opt - step)
        s_high = min(s_highest, s_opt + step)
        m_low = max(m_lowest, m_opt - step)
        m_high = min(m_highest, m_opt + step)
        
        step = int(step / 2)
    
    print ("BEST FOUND PARAMETERS:", n, q, m_opt, h, sigma, secret_structure, bkzmodel, nr_rotations, scaling, 2**(-s_opt), [k_opt], print_results, preci, verbose)   
    print ("Results for [bitsec, opt_bs, opt_r, m, round(log(size_S, 2))]")
    return find_r_quantum_improved(n, q, m_opt, h, sigma, secret_structure, bkzmodel, nr_rotations, scaling, 2**(-s_opt), [k_opt], print_results, preci, verbose)


for paramSet in r5paramSets:
    print ("Starting computation of parameter set " + paramSet.name)
    d = paramSet.d
    q = paramSet.q
    p = paramSet.p
    h = paramSet.h

    verbose = False

    sigma = sqrt(1/12*((q*q)/(p*p) -1 ))
    secret_structure = 'hweight_trinary'
    bkzmodel = 'NIST_core_enum'
    nr_rotations = 1
    scaling = True
    print_results = False
    preci = 1
    print (d, q, p, h, sigma)

    print (find_r_quantum_improved_optimizer(d, q, h, sigma, secret_structure, bkzmodel, nr_rotations, scaling, print_results, preci, verbose))
    print ("The above is " + paramSet.name)







