#!/usr/bin/env python
from __future__ import division
import sys, os
import pickle
from analysis import hybrid_attack_cost_r_nist
from analysis import optimizeattack_b_m
from analysis import svp_classical_core_sieving
from analysis import svp_quantum_core_sieving
from analysis import svp_classical_enumeration
from analysis import svp_quantum_enumeration
from analysis import build_sparse_ternary_law
from analysis import round5_fp_fast
from analysis import optimize_guessing_dual
#from compute_maxH_ct import compute_hi

#
# Setup of analysis parameters
#

# Oracles considered
svp_oracles = [ svp_classical_core_sieving, svp_quantum_core_sieving, svp_classical_enumeration, svp_quantum_enumeration ]


#
# Convenience functions
#

def ceildiv(a, b):
    """
    Calculates the division of two integers, rounded up.
    """
    return ((a + b - 1) // b)

def bits2bytes(b):
    """
    Calculates the number of bytes required to represent the given number of bits.
    """
    return ((b+7)//8)

#
# A class representing a Round5 parameter set
#
class R5_paramSet(object):

    def __init__(self,
               name = "R5ND_1KEM_5d",
               d = 490,                     # lattice dimension
               n = 490,                     # ring dimension
               h = 162,                     # hamming weight of secrets
               q_bits = 10,                # log_2(GLWR large modulus)
               p_bits = 7,                # log_2(GLWR rounding modulus)
               t_bits = 3,                # log_2(GLWR compression modulus)
               b_bits = 1,                # log_2(Bits extracted from ciphertext symbols)
               n_bar = 1,                 # Alice's instances
               m_bar = 1,                 # Bob's instances
               kappa_bytes = 16,           # kappa (in bytes)
               f = 5,                     # Bit errors corrected by XEf (-1 to indicate N(x) even if there is no error correction)
               xe = 190,                    # Number of bits for error correction
               switchrings = True,           # Ciphertext computed modulo N(x)=x^(n+1)-1 or Phi(x).
               expr_param=False,           # Experimental/test parameter set
               repetition_code=1):         # usage of a repetition code

        
        
        # Initalize parameters
        # Defaults in case of
        # experimental parameter
        self.name = name
        if expr_param:
            encrypt = 0
            self.suffix = ""
            self.type = "ND"
            self.level = 0
            self.protocol = 'KEM'
        else:
            encrypt = (name[6:9] == 'PKE')
            self.suffix = name[11:]
            self.type = name[2:4]
            self.level = int(name[5])
            self.protocol = name[6:9]
        
        self.d = d
        self.n = n
        self.h = h
        #self.hmax = 9999
        self.q_bits = q_bits
        self.q = 1<<q_bits
        self.p_bits = p_bits
        self.p = 1<<p_bits
        self.t_bits = t_bits
        self.t = 1<<t_bits
        self.z = max(self.p, self.t*self.q/self.p)
        self.b_bits = b_bits
        self.b = 1<<b_bits
        self.n_bar = n_bar
        self.m_bar = m_bar
        self.f = f
        self.kappa = 8*kappa_bytes
        self.kappa_bytes = kappa_bytes
        self.xe = xe
        self.mu = ceildiv(repetition_code*(self.kappa + xe ), b_bits)
        
        self.k = d // n
        self.sk_size = self.kappa_bytes
        self.pk_size = self.kappa_bytes + bits2bytes(d * n_bar * p_bits)
        self.ct_size = bits2bytes(d * m_bar * p_bits) + bits2bytes(self.mu * t_bits)
        self.CRYPTO_SECRETKEYBYTES = self.sk_size + self.kappa_bytes + self.pk_size if encrypt else self.sk_size
        self.CRYPTO_PUBLICKEYBYTES = self.pk_size
        self.CRYPTO_BYTES = self.ct_size + self.kappa_bytes + 16 if encrypt else self.kappa_bytes
        self.CRYPTO_CIPHERTEXTBYTES = 0 if encrypt else self.ct_size
        self.bandwidth = self.CRYPTO_PUBLICKEYBYTES + (self.CRYPTO_BYTES if encrypt else self.CRYPTO_CIPHERTEXTBYTES)
        self.switchrings = switchrings
        self.fail = 0
        self.repetition_code = repetition_code
 
        self.primal = {'svp_classical_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_quantum_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_classical_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_quantum_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0}}
        self.dual = {'svp_classical_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_quantum_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_classical_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_quantum_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0}}
        self.hybrid = {'svp_classical_core_sieving':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_quantum_core_sieving':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_classical_enumeration':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_quantum_enumeration':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0}}
        self.guess_dual = {'svp_classical_core_sieving':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_quantum_core_sieving':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_classical_enumeration':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_quantum_enumeration':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0}}
        

    def __run_primal(self, svp_oracle, fastsearch=True):
        """
        Summarize security level against the
        Primal lattice-reduction based attack,
        rescaled to exploit sparse-ternary secret keys.
        """

        # Max. number of GLWR samples
        # available to an attacker
        max_m = self.d
        if(self.d==self.n):
            # The ring case
            ring = 1
            max_m += self.mu
        else:
            ring = 0
            # the non ring case
            max_m += max(self.n_bar, self.m_bar)

        use_lwr_reduction = False    # Not using a LWE->LWR reduction to assess concrete security

        (b_opt,m_opt,cost,gap) = optimizeattack_b_m(self.q, self.p, self.d, max_m,
                                                        self.h/1./self.d,
                                                        # secret-key distro is always sparse-ternary
                                                        False,False,False,False,False,True,
                                                        1,                   # attack type is primal
                                                        ring, use_lwr_reduction,
                                                        fastsearch, svp_oracle)
        return {'cost': cost, 'bkz': b_opt, 'samples': self.d + m_opt}

    def __run_dual(self, svp_oracle, fastsearch=True):
        """
        Summarize security level against the
        Dual lattice-reduction based attack,
        rescaled to exploit sparse-ternary secret keys.
        """

        # Max. number of GLWR samples
        # available to an attacker
        max_m = self.d
        if(self.d==self.n):
            # The ring case
            ring = 1
            max_m += self.mu
        else:
            ring = 0
            # the non ring case
            max_m += max(self.n_bar, self.m_bar)

        use_lwr_reduction = False   # Not using a LWE->LWR reduction to assess concrete security

        (b_opt,m_opt,cost,gap) = optimizeattack_b_m(self.q, self.p, self.d, max_m,
                                                        self.h/1./self.d,
                                                        # secret-key distro is always sparse-ternary
                                                        False,False,False,False,False,True,
                                                        2,                   # attack type is dual
                                                        ring, use_lwr_reduction,
                                                        fastsearch, svp_oracle)
        return {'cost': cost, 'bkz': b_opt, 'samples': self.d + m_opt}

    def __run_hybrid(self, svp_oracle):
        """
        Summarize security level against the Hybrid
        lattice reduction and Meet-in-the-Middle attack
        """
        # Max. number of GLWR samples
        # available to an attacker
        max_m = self.d
        if(self.d==self.n):
            # The ring case
            ring = 1
            max_m += self.mu
        else:
            ring = 0
            # the non ring case
            max_m += max(self.n_bar, self.m_bar)

        # Secret-key distro
        secretkey_distro = build_sparse_ternary_law
        sparseness = self.h/1./self.d

        (mitm_cost, lr_cost, mitm_dim, blocksz_opt) = hybrid_attack_cost_r_nist(self.q, self. p, self.d, max_m,
                                                                       sparseness, ring,
                                                                       # secret-key distro is always sparse-ternary
                                                                       False,False,False,False,False,True,
                                                                       svp_oracle, secretkey_distro)
        return {'cost': lr_cost, 'mitm_cost': mitm_cost, 'mitm_dim': mitm_dim, 'bkz': blocksz_opt}

    def __run_guessing_dual(self,svp_oracle):
        # Summarize security level against guess_dual form of Albrecht's dual attack
        max_m = self.d
        if (self.d == self.n):
            # The ring case
            ring = 1
            max_m += self.mu
        else:
            ring = 0
            # the non ring case
            max_m += max(self.n_bar, self.m_bar)
        
        (mincost,bopt,kopt,lopt,mopt,log_guess_cost,bkzcost,ropt,cost_of_svp) = optimize_guessing_dual(60,self.q,self.p,self.d,max_m,
                                                                                                                       self.h/1./self.d,False,False,False,False,
                                                                                                                       False,True,
                                                                                                                       ring,False,svp_oracle, EDA_BKZ_LLL=False)
        return{'cost': mincost, 'columns dropped': kopt, 'Max weight': lopt, 'bkz':bkzcost, 'Guess cost': log_guess_cost, 'repetitions': ropt, 'BKZ rep': bkzcost-cost_of_svp}

    
    # this function determines whether a parameter set fulfils a given classical security level and failure probability assuming a sieving cost model.
    # if the conditions are not met at any point, then the computation is interrumpted.
    def analyze_with_constrains(self, c_sec_threshold, fail_threshold, margin, run_guess_dual):
        # Failure rate
        self.fail = round5_fp_fast(self.d, self.h, self.q, self.p, self.t, self.f, self.mu,  self.b_bits,self.repetition_code, not self.switchrings,)
        #print self.fail
        if self.fail > fail_threshold:
            return False
        self.primal = {'svp_classical_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_quantum_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_classical_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_quantum_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0}}
        self.dual = {'svp_classical_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_quantum_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_classical_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_quantum_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0}}
        self.hybrid = {'svp_classical_core_sieving':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_quantum_core_sieving':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_classical_enumeration':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_quantum_enumeration':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0}}
        self.guess_dual = {'svp_classical_core_sieving':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_quantum_core_sieving':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_classical_enumeration':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_quantum_enumeration':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0}}

        oracleName = 'svp_classical_core_sieving' #svp_oracle[0].__name__ # classical sieving
        self.primal[oracleName] = self.__run_primal(svp_oracles[0])
        #print self.primal[oracleName]
        if self.primal[oracleName]['cost'] < c_sec_threshold + margin:
            return False
        self.dual[oracleName] = self.__run_dual(svp_oracles[0])
        if self.dual[oracleName]['cost'] < c_sec_threshold + margin:
            return False
        self.hybrid[oracleName] = self.__run_hybrid(svp_oracles[0])
        if self.hybrid[oracleName]['cost'] < c_sec_threshold + margin:
            return False
        # only run if needed, otherwise, too slow
        if run_guess_dual:
            self.guess_dual[oracleName] = self.__run_guessing_dual(svp_oracles[0])
            if self.guess_dual[oracleName]['cost'] < c_sec_threshold: # do not use margin in this case.
                return False
        return True
  
    def analyze(self):
        # number of times to call secret key generation algorithm to ensure at least h different positions are filled in with probabiliy 1 - 2^{-kappa}
        #self.hmax = compute_hi(self.d, self.h, self.name, self.kappa)
        
        # Failure rate
        self.fail = round5_fp_fast(self.d, self.h, self.q, self.p, self.t, self.f, self.mu, self.b_bits, self.repetition_code, not self.switchrings, )
    
                                                            
        # initialize the dictionaries for the values
        self.primal = {'svp_classical_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_quantum_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_classical_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0},
                       'svp_quantum_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0}}
        self.dual = {'svp_classical_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_quantum_core_sieving':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_classical_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0},
                     'svp_quantum_enumeration':{'cost': 0, 'bkz': 0, 'samples': 0}}
        self.hybrid = {'svp_classical_core_sieving':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_quantum_core_sieving':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_classical_enumeration':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0},
                       'svp_quantum_enumeration':{'cost': 0, 'mitm_cost': 0, 'mitm_dim': 0, 'bkz': 0}}
        self.guess_dual = {'svp_classical_core_sieving':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_quantum_core_sieving':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_classical_enumeration':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0},
                           'svp_quantum_enumeration':{'cost': 0, 'columns dropped': 0, 'Max weight': 0, 'bkz':0, 'Guess cost': 0, 'repetitions': 0, 'BKZ rep': 0}}
        

        for svp_oracle in svp_oracles:
            oracleName = svp_oracle.__name__
            self.primal[oracleName] = self.__run_primal(svp_oracle)
            self.dual[oracleName] = self.__run_dual(svp_oracle)
            self.hybrid[oracleName] = self.__run_hybrid(svp_oracle)
            self.guess_dual[oracleName] = self.__run_guessing_dual(svp_oracle)


    # this function either loads from a file or computes a parameter set.
    def get(self, path = '.ouf/default/results/'):

        print( path)
        if not (os.path.exists(path)):
            try:
                os.makedirs(path)
            except OSError:
                pass
    
        total_path = path + str(self.name)
        
        if os.path.isfile( total_path ):
            print( "---- Reading analysis for %s" % self.name)
            f = open(total_path,"r")
            for svp_oracle in svp_oracles:
                oracleName = svp_oracle.__name__
                self.primal[oracleName]['cost'] = pickle.load(f)
                self.primal[oracleName]['bkz'] = pickle.load(f)
                self.primal[oracleName]['samples'] = pickle.load(f)

                self.dual[oracleName]['cost'] = pickle.load(f)
                self.dual[oracleName]['bkz'] = pickle.load(f)
                self.dual[oracleName]['samples'] = pickle.load(f)

                self.hybrid[oracleName]['cost'] = pickle.load(f)
                self.hybrid[oracleName]['mitm_cost'] = pickle.load(f)
                self.hybrid[oracleName]['mitm_dim'] = pickle.load(f)
                self.hybrid[oracleName]['bkz'] = pickle.load(f)

                self.guess_dual[oracleName]['cost'] = pickle.load(f)
                self.guess_dual[oracleName]['columns_dropped'] = pickle.load(f)
                self.guess_dual[oracleName]['Max_weight'] = pickle.load(f)
                self.guess_dual[oracleName]['bkz'] = pickle.load(f)
                self.guess_dual[oracleName]['Guess cost'] = pickle.load(f)
                self.guess_dual[oracleName]['BKZ rep'] = pickle.load(f)
                self.guess_dual[oracleName]['repetitions'] = pickle.load(f)
            self.fail = pickle.load(f)
            self = pickle.load(f)
            f.close()
        else:
            print( "---- Running analysis for %s" % self.name)
            f = open(total_path,"w")
            self.analyze()
            
            for svp_oracle in svp_oracles:
                oracleName = svp_oracle.__name__
                pickle.dump(self.primal[oracleName]['cost'], f)
                pickle.dump(self.primal[oracleName]['bkz'], f)
                pickle.dump(self.primal[oracleName]['samples'], f)

                pickle.dump(self.dual[oracleName]['cost'], f)
                pickle.dump(self.dual[oracleName]['bkz'], f)
                pickle.dump(self.dual[oracleName]['samples'], f)

                pickle.dump(self.hybrid[oracleName]['cost'], f)
                pickle.dump(self.hybrid[oracleName]['mitm_cost'], f)
                pickle.dump(self.hybrid[oracleName]['mitm_dim'], f)
                pickle.dump(self.hybrid[oracleName]['bkz'], f)

                pickle.dump(self.guess_dual[oracleName]['cost'], f)
                pickle.dump(self.guess_dual[oracleName]['columns dropped'], f)
                pickle.dump(self.guess_dual[oracleName]['Max weight'], f)
                pickle.dump(self.guess_dual[oracleName]['bkz'], f)
                pickle.dump(self.guess_dual[oracleName]['Guess cost'], f)
                pickle.dump(self.guess_dual[oracleName]['BKZ rep'], f)
                pickle.dump(self.guess_dual[oracleName]['repetitions'], f)
            pickle.dump(self.fail, f)
            pickle.dump(self, f)
            f.close()

 # this function either loads from a file or computes a parameter set, and then stores it again.
    def reload(self, path = '.ouf/default/results/'):

        if not (os.path.exists(path)):
            try:
                os.makedirs(path)
            except OSError:
                pass
    
        total_path = path + str(self.name)
        
        if os.path.isfile( total_path ):
            print( "---- Reading analysis for %s" % self.name)
            f = open(total_path,"r")
            for svp_oracle in svp_oracles:
                oracleName = svp_oracle.__name__
                self.primal[oracleName]['cost'] = pickle.load(f)
                self.primal[oracleName]['bkz'] = pickle.load(f)
                self.primal[oracleName]['samples'] = pickle.load(f)

                self.dual[oracleName]['cost'] = pickle.load(f)
                self.dual[oracleName]['bkz'] = pickle.load(f)
                self.dual[oracleName]['samples'] = pickle.load(f)

                self.hybrid[oracleName]['cost'] = pickle.load(f)
                self.hybrid[oracleName]['mitm_cost'] = pickle.load(f)
                self.hybrid[oracleName]['mitm_dim'] = pickle.load(f)
                self.hybrid[oracleName]['bkz'] = pickle.load(f)

                self.guess_dual[oracleName]['cost'] = pickle.load(f)
                self.guess_dual[oracleName]['columns_dropped'] = pickle.load(f)
                self.guess_dual[oracleName]['Max_weight'] = pickle.load(f)
                self.guess_dual[oracleName]['bkz'] = pickle.load(f)
                self.guess_dual[oracleName]['Guess cost'] = pickle.load(f)
                self.guess_dual[oracleName]['BKZ rep'] = pickle.load(f)
                self.guess_dual[oracleName]['repetitions'] = pickle.load(f)
            self.fail = pickle.load(f)
            self = pickle.load(f)

            f.close()
        
            print( "---- Running analysis for %s" % self.name)
            f = open(total_path,"w")
            self.analyze_with_constrains(16, -300, 0, False)
            
            for svp_oracle in svp_oracles:
                oracleName = svp_oracle.__name__
                pickle.dump(self.primal[oracleName]['cost'], f)
                pickle.dump(self.primal[oracleName]['bkz'], f)
                pickle.dump(self.primal[oracleName]['samples'], f)

                pickle.dump(self.dual[oracleName]['cost'], f)
                pickle.dump(self.dual[oracleName]['bkz'], f)
                pickle.dump(self.dual[oracleName]['samples'], f)

                pickle.dump(self.hybrid[oracleName]['cost'], f)
                pickle.dump(self.hybrid[oracleName]['mitm_cost'], f)
                pickle.dump(self.hybrid[oracleName]['mitm_dim'], f)
                pickle.dump(self.hybrid[oracleName]['bkz'], f)

                pickle.dump(self.guess_dual[oracleName]['cost'], f)
                pickle.dump(self.guess_dual[oracleName]['columns dropped'], f)
                pickle.dump(self.guess_dual[oracleName]['Max weight'], f)
                pickle.dump(self.guess_dual[oracleName]['bkz'], f)
                pickle.dump(self.guess_dual[oracleName]['Guess cost'], f)
                pickle.dump(self.guess_dual[oracleName]['BKZ rep'], f)
                pickle.dump(self.guess_dual[oracleName]['repetitions'], f)
            pickle.dump(self, f)
            pickle.dump(self.fail, f)

            f.close()

    # this function either loads from a file or computes a parameter set,
    # exchanges hybrid values by the ones of Yongha, and then stores it again.
    def reload_hybrid(self, qshy, qehy, cshy, cehy, path = '.ouf/default/results/'):

        if not (os.path.exists(path + "hybrid_reloaded/")):
            try:
                os.makedirs(path + "hybrid_reloaded/")
            except OSError:
                pass

        reloaded_path = path + "hybrid_reloaded/" + str(self.name)

        print( "---- Storing for %s" % self.name)
        f = open(reloaded_path,"w")

        for svp_oracle in svp_oracles:
            oracleName = svp_oracle.__name__
            pickle.dump(self.primal[oracleName]['cost'], f)
            pickle.dump(self.primal[oracleName]['bkz'], f)
            pickle.dump(self.primal[oracleName]['samples'], f)

            pickle.dump(self.dual[oracleName]['cost'], f)
            pickle.dump(self.dual[oracleName]['bkz'], f)
            pickle.dump(self.dual[oracleName]['samples'], f)

            if oracleName == 'svp_classical_core_sieving':
                #pickle.dump(self.hybrid[oracleName]['cost'], f)
                pickle.dump(cshy, f)
            elif oracleName == 'svp_quantum_core_sieving':
                #pickle.dump(self.hybrid[oracleName]['cost'], f)
                pickle.dump(qshy, f)
            elif oracleName == 'svp_classical_enumeration':
                #pickle.dump(self.hybrid[oracleName]['cost'], f)
                pickle.dump(cehy, f)
            elif oracleName == 'svp_quantum_enumeration':
                #pickle.dump(self.hybrid[oracleName]['cost'], f)
                pickle.dump(qehy, f)
            else:
                pickle.dump(-1, f)

            
            pickle.dump(self.hybrid[oracleName]['mitm_cost'], f)
            pickle.dump(self.hybrid[oracleName]['mitm_dim'], f)
            pickle.dump(self.hybrid[oracleName]['bkz'], f)

            pickle.dump(self.guess_dual[oracleName]['cost'], f)
            pickle.dump(self.guess_dual[oracleName]['columns dropped'], f)
            pickle.dump(self.guess_dual[oracleName]['Max weight'], f)
            pickle.dump(self.guess_dual[oracleName]['bkz'], f)
            pickle.dump(self.guess_dual[oracleName]['Guess cost'], f)
            pickle.dump(self.guess_dual[oracleName]['BKZ rep'], f)
            pickle.dump(self.guess_dual[oracleName]['repetitions'], f)

        pickle.dump(self.fail, f)
        pickle.dump(self, f)

        f.close()

    def print_summary(self):
        res = " "
        res += str(self.name)
        res += " d:"
        res += str(self.d)
        res += " h:"
        res += str(self.h)
#        res += " h_max:"
#        res += str(self.hmax)
        res += " d*h:"
        res += str(self.d*self.h)
        res += " bw:"
        res += str(self.bandwidth)
        res += " ct:"
        res += str(self.kappa_bytes + self.m_bar*(self.p_bits*self.d + self.t_bits * (self.xe + self.kappa_bytes * 8))/ 8 )
        res += " fp:"
        res += str(self.fail)
        res += " primal:"
        res += str(self.primal['svp_classical_core_sieving']['cost'])
        res += " dual:"
        res += str(self.dual['svp_classical_core_sieving']['cost'])
        res += " hybrid:"
        res += str(self.hybrid['svp_classical_core_sieving']['cost'])
        res += " guess+dual:"
        res += str(self.guess_dual['svp_classical_core_sieving']['cost'])
        #res += "\n"
        print( res)
        return res
    
    
    # this function stores the best parameter found for the current search in a file
    def store(self, path='./.out/default/completeSearch/', name='summary_search'):
        
        if not (os.path.exists(path)):
            try:
                os.makedirs(path)
            except OSError:
                pass
        
        f_name = path + name + '.py'

        if not (os.path.isfile(f_name)):
            f = open(f_name,"w")
            f.write("# generated file \n")
            f.write("import sys, inspect, os \n")
            f.write("currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) \n")
            f.write("parentdir = os.path.dirname(currentdir) \n")
            f.write("sys.path.insert(0,parentdir) \n")
            f.write("from r5_parameter_set import R5_paramSet \n")
            f.close()
        f = open(f_name,"a+")
        s = ''
        s+= " R5_paramSet(\""
        s+= self.name
        s+= '\",' + str(self.d)
        s+= ',' + str(self.n)
        s+= ',' + str(self.h)
        s+= ',' + str(self.q_bits)
        s+= ',' + str(self.p_bits)
        s+= ',' + str(self.t_bits)
        s+= ',' + str(self.b_bits)
        s+= ',' + str(self.n_bar)
        s+= ',' + str(self.m_bar)
        s+= ',' + str(self.kappa_bytes)
        s+= ',' + str(self.f)
        s+= ',' + str(self.xe)
        s+= ',' + str(self.switchrings)
        s+= ') \n'
        f.write(s)
        s = '#'
        s+= self.print_summary()
        f.write(s)
        f.close()

    def print_r5_parameter_set_c(self):
    
        res = ""
        res += "#elif defined("
        res +=  self.name
        res += ") \n"
        res += "#define PARAMS_KAPPA_BYTES "
        res += str(self.kappa_bytes)
        res += "\n"
        res += "#define PARAMS_D          "
        res += str(self.d)
        res += "\n"
        res += "#define PARAMS_N          "
        res += str(self.n)
        res += "\n"
        res += "#define PARAMS_H          "
        res += str(self.h)
#        res += "#define PARAMS_H_MAX          "
#        res += str(self.hmax)
        res += "\n"
        res += "#define PARAMS_Q_BITS     "
        res += str(self.q_bits)
        res += "\n"
        res += "#define PARAMS_P_BITS     "
        res += str(self.p_bits)
        res += "\n"
        res += "#define PARAMS_T_BITS     "
        res += str(self.t_bits)
        res += "\n"
        res += "#define PARAMS_B_BITS     "
        res += str(self.b_bits)
        res += "\n"
        res += "#define PARAMS_N_BAR      "
        res += str(self.n_bar)
        res += "\n"
        res += "#define PARAMS_M_BAR      "
        res += str(self.m_bar)
        res += "\n"
        res += "#define PARAMS_F          "
        res += str(self.f)
        res += "\n"
        res += "#define PARAMS_XE         "
        res += str(self.xe)
        res += "\n"
        res += "#define CRYPTO_ALGNAME    \""
        res += self.name
        res += "\" \n"
        return res
    
    def __repr__(self):
        res = ""
        res += "%.80s\n" % ("------ %s Parameters -----------------------------------------------------------------------" % self.name)
        res += "d:            %d\n" % self.d
        res += "n:            %d\n" % self.n
        res += "h:            %d\n" % self.h
#        res += "hmax:         %d\n" % self.hmax
        res += "q_bits:       %d\n" % self.q_bits
        res += "p_bits:       %d\n" % self.p_bits
        res += "t_bits:       %d\n" % self.t_bits
        res += "b_bits        %d\n" % self.b_bits
        res += "q:            %d\n" % self.q
        res += "p:            %d\n" % self.p
        res += "t:            %d\n" % self.t
        res += "b             %d\n" % self.b
        res += "n_bar:        %d\n" % self.n_bar
        res += "m_bar:        %d\n" % self.m_bar
        res += "f:            %d\n" % self.f
        res += "xe:           %d\n" % self.xe
        res += "mu:           %d\n" % self.mu
        res += "kappa:        %d\n" % self.kappa
        res += "kappa_bytes:  %d\n" % self.kappa_bytes
        res += "RingSwitched: %s\n" % (self.switchrings if (self.n == self.d) else "N/A")
        res += "\n"
        res += "CRYPTO_SECRETKEYBYTES:  %d\n" % self.CRYPTO_SECRETKEYBYTES
        res += "CRYPTO_PUBLICKEYBYTES:  %d\n" % self.CRYPTO_PUBLICKEYBYTES
        res += "CRYPTO_BYTES:           %d\n" % self.CRYPTO_BYTES
        res += "CRYPTO_CIPHERTEXTBYTES: %d\n" % self.CRYPTO_CIPHERTEXTBYTES
        res += "\n"
        res += "Bandwidth:              %d\n" % self.bandwidth
        res += "Failure-rate less than: 2^(%f)\n" % self.fail
        res += "\n"
        for svp_oracle in svp_oracles:
            oracleName = svp_oracle.__name__
            res += "%.80s\n" % ("------ %s Analysis (%s) -----------------------------------------------------------------------" % (self.name, oracleName))
            res += "Primal attack (cost):                                  %.3f\n" % self.primal[oracleName]['cost']
            res += "Primal attack (blocksize BKZ):                         %d\n" % self.primal[oracleName]['bkz']
            res += "Primal attack (GLWR samples):                          %d\n" % self.primal[oracleName]['samples']
            res += "\n"
            res += "Dual attack (cost):                                    %.3f\n" % self.dual[oracleName]['cost']
            res += "Dual attack (blocksize BKZ):                           %d\n" % self.dual[oracleName]['bkz']
            res += "Dual attack (GLWR samples):                            %d\n" % self.dual[oracleName]['samples']
            res += "\n"
            res += "Hybrid attack cost (lattice-reduction):                %.3f\n" % self.hybrid[oracleName]['cost']
            res += "Hybrid attack (meet in the middle cost):               %.3f\n" % self.hybrid[oracleName]['mitm_cost']
            res += "Hybrid attack (meet in the middle dimension):          %.3f\n" % self.hybrid[oracleName]['mitm_dim']
            res += "Hybrid attack (optimal BKZ blocksize):                 %d\n" % self.hybrid[oracleName]['bkz']
            res += "\n"
            res += "Albrecht's guess_dual attack (cost):                     %.3f\n" % self.guess_dual[oracleName]['cost']
            res += "Albrecht's guess_dual attack BKZ cost                    %.3f\n" % self.guess_dual[oracleName]['bkz']
            res += "Albrecht's guess_dual attack: log(no of repetitions)     %.3f\n" % self.guess_dual[oracleName]['repetitions']
        return res


def getSvpOracleNames():
    return [svp_oracle.__name__ for svp_oracle in svp_oracles]
