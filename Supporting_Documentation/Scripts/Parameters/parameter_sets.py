#!/usr/bin/env python
from __future__ import division
import sys
from analysis import hybrid_attack_cost_r_nist, optimizeattack_b_m, optimize_albrecht_spattack, svp_classical_core_sieving, svp_quantum_core_sieving, svp_classical_enumeration, svp_quantum_enumeration, build_sparse_ternary_law, round5_fp_fast

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
               name,
               d,                     # lattice dimension
               n,                     # ring dimension
               h,                     # hamming weight of secrets
               q_bits,                # log_2(GLWR large modulus)
               p_bits,                # log_2(GLWR rounding modulus)
               t_bits,                # log_2(GLWR compression modulus)
               b_bits,                # log_2(Bits extracted from ciphertext symbols)
               n_bar,                 # Alice's instances
               m_bar,                 # Bob's instances
               kappa_bytes,           # kappa (in bytes)
               f,                     # Bit errors corrected by XEf (-1 to indicate N(x) even if there is no error correction)
               xe,                    # Number of bits for error correction
               switchrings,           # Ciphertext computed modulo N(x)=x^(n+1)-1 or Phi(x).
               expr_param=False):     # Experimental/test parameter set

        # Sanity checks to see if name is valid and matches parameters
        # Skipped, if the parameter is an experimental one
        if not expr_param:
            if (len(name) < 11 or name[0:2] != 'R5' or name[4] != '_' or (name[6:9] != 'PKE' and name[6:9] != 'KEM') or name[9] != '_' or (name[10] == 'x' and f != -1) or (name[10] != 'x' and name[10] != str(f))):
                print "Got %s %d" % (name[10], f)
                raise ValueError("Invalid name: %s" % name)
            if (name[2:4] == 'ND'):
                if (n != d):
                    raise ValueError("n(%d) != d(%d) for ring: %s" % (n, d, name))
            elif (name[2:4] == 'N1'):
                if (n != 1):
                    raise ValueError("n(%d) != 1 for non-ring: %s" % (n, name))
            else:
                raise ValueError("Invalid name: %s" % name)
            if (name[5] != '0' and name[5] != '1' and name[5] != '3' and name[5] != '5'):
                raise ValueError("Invalid level in name: %s" % name)

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
        self.mu = ceildiv(self.kappa + xe, b_bits)
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

    def __run_albrechtattack(self, svp_oracle, fastsearch=True):
        """
        Summarize security level against Martin
        Albrecht's attack on sparse secret-keys
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

        (cost, colsdrop) = optimize_albrecht_spattack(self.q, self.p, self.d, max_m,
                                                          self.h/1./self.d,
                                                          ring, use_lwr_reduction,
                                                          # secret-key distro is always sparse-ternary
                                                          False,False,False,False,False,True,
                                                          fastsearch,svp_oracle)
        return {'cost': cost, 'columns': colsdrop}

    def analyze(self):
        # Failure rate
        self.fail = round5_fp_fast(self.d, self.n, [self.h],
                                       self.q, self.p, self.t,
                                       self.f, self.mu, self.b_bits,
                                       not self.switchrings)     # if ring switching is done, then
                                                                 # convolution lengths must be adapted

        # initialize the dictionaries for the values
        self.primal = {}
        self.dual = {}
        self.hybrid = {}
        self.albrecht = {}
        for svp_oracle in svp_oracles:
            oracleName = svp_oracle.__name__
            self.primal[oracleName] = self.__run_primal(svp_oracle)
            self.dual[oracleName] = self.__run_dual(svp_oracle)
            self.hybrid[oracleName] = self.__run_hybrid(svp_oracle)
            self.albrecht[oracleName] = self.__run_albrechtattack(svp_oracle)

    def __repr__(self):
        res = ""
        res += "%.80s\n" % ("------ %s Parameters -----------------------------------------------------------------------" % self.name)
        res += "d:            %d\n" % self.d
        res += "n:            %d\n" % self.n
        res += "h:            %d\n" % self.h
        res += "q_bits:       %d\n" % self.q_bits
        res += "p_bits:       %d\n" % self.p_bits
        res += "t_bits:       %d\n" % self.t_bits
        res += "q:            %d\n" % self.q
        res += "p:            %d\n" % self.p
        res += "t:            %d\n" % self.t
        res += "n_bar:        %d\n" % self.n_bar
        res += "m_bar:        %d\n" % self.m_bar
        res += "b_bits        %d\n" % self.b_bits
        res += "b             %d\n" % self.b
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
            res += "Albrecht's sparse-secrets attack (cost):               %d\n" % self.albrecht[oracleName]['cost']
            res += "Albrecht's sparse-secrets attack (cols A drop/ignore): %d\n" % self.albrecht[oracleName]['columns']

        return res

################################################################################


# Round5 ring parameter sets without error correction.
r5nd_1kem_0 = R5_paramSet("R5ND_1KEM_0d",618,618,104,11,8,4,1,1,1,16,0,0,False)
r5nd_3kem_0 = R5_paramSet("R5ND_3KEM_0d",786,786,384,13,9,4,1,1,1,24,0,0,False)
r5nd_5kem_0 = R5_paramSet("R5ND_5KEM_0d",1018,1018,428,14,9,4,1,1,1,32,0,0,False)

r5nd_1pke_0 = R5_paramSet("R5ND_1PKE_0d",586,586,182,13,9,4,1,1,1,16,0,0,False)
r5nd_3pke_0 = R5_paramSet("R5ND_3PKE_0d",852,852,212,12,9,5,1,1,1,24,0,0,False)
r5nd_5pke_0 = R5_paramSet("R5ND_5PKE_0d",1170,1170,222,13,9,5,1,1,1,32,0,0,False)

# Round5 ring parameter sets with error correction.
r5nd_1kem_5 = R5_paramSet("R5ND_1KEM_5d",490,490,162,10,7,3,1,1,1,16,5,190,True)
r5nd_3kem_5 = R5_paramSet("R5ND_3KEM_5d",756,756,242,12,8,2,1,1,1,24,5,218,True)
r5nd_5kem_5 = R5_paramSet("R5ND_5KEM_5d",940,940,414,12,8,2,1,1,1,32,5,234,True)

r5nd_1pke_5 = R5_paramSet("R5ND_1PKE_5d",508,508,136,10,7,4,1,1,1,16,5,190,True)
r5nd_3pke_5 = R5_paramSet("R5ND_3PKE_5d",756,756,242,12,8,3,1,1,1,24,5,218,True)
r5nd_5pke_5 = R5_paramSet("R5ND_5PKE_5d",946,946,388,11,8,5,1,1,1,32,5,234,True)

# Round5 non-ring parameter sets
# Note: setting the "switchrings" parameter=True for non-ring variants
# due to absence of correlation issue associated with reductions
# modulo the prime-order cyclotomic polynomial
r5n1_1kem_0 = R5_paramSet("R5N1_1KEM_0d",594,1,238,13,10,7,3,7,7,16,0,0,True)
r5n1_3kem_0 = R5_paramSet("R5N1_3KEM_0d",881,1,238,13,10,7,3,8,8,24,0,0,True)
r5n1_5kem_0 = R5_paramSet("R5N1_5KEM_0d",1186,1,712,15,12,7,4,8,8,32,0,0,True)

r5n1_1pke_0 = R5_paramSet("R5N1_1PKE_0d",636,1,114,12,9,6,2,8,8,16,0,0,True)
r5n1_3pke_0 = R5_paramSet("R5N1_3PKE_0d",876,1,446,15,11,7,3,8,8,24,0,0,True)
r5n1_5pke_0 = R5_paramSet("R5N1_5PKE_0d",1217,1,462,15,12,9,4,8,8,32,0,0,True)

# Round5 parameters for specific use cases
r5nd_0kem_2iot     = R5_paramSet("R5ND_0KEM_2iot",372,372,178,11,7,3,1,1,1,16,2,53,True)
r5nd_1kem_4longkey = R5_paramSet("R5ND_1KEM_4longkey",490,490,162,10,7,3,1,1,1,24,4,163,True)
r5n1_3pke_0smallct = R5_paramSet("R5N1_3PKE_0smallCT",757,1,378,14,9,4,1,192,1,24,0,0,True)

# Some other alternatives that are not submitted, but that could be interesting.
r5nd_5kem_5td   = R5_paramSet("R5ND_5KEM_5td",946,946,388,12,8,2,1,1,1,32,5,234,True)
r5nd_1pke_5td   = R5_paramSet("R5ND_1PKE_5td",522,522,124,10,7,4,1,1,1,16,5,190,True)
r5nd_5pke_5td   = R5_paramSet("R5ND_5PKE_5td",946,946,388,12,8,3,1,1,1,32,5,234,True)
r5nd_0kem_3iots = R5_paramSet("R5ND_0KEM_3iots",372,372,178,11,7,3,1,1,1,16,3,91,True)
r5nd_0pke_4iot  = R5_paramSet("R5ND_0PKE_4iot",418,418,92,11,7,3,1,1,1,16,4,149,True)

# Experimental "High-failure" Round5 parameters, to intentionally force decryption
# errors, and cross-check with statistically computed failure rates obtained from actual
# simulations of Round5 (see Correctness/r5sim.c)
r5nd_0kem_0fail_phi_0 = R5_paramSet("R5ND_0KEM_0fail_phi_0",800,800,170,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_1 = R5_paramSet("R5ND_0KEM_0fail_phi_1",800,800,180,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_2 = R5_paramSet("R5ND_0KEM_0fail_phi_2",800,800,200,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_3 = R5_paramSet("R5ND_0KEM_0fail_phi_3",800,800,220,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_4 = R5_paramSet("R5ND_0KEM_0fail_phi_4",800,800,250,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_5 = R5_paramSet("R5ND_0KEM_0fail_phi_5",800,800,270,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_6 = R5_paramSet("R5ND_0KEM_0fail_phi_6",800,800,300,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_7 = R5_paramSet("R5ND_0KEM_0fail_phi_7",800,800,320,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_8 = R5_paramSet("R5ND_0KEM_0fail_phi_8",800,800,350,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_9 = R5_paramSet("R5ND_0KEM_0fail_phi_9",800,800,370,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_10 = R5_paramSet("R5ND_0KEM_0fail_phi_10",800,800,400,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_11 = R5_paramSet("R5ND_0KEM_0fail_phi_11",800,800,420,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_12 = R5_paramSet("R5ND_0KEM_0fail_phi_12",800,800,440,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_13 = R5_paramSet("R5ND_0KEM_0fail_phi_13",800,800,450,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_14 = R5_paramSet("R5ND_0KEM_0fail_phi_14",800,800,470,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_15 = R5_paramSet("R5ND_0KEM_0fail_phi_15",800,800,500,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_16 = R5_paramSet("R5ND_0KEM_0fail_phi_16",800,800,520,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_17 = R5_paramSet("R5ND_0KEM_0fail_phi_17",800,800,540,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_18 = R5_paramSet("R5ND_0KEM_0fail_phi_18",800,800,550,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_19 = R5_paramSet("R5ND_0KEM_0fail_phi_19",800,800,570,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_20 = R5_paramSet("R5ND_0KEM_0fail_phi_20",800,800,590,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_21 = R5_paramSet("R5ND_0KEM_0fail_phi_21",800,800,600,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_22 = R5_paramSet("R5ND_0KEM_0fail_phi_22",800,800,620,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_23 = R5_paramSet("R5ND_0KEM_0fail_phi_23",800,800,640,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_24 = R5_paramSet("R5ND_0KEM_0fail_phi_24",800,800,650,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_25 = R5_paramSet("R5ND_0KEM_0fail_phi_25",800,800,670,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_26 = R5_paramSet("R5ND_0KEM_0fail_phi_26",800,800,700,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_27 = R5_paramSet("R5ND_0KEM_0fail_phi_27",800,800,720,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_28 = R5_paramSet("R5ND_0KEM_0fail_phi_28",800,800,740,11,7,4,1,1,1,16,0,0,False)
r5nd_0kem_0fail_phi_29 = R5_paramSet("R5ND_0KEM_0fail_phi_29",800,800,750,11,7,4,1,1,1,16,0,0,False)

r5nd_0kem_xfail_ntru_0 = R5_paramSet("R5ND_0KEM_xfail_ntru_0",800,800,170,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_1 = R5_paramSet("R5ND_0KEM_xfail_ntru_1",800,800,180,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_2 = R5_paramSet("R5ND_0KEM_xfail_ntru_2",800,800,200,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_3 = R5_paramSet("R5ND_0KEM_xfail_ntru_3",800,800,220,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_4 = R5_paramSet("R5ND_0KEM_xfail_ntru_4",800,800,250,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_5 = R5_paramSet("R5ND_0KEM_xfail_ntru_5",800,800,270,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_6 = R5_paramSet("R5ND_0KEM_xfail_ntru_6",800,800,300,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_7 = R5_paramSet("R5ND_0KEM_xfail_ntru_7",800,800,320,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_8 = R5_paramSet("R5ND_0KEM_xfail_ntru_8",800,800,350,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_9 = R5_paramSet("R5ND_0KEM_xfail_ntru_9",800,800,370,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_10 = R5_paramSet("R5ND_0KEM_xfail_ntru_10",800,800,400,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_11 = R5_paramSet("R5ND_0KEM_xfail_ntru_11",800,800,420,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_12 = R5_paramSet("R5ND_0KEM_xfail_ntru_12",800,800,440,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_13 = R5_paramSet("R5ND_0KEM_xfail_ntru_13",800,800,450,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_14 = R5_paramSet("R5ND_0KEM_xfail_ntru_14",800,800,470,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_15 = R5_paramSet("R5ND_0KEM_xfail_ntru_15",800,800,500,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_16 = R5_paramSet("R5ND_0KEM_xfail_ntru_16",800,800,520,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_17 = R5_paramSet("R5ND_0KEM_xfail_ntru_17",800,800,540,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_18 = R5_paramSet("R5ND_0KEM_xfail_ntru_18",800,800,550,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_19 = R5_paramSet("R5ND_0KEM_xfail_ntru_19",800,800,570,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_20 = R5_paramSet("R5ND_0KEM_xfail_ntru_20",800,800,590,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_21 = R5_paramSet("R5ND_0KEM_xfail_ntru_21",800,800,600,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_22 = R5_paramSet("R5ND_0KEM_xfail_ntru_22",800,800,620,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_23 = R5_paramSet("R5ND_0KEM_xfail_ntru_23",800,800,640,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_24 = R5_paramSet("R5ND_0KEM_xfail_ntru_24",800,800,650,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_25 = R5_paramSet("R5ND_0KEM_xfail_ntru_25",800,800,670,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_26 = R5_paramSet("R5ND_0KEM_xfail_ntru_26",800,800,700,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_27 = R5_paramSet("R5ND_0KEM_xfail_ntru_27",800,800,720,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_28 = R5_paramSet("R5ND_0KEM_xfail_ntru_28",800,800,740,11,7,4,1,1,1,16,-1,0,True)
r5nd_0kem_xfail_ntru_29 = R5_paramSet("R5ND_0KEM_xfail_ntru_29",800,800,750,11,7,4,1,1,1,16,-1,0,True)

#
# The parameter sets to use
#
r5paramSets = [
    r5nd_1kem_0,
    r5nd_3kem_0,
    r5nd_5kem_0,
    r5nd_1pke_0,
    r5nd_3pke_0,
    r5nd_5pke_0,

    r5nd_1kem_5,
    r5nd_3kem_5,
    r5nd_5kem_5,
    r5nd_1pke_5,
    r5nd_3pke_5,
    r5nd_5pke_5,

    r5n1_1kem_0,
    r5n1_3kem_0,
    r5n1_5kem_0,
    r5n1_1pke_0,
    r5n1_3pke_0,
    r5n1_5pke_0,

    # Specific use cases
    r5nd_0kem_2iot,
    r5nd_1kem_4longkey,
    r5n1_3pke_0smallct,

    # Alternatives
    # r5nd_0kem_3iots,
    # r5nd_0pke_4iot,
    # r5nd_5kem_5td,
    # r5nd_1pke_5td,
    # r5nd_5pke_5td,

    # High failure rate testing
    r5nd_0kem_0fail_phi_0,
    r5nd_0kem_0fail_phi_1,
    r5nd_0kem_0fail_phi_2,
    r5nd_0kem_0fail_phi_3,
    r5nd_0kem_0fail_phi_4,
    r5nd_0kem_0fail_phi_5,
    r5nd_0kem_0fail_phi_6,
    r5nd_0kem_0fail_phi_7,
    r5nd_0kem_0fail_phi_8,
    r5nd_0kem_0fail_phi_9,
    r5nd_0kem_0fail_phi_10,
    r5nd_0kem_0fail_phi_11,
    r5nd_0kem_0fail_phi_12,
    r5nd_0kem_0fail_phi_13,
    r5nd_0kem_0fail_phi_14,
    r5nd_0kem_0fail_phi_15,
    r5nd_0kem_0fail_phi_16,
    r5nd_0kem_0fail_phi_17,
    r5nd_0kem_0fail_phi_18,
    r5nd_0kem_0fail_phi_19,
    r5nd_0kem_0fail_phi_20,
    r5nd_0kem_0fail_phi_21,
    r5nd_0kem_0fail_phi_22,
    r5nd_0kem_0fail_phi_23,
    r5nd_0kem_0fail_phi_24,
    r5nd_0kem_0fail_phi_25,
    r5nd_0kem_0fail_phi_26,
    r5nd_0kem_0fail_phi_27,
    r5nd_0kem_0fail_phi_28,
    r5nd_0kem_0fail_phi_29,
    r5nd_0kem_xfail_ntru_0,
    r5nd_0kem_xfail_ntru_1,
    r5nd_0kem_xfail_ntru_2,
    r5nd_0kem_xfail_ntru_3,
    r5nd_0kem_xfail_ntru_4,
    r5nd_0kem_xfail_ntru_5,
    r5nd_0kem_xfail_ntru_6,
    r5nd_0kem_xfail_ntru_7,
    r5nd_0kem_xfail_ntru_8,
    r5nd_0kem_xfail_ntru_9,
    r5nd_0kem_xfail_ntru_10,
    r5nd_0kem_xfail_ntru_11,
    r5nd_0kem_xfail_ntru_12,
    r5nd_0kem_xfail_ntru_13,
    r5nd_0kem_xfail_ntru_14,
    r5nd_0kem_xfail_ntru_15,
    r5nd_0kem_xfail_ntru_16,
    r5nd_0kem_xfail_ntru_17,
    r5nd_0kem_xfail_ntru_18,
    r5nd_0kem_xfail_ntru_19,
    r5nd_0kem_xfail_ntru_20,
    r5nd_0kem_xfail_ntru_21,
    r5nd_0kem_xfail_ntru_22,
    r5nd_0kem_xfail_ntru_23,
    r5nd_0kem_xfail_ntru_24,
    r5nd_0kem_xfail_ntru_25,
    r5nd_0kem_xfail_ntru_26,
    r5nd_0kem_xfail_ntru_27,
    r5nd_0kem_xfail_ntru_28,
    r5nd_0kem_xfail_ntru_29,
    ]

def getParameterSets():
    return r5paramSets

def getSvpOracleNames():
    return [svp_oracle.__name__ for svp_oracle in svp_oracles]
