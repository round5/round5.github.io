#!/usr/bin/env python
from __future__ import division
import sys, os
import pickle

#
# Setup of analysis parameters
#

# Oracles considered



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
               expr_param=False):     # Experimental/test parameter set

#        # Sanity checks to see if name is valid and matches parameters
#        # Skipped, if the parameter is an experimental one
#        if not expr_param:
#            if (len(name) < 11 or name[0:2] != 'R5' or name[4] != '_' or (name[6:9] != 'PKE' and name[6:9] != 'KEM') or name[9] != '_' or (name[10] == 'x' and f != -1) or (name[10] != 'x' and name[10] != str(f))):
#                print "Got %s %d" % (name[10], f)
#                raise ValueError("Invalid name: %s" % name)
#            if (name[2:4] == 'ND'):
#                if (n != d):
#                    raise ValueError("n(%d) != d(%d) for ring: %s" % (n, d, name))
#            elif (name[2:4] == 'N1'):
#                if (n != 1):
#                    raise ValueError("n(%d) != 1 for non-ring: %s" % (n, name))
#            else:
#                raise ValueError("Invalid name: %s" % name)
#            if (name[5] != '0' and name[5] != '1' and name[5] != '3' and name[5] != '5'):
#                raise ValueError("Invalid level in name: %s" % name)
#
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
        self.fail = 0
#        self.primal = {}
#        self.dual = {}
#        self.hybrid = {}
#        self.albrecht = {}
#        self.guess_dual = {}
       