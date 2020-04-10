Round5 Parameter Summary
========================

This folder contains scripts used to obtain security, failure probability estimates, and certain Round5 parameters.

The main script is in the `main` folder. The `additional` folder includes some additional scripts.

1. Main Scripts
==============


The script `summarize.py` can be used to summarize the properties of all parameters proposed as part of the Round5 cryptosystem including failure probability estimates as well as primal, dual, hybrid and guessing+dual attacks.  

To obtain this summary, open terminal and executte:

``  python summarize.py `` 

Note that when running above script, the security estimates for hybrid are still as submitted to the second round. The Round5 team has incorporated a more recent analysis that is located in the `additional` folder. 

In total, **18 main parameters** are summarized, along with 3 additional “specific use-case”
parameters that demonstrate the flexibility of Round5 and applicability to 3
different  specialized usage scenarios.

The parameters themselves are summarized, along with security levels and 
performance details, including bandwidth requirements and decryption failure
rates.

Note that the script also includes other parameter sets of interest. 

Note that the parameters in this README are named according to the update in February 2020, i.e., as:

				`R5N*_*{CPA/CCA}_ *d`
				
but the parameter sets in the scripts are still named as:


				`R5N*_*{KEM/PKE}_ *d`

Main Round5 parameters
----------------------

* `R5ND_{1,3,5}CPA_0d`:

    Parameters for the ring variant of INDCPA-secure Round5 key-encapsulation 
    mechanism,
    for NIST security levels 1, 3 and 5. No forward error correction is used,
    hence the `0` at the end of the parameter designator. All polynomial
    multiplications are done modulo the prime-order cyclotomic polynomial. 

* `R5ND_{1,3,5}CCA_0d`:

   Parameters for the ring variant of INDCCA-secure Round5 key-encapsulation mechanism and public-key encryption,
   for NIST security levels 1, 3 and 5. No forward error correction is used,
   hence the `0` at the end of the parameter designator. All polynomial
   multiplications are done modulo the prime-order cyclotomic polynomial. 

* `R5ND_{1,3,5}CPA_5d`:

   Merged parameters for the ring variant of INDCPA-secure Round5 key-encapsulation 
   mechanism,
   for NIST security levels 1, 3 and 5, resulting from the merge of the NIST PQC
   first round candidate cryptosystems _Round2_ and _HILA5_. XE5 forward error
   correction is used to decrease decryption failure rates (hence the `5` at the 
   end of the parameter designator), and improve bandwidth and security. 

* `R5ND_{1,3,5}CCA_5d`:

   Merged parameters for the ring variant of INDCCA-secure Round5 key-encapsulation mechanism and public-key encryption,
   for NIST security levels 1, 3 and 5, resulting from the merge of the NIST PQC
   first round candidate cryptosystems _Round2_ and _HILA5_. XE5 forward error
   correction is used to decrease decryption failure rates (hence the `5` at the 
   end of the parameter designator), and improve bandwidth
   and security. 
   
* `R5N1_{1,3,5}CPA_0d`:

   Parameters for the non-ring/unstructured variant of INDCPA-secure Round5 key-encapsulation 
   mechanism, for NIST security levels 1, 3 and 5. No forward error correction is used,
   hence the `0` at the end of the parameter designator.

* `R5N1_{1,3,5}CCA_0d`:

   Parameters for the non-ring/unstructured variant of INDCCA-secure Round5 key-encapsulation mechanism and public-key 
   encrpytion, 
   for NIST security levels 1, 3 and 5. No forward error correction is used,
   hence the `0` at the end of the parameter designator.

Round5 parameters for specific use-cases
----------------------------------------

* `R5ND_0CPA_2iot`:

   A _small_ INDCPA-secure key-encapsulation parameter set targeting the Internet of Things 
   use-cases, with very low bandwidth
   and computational footprints, yet still providing 88 bits of (quantum) security. XE2
   forward error correction is used to improve failure rate, bandwidth and security.

* `R5ND_1CPA_4longkey`:

   An alternative to the NIST level 3, ring variant INDCPA-secure key-encapsulation parameter 
   set `R5ND_1CPA` 
   that encapsulates a 192-bit key despite targeting the NIST security level 1, to ensure that 
   the (quantum) cost of atacking the encapsulated key (e.g., by using Grover's quantum search
   algorithm) is as much as the (quantum) cost of  attacking the underlying cryptosystem, i.e.,
   Round5.

* `R5N1_3CCA_0smallCT`:

   An alternative to the NIST level 3, non-ring variant INDCCA-secure public-key encryption 
   parameter set
   `R5N1_3CCA_0d` that has exceptionally small ciphertexts, targeting usage scenarios where the
   public-key is static and hence exchanged rarely, implying that bandwidth footprint depends
   on the size of the ciphertext. Hence this parameter set, despite enjoying the more conservative
   security assumption based on unstructured lattices, has a bandwidth requirement comparable to
   ring or structured variants.
   
2. Additional Scripts
==============

The `additional` folder includes:

* Subfolder `hybrid` contains scripts to compute estimates for the hybrid attack following Wunderer's methodology.
* Subfolder `HMAX` contains a script to compute the minimum number of calls to obtain at least `h` different values in `[0,d-1]` with a likelihood at least `1 - 2^kappa`.
* Subfolder `malformedBcheck` contains a script to compute the thresholds for the Chi2 and binomial tests. These tests are used to check whether exchange parameters B and U are malformed.
