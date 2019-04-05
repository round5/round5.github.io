Round5 Parameter Summary
========================

This folder contains the script `summarize.py` that can be used to summarize
**all** parameters proposed as part of the Round5 cryptosystem. 
Open terminal and run:

``  python summarize.py `` 

In total, **18 main parameters** are summarized, along with 3 additional “specific use-case”
parameters that demonstrate the flexibility of Round5 and applicability to 3
different  specialized usage scenarios.

The parameters themselves are summarized, along with security levels and 
performance details, including bandwidth requirements and decryption failure
rates.

Note that the script also includes the "High-Failure" parameters that are used to verify
failure estimates and error correlations. 

Main Round5 parameters
----------------------

* `R5ND_{1,3,5}KEM_0d`:

    Parameters for the ring variant of INDCPA-secure Round5 key-encapsulation 
    mechanism,
    for NIST security levels 1, 3 and 5. No forward error correction is used,
    hence the `0` at the end of the parameter designator. All polynomial
    multiplications are done modulo the prime-order cyclotomic polynomial. 

* `R5ND_{1,3,5}PKE_0d`:

   Parameters for the ring variant of INDCCA-secure Round5 public-key encryption,
   for NIST security levels 1, 3 and 5. No forward error correction is used,
   hence the `0` at the end of the parameter designator. All polynomial
   multiplications are done modulo the prime-order cyclotomic polynomial. 

* `R5ND_{1,3,5}KEM_5d`:

   Merged parameters for the ring variant of INDCPA-secure Round5 key-encapsulation 
   mechanism,
   for NIST security levels 1, 3 and 5, resulting from the merge of the NIST PQC
   first round candidate cryptosystems _Round2_ and _HILA5_. XE5 forward error
   correction is used to decrease decryption failure rates (hence the `5` at the 
   end of the parameter designator), and improve bandwidth and security. 

* `R5ND_{1,3,5}PKE_5d`:

   Merged parameters for the ring variant of INDCCA-secure Round5 public-key encryption,
   for NIST security levels 1, 3 and 5, resulting from the merge of the NIST PQC
   first round candidate cryptosystems _Round2_ and _HILA5_. XE5 forward error
   correction is used to decrease decryption failure rates (hence the `5` at the 
   end of the parameter designator), and improve bandwidth
   and security. 
   
* `R5N1_{1,3,5}KEM_0d`:

   Parameters for the non-ring/unstructured variant of INDCPA-secure Round5 key-encapsulation 
   mechanism, for NIST security levels 1, 3 and 5. No forward error correction is used,
   hence the `0` at the end of the parameter designator.

* `R5N1_{1,3,5}PKE_0d`:

   Parameters for the non-ring/unstructured variant of INDCCA-secure Round5 public-key 
   encrpytion, 
   for NIST security levels 1, 3 and 5. No forward error correction is used,
   hence the `0` at the end of the parameter designator.

Round5 parameters for specific use-cases
----------------------------------------

* `R5ND_0KEM_2iot`:

   A _small_ INDCPA-secure key-encapsulation parameter set targeting the Internet of Things 
   use-cases, with very low bandwidth
   and computational footprints, yet still providing 88 bits of (quantum) security. XE2
   forward error correction is used to improve failure rate, bandwidth and security.

* `R5ND_1KEM_4longkey`:

   An alternative to the NIST level 3, ring variant INDCPA-secure key-encapsulation parameter 
   set `R5ND_1KEM` 
   that encapsulates a 192-bit key despite targeting the NIST security level 1, to ensure that 
   the (quantum) cost of atacking the encapsulated key (e.g., by using Grover's quantum search
   algorithm) is as much as the (quantum) cost of  attacking the underlying cryptosystem, i.e.,
   Round5.

* `R5N1_3PKE_0smallCT`:

   An alternative to the NIST level 3, non-ring variant INDCCA-secure public-key encryption 
   parameter set
   `R5N1_3PKE_0d` that has exceptionally small ciphertexts, targeting usage scenarios where the
   public-key is static and hence exchanged rarely, implying that bandwidth footprint depends
   on the size of the ciphertext. Hence this parameter set, despite enjoying the more conservative
   security assumption based on unstructured lattices, has a bandwidth requirement comparable to
   ring or structured variants.
