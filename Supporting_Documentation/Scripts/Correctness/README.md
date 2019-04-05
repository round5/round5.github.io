# Round5 decryption failure simulation and analysis

Code and scripts for the analysis of operational correctness and decryption
failure rates in Round5, correspondng to Sec. "Correctness of Round5", and
subsection "Theoretical model and experimental results" in the specification.

## Contents

This folder contains:

* ``correctness.sage``: goal is to present and verify the model used to obtain
  failure probability estimates.

* ``runFailureAnalysis.sh`` and Round5 C code: goal is to show that the model
  used to obtain failure probability estimates matches actual experiments using
  Round5 C code.

* ``r5sim.c``: goal is to analyze decryption errors and conditional probability
  of them depending on

    - the type of reduction polynomials NTRU(x)=x^(n+1) - 1 or Phi(x) = NTRU(x)/(x-1)

    - the type of error used: pure rounding as in Round5, or idealized errors, namely uniform and gaussian.

More detailed contents are as follows:

## correctness.sage:

Prerequisite: SAGE Jupyter Notebook.

This is a SAGE script to analyze decryption failure probabilities in Round5
parameters, both those that use error correction and those that do not. More
explicitly, this script makes explicit the assumptions made in the correctness
analysis of Round5, specifically in Section "Theoretical model and experimental
results" of the Round5 specification.

The effect of dependence between polynomial coefficients in Round5 ring
parameters when polynomial multiplication is done modulo Phi(x) =
x^n+x^(n-1)+...+1 is analyzed.  Also it is analyzed how polynomial
multiplication modulo NTRU(x) = x^(n+1)-1 combined with balanced secrets
effectively solves this dependence issue and the presence of correlated errors,
thereby allowing the use of error correction.

The results of this analysis are shown by, plotting the probabilities of at
least one and two errors occurring for different parameter sets. For a more
detailed explanation of the results and how they justify the assumptions made
in the scheme's correctness analysis, please refer to the specification.

### Running correctness.sage:

To run, copy and paste this SAGE code into the SAGE Jupyter Notebook.

The primary function is `showfail()` that computes and displays the
failure rate of submitted Round5 parameters. This can be found inside the
correctness.sage script. In order to compute and display the failure rate of all submitted
parameters, please uncomment the corresponding parameter in the code.

Next, a number of "high failure" or "scaled down" Round5 parameters are
provided. The code provided computes:

* the (estimated) probability of at least *one* error occurring, and

* the (estimated) probability of at least *two* errors occurring in these
  parameters,

and plots them. This is done specifically for Round5 ring parameters that do
and do not use forward error correction, in order to justify the independence
assumption in the theoretical analysis for the latter case that allows the use
of error correction. Please see the specification for more details.


# runFailureAnalysis.sh

A shell script to test the actual failure rates of ``high-failure'' parameter
sets. Testing is done by using the actual Round5 C code and that is located in
folder src.

The only argument of the script is the number of times each algorithm variant
is tested (default: 1000 times).

The user can obtain the failure rates for high-failure parameter sets and
verify that these experimental results using Round5 C code fits with the output
coming from correctness.sage.

### Prerequisites

To be able to build and run Round5 code, the following conditions must be met:

* The OpenSSL library (preferably 1.1.1, or later) must be installed.  Use
  `apt-get install libssl-dev` for most Linux distributions.  On a Mac, an easy
  way is to use [brew](https://brew.sh), install it with `brew install
  openssl@1.1` and then add it to the `CPATH` and `LIBRARY_PATH` environment
  variables:

  ```
  export CPATH=${CPATH+$CPATH:}/usr/local/opt/openssl@1.1/include
  export LIBRARY_PATH=${LIBRARY_PATH+$LIBRARY_PATH:}/usr/local/opt/openssl@1.1/openssl/lib
  ```

* The Keccak library must be installed in your system.  This is done as
  follows:

  1. Linux: Install xsltproc (e.g.  `sudo apt-get install xsltproc`, on a Mac
     it should already be installed).

  2. Clone the [XKCP git repository](https://github.com/XKCP/XKCP.git).

  3. Build the library, e.g. using `make generic64/libkeccak.a` (for
     alternatives, see the github page).

  4. Add the library to your `CPATH` and `LIBRARY_PATH` environment variables:

    ```
    export CPATH=${CPATH+$CPATH:}<XKCP_DIRECTORY>/bin/generic64
    export LIBRARY_PATH=${LIBRARY_PATH+$LIBRARY_PATH:}<XKCP_DIRECTORY>/bin/generic64
    ```

## r5sim.c

This file contains C code to run a (simplified) version of the ring variant of
the Round5 cpa public-key encryption protocol in order to simulate decryption
failures in Round5, and statistically derive the failure rate.

While with "runFailureAnalysis.sh", we analyze Round5 and it exact
specification, r5sim.c code is more generic in the sense that it allows
incorporating different types of noise, namely rounding, uniform and
gaussian. This is done to identify the potential influence of the type of used
noise/reduction polynomial in decryption failures or error correlations.

The user can run the script and obtain decryption failure rates for the
different types of noise and the different types of reduction polynomials,
namely NTRU(X) and PHI(X) to validate the rationale in Round5's submission,
namely:

* rounding errors do not introduce major correlations compared with a scheme
  based on gaussian errors,

* Round5 parameter sets with error correction do not suffer of severe failure
  correlations.

### Round5 variants simulated in r5sim.c

r5sim.c is a very simplified version of Round5 that computes the differences
between the ``raw keys'' computed by the initiator/Alice/decryptor and the
responder/Bob/encryptor.

In r5sim.c some Round5 functions are not implemented, e.g., Sample \mu function
is ignored in the simulation.  Errors are computed over the entire `n`-long
polynomial. All rounding is to the closest integer.

There are three variants of Round5 simulated, for three different kinds of
noise, in order to demonstrate different error resilience and correlation
behaviors for different kinds of noise.

* Round5 with actual rounding noise.

* Round5, with continuous Gaussian noise, to demonstrate how correlation
  between rounding noise and the secret-keys affects the final failure rate, in
  the most ideal version of the protocol.

* Round5, with uniform additive noise, to demonstrate how correlation between
  rounding noise and the secret-keys affects the final failure rate.

### Output of r5sim.c

Output has three error probabilities in log2:

* the first is the probability of 1 bit error occurring.

* the second is the conditional probability of 1 bit error occurring assuming
  that *one* other bit error has occurred (this shows the effect of correlated
  errors in the prime cyclotomic ring).

* the third is the conditional probability of 1 bit error occurring assuming
  that *two* other bit errors have occurred; this further shows how the issue
  of correlated errors is effectively solved in Round5 ring parameters that use
  error correction, due to polynomial multiplication modulo x^(n+1)-1 combined
  with balanced secrets.

### Running r5sim.c

There are no prerequisites. Open the Makefile to view or edit the Round5
parameters that will be simulated -- these include:

* the ring dimension `n`,

* the Hamming weight of Round5 secret-key vectors `h`,

* the exponent of the Round5 large modulus `q` (a power of 2),

* the exponent of the Round5 rounding modulus `p` (a power of 2),

* the exponent of the Round5 ciphertext compression modulus `t` (a power of 2),

* the number of simulations `runs` that will be run, this must be large in
  order to derive any meaningful statistics, on the order of 10^6 or more.

* the number of `intermediate_runs` after which data is periodically written
  into/saved to the simulation output `.out` file. Defaults to the same value
  as `runs`.

Run `make` to build the simulation source and run the simulation executable.
The simulation results are written to an `.out` file, that can be viewed with
any text editor.

The following make targets are available:

* `build`:

   Builds the simulation source into an executable.

* `r5_rounding`:

   Runs a Round5 simulation executable where the noise is generated by rounding
   down from a larger to smaller modulus, passing the above Round5 parameters
   to it as arguments.

* `r5_gaussiannoise`:

   Runs a Round5 simulation executable where the noise is sampled from a
   continuous Gaussian distribution, passing the above Round5 parameters to it
   as arguments. The variance of the distribution is derived from the Round5
   moduli `q`, `p` and `t`, so as to be balanced with that of equivalent
   rounding noise for said moduli.

* `r5_uniformnoise`:

   Runs a Round5 simulation executable where the noise is generated from a
   uniform, bounded distribution with parameters derived from the Round5
   moduli.

* `all`:

   Builds the simulation executable and runs it. This is the default target.

* `clean`:

   Remove all build artifacts.

* `cleanfull`:

   Remove all build artifacts, and additionally all existing simulation
   results.
