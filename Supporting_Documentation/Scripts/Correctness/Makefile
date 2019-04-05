# Makefile for simulating a simplified version of the
# Round5 CPA public-key encryption building block,
# more specifically, encryption and decryption
# by the two parties and differences in the ``shared 
# secrets'' computed by each party as part of this. The goal
# is to collect statistical information on decryption
# failures in Round5, specifically:
#
# a) bit error probabilities,
#
# b) conditional error probabilities, assuming errors
#    have already occurred, in order to demonstrate 
#    the effect of correlated errors on the overall 
#    error rate and the usage of forward error correction.
#
# The above is computed for two instances: 
#
# one in which
# (a) the ciphertext component is computed via polynomial
# multiplication modulo  N_{n+1}(x) = x^{n+1}(x),
# and 
# 
# (b) one in which it is computed modulo \Phi_{n+1}(x) =
# [ N_{n+1}(x) ] / (x-1).
#
# This is done to indicate that correlated errors in the
# N_{n+1}(x) is unlikely and justify the use of error
# correction in this case, by showing the different behavior of 
# conditional error probabilities in the two instances.
# See specification for more details.
#
# Apart from the two Round5 variants above with differing
# reduction polynomials, different noise distributions can also be 
# simulated, in order to compare their effect on the decryption failure
# and correlated errors. 
#
# (a) Rounding noise,
# (b) Independent, additive uniform noise,
# (c) Additive, continuous Gaussian noise.


# Round5 parameters to simulate
n                  = 800
h                  = 300
qexp               = 12
pexp               = 7
texp               = 5


# number of simulation iterations
runs               = 10000
# controls after how many intermediate
# runs output is flushed to file
intermediate_runs := $(runs)


# Files
outfile           := $(CURDIR)/R5sim_$(n)-$(h)-$(qexp)-$(pexp)-$(texp)
outfilepref       := $(CURDIR)/R5sim_
now               := $(shell /bin/date "+%m%d%H%M%Y%S")   # timestamp
r5outfile         := $(addsuffix .out, $(addsuffix _$(now),$(outfile)) )


# Compiler setup and flags
CC                 = gcc
CFLAGS            += -std=c99 -Wall -Wextra
LIBS               = -lm
SRC                = r5sim.c
OBJ                = r5sim.o
# Timestamping the executable name so we can run
# multiple parallel instances of this simulation
#
r5simobj           = r5simobj
BIN               := $(CURDIR)/$(r5simobj)$(now)


# Make targets
all: build r5_rounding r5_gaussiannoise r5_uniformnoise

build:
	$(CC) $(CFLAGS) -o $(BIN) $(SRC) $(LIBS)

# Run simulations
# 1st parameter: Whether XEf forward error correction is used or not.
# This determines also whether the ciphertext component v is computed 
# via polynomial multiplications modulo x^(n+1)-1 or modulo 
# 1+x+...+x^n.

r5_rounding: build
	@# Actual Round5
	@#                                                              Gaussian noise Uniform noise Rounding noise
	@$(BIN) 0 $(n) $(h) $(qexp) $(pexp) $(texp) $(runs) $(r5outfile)             0             0              1 $(intermediate_runs)
	@$(BIN) 1 $(n) $(h) $(qexp) $(pexp) $(texp) $(runs) $(r5outfile)             0             0              1 $(intermediate_runs)

r5_gaussiannoise: build
	@# Additive Gaussian noise Round5
	@#                                                              Gaussian noise Uniform noise Rounding noise
	@$(BIN) 0 $(n) $(h) $(qexp) $(pexp) $(texp) $(runs) $(r5outfile)             1             0              0 $(intermediate_runs)
	@$(BIN) 1 $(n) $(h) $(qexp) $(pexp) $(texp) $(runs) $(r5outfile)             1             0              0 $(intermediate_runs)

r5_uniformnoise: build
	@# Additive uniform noise Round5
	@#                                                              Gaussian noise Uniform noise Rounding noise
	@$(BIN) 0 $(n) $(h) $(qexp) $(pexp) $(texp) $(runs) $(r5outfile)             0             1              0 $(intermediate_runs)
	@$(BIN) 1 $(n) $(h) $(qexp) $(pexp) $(texp) $(runs) $(r5outfile)             0             1              0 $(intermediate_runs)

clean:
	@rm -f $(OBJ)
	@rm -f $(wildcard $(r5simobj)*)

cleanall: clean
	$(info Warning! Removing all simulation result files [${wildcard $(outfilepref)*.out}])
	@rm -f $(wildcard $(outfilepref)*.out)

.PHONY: all r5_rounding r5_gaussiannoise r5_uniformnoise clean cleanall
