#!/bin/bash

# The variants to analyse
VARIANTS="R5ND_0KEM_0fail_phi_0 R5ND_0KEM_0fail_phi_1 R5ND_0KEM_0fail_phi_2 R5ND_0KEM_0fail_phi_3 R5ND_0KEM_0fail_phi_4 R5ND_0KEM_0fail_phi_5 R5ND_0KEM_0fail_phi_6 R5ND_0KEM_0fail_phi_7 R5ND_0KEM_0fail_phi_8 R5ND_0KEM_0fail_phi_9 R5ND_0KEM_0fail_phi_10 R5ND_0KEM_0fail_phi_11 R5ND_0KEM_0fail_phi_12 R5ND_0KEM_0fail_phi_13 R5ND_0KEM_0fail_phi_14 R5ND_0KEM_0fail_phi_15 R5ND_0KEM_0fail_phi_16 R5ND_0KEM_0fail_phi_17 R5ND_0KEM_0fail_phi_18 R5ND_0KEM_0fail_phi_19 R5ND_0KEM_0fail_phi_20 R5ND_0KEM_0fail_phi_21 R5ND_0KEM_0fail_phi_22 R5ND_0KEM_0fail_phi_23 R5ND_0KEM_0fail_phi_24 R5ND_0KEM_0fail_phi_25 R5ND_0KEM_0fail_phi_26 R5ND_0KEM_0fail_phi_27 R5ND_0KEM_0fail_phi_28 R5ND_0KEM_0fail_phi_29 R5ND_0KEM_xfail_ntru_0 R5ND_0KEM_xfail_ntru_1 R5ND_0KEM_xfail_ntru_2 R5ND_0KEM_xfail_ntru_3 R5ND_0KEM_xfail_ntru_4 R5ND_0KEM_xfail_ntru_5 R5ND_0KEM_xfail_ntru_6 R5ND_0KEM_xfail_ntru_7 R5ND_0KEM_xfail_ntru_8 R5ND_0KEM_xfail_ntru_9 R5ND_0KEM_xfail_ntru_10 R5ND_0KEM_xfail_ntru_11 R5ND_0KEM_xfail_ntru_12 R5ND_0KEM_xfail_ntru_13 R5ND_0KEM_xfail_ntru_14 R5ND_0KEM_xfail_ntru_15 R5ND_0KEM_xfail_ntru_16 R5ND_0KEM_xfail_ntru_17 R5ND_0KEM_xfail_ntru_18 R5ND_0KEM_xfail_ntru_19 R5ND_0KEM_xfail_ntru_20 R5ND_0KEM_xfail_ntru_21 R5ND_0KEM_xfail_ntru_22 R5ND_0KEM_xfail_ntru_23 R5ND_0KEM_xfail_ntru_24 R5ND_0KEM_xfail_ntru_25 R5ND_0KEM_xfail_ntru_26 R5ND_0KEM_xfail_ntru_27 R5ND_0KEM_xfail_ntru_28 R5ND_0KEM_xfail_ntru_29"

# Hamming weights
PHI_HAMMINGWEIGHTS="[ 170, 180, 200, 220, 250, 270, 300, 320, 350, 370, 400, 420, 440, 450, 470, 500, 520, 540, 550, 570, 590, 600, 620, 640, 650, 670, 700, 720, 740, 750 ]"
NTRU_HAMMINGWEIGHTS="[ 170, 180, 200, 220, 250, 270, 300, 320, 350, 370, 400, 420, 440, 450, 470, 500, 520, 540, 550, 570, 590, 600, 620, 640, 650, 670, 700, 720, 740, 750 ]"

# The number of times the tests should be repeated
REPEAT=$1
REPEAT=${REPEAT:-1000}

# Compiler options
GCCOPTIONS="-march=native -mtune=native -O3 -fomit-frame-pointer -fwrapv"

# Function to join arguments with a comma
function join { printf "$1"; shift; printf "%s" "${@/#/, }"; }

RESULT=""
# Compile and analyse all variants
for variant in $VARIANTS; do
    # Compile code
    gcc $GCCOPTIONS -Isrc -D$variant src/*.c -lm -lcrypto -lkeccak -o failure_analysis
    # Run code
    RES=`./failure_analysis $REPEAT`
    echo "$RES"
    RESULT+="$RES
"
done

echo
echo "================================================================================"
echo
echo "repeat=$REPEAT"
echo "phiHammingWeights=$PHI_HAMMINGWEIGHTS"
echo "phiFailures=[" $(join `echo "$RESULT" | awk '/phi/ { printf $3 " " }'`) "]"
echo "ntruHammingWeights=$PHI_HAMMINGWEIGHTS"
echo "ntruFailures=[" $(join `echo "$RESULT" | awk '/ntru/ { printf $3 " " }'`) "]"

# Remove executable
rm -f failure_analysis
