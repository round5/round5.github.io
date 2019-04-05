// Failure analysis test application

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include "r5_parameter_sets.h"
#include "r5_cpa_pke.h"
#include "misc.h"
#include "rng.h"

#if PARAMS_TAU == 1
#include "a_fixed.h"
#endif

/* Counts number of bit errors */
static unsigned count_bit_errors(unsigned char m_r[PARAMS_KAPPA_BYTES], unsigned char m_i[PARAMS_KAPPA_BYTES]) {
    unsigned errors = 0;
    for (int i = 0; i < PARAMS_KAPPA_BYTES; ++i){
        for (int j = 0; j < 8 ; ++j){
            if (((m_r[i] >> j) & 0x1) != ((m_i[i] >> j) & 0x1)) {
                ++errors;
            }
        }
    }
    return errors;
}

int main(int argc, char **argv) {
    long repeat = 0;
    if (argc > 1) {
         repeat = strtol(argv[1], NULL, 10);
    }
    if (repeat <= 0) {
         repeat = 1000;
    }

    unsigned long nr_failed = 0;
    unsigned long bit_error_counts[8] = {0};

    /* Initialize random bytes RNG */
    unsigned char entropy_input[48];
    int i;
    for (i = 0; i < 48; i++) {
        entropy_input[i] = (unsigned char) i;
    }
    randombytes_init(entropy_input, NULL, 256);

#if PARAMS_TAU == 1
    uint8_t seed[PARAMS_KAPPA_BYTES];
    randombytes(seed, PARAMS_KAPPA_BYTES);
    create_A_fixed(seed);
#endif

    unsigned char pk[PARAMS_PK_SIZE];
    unsigned char sk[PARAMS_KAPPA_BYTES];
    unsigned char ct[PARAMS_CT_SIZE];
    uint8_t m_r[PARAMS_KAPPA_BYTES];
    uint8_t m_i[PARAMS_KAPPA_BYTES];
    uint8_t rho[PARAMS_KAPPA_BYTES];

    for (long r = 0; r < repeat; ++r) {
        randombytes(m_r, PARAMS_KAPPA_BYTES);
        randombytes(rho, PARAMS_KAPPA_BYTES);
        r5_cpa_pke_keygen(pk, sk);
        r5_cpa_pke_encrypt(ct, pk, m_r, rho);
        r5_cpa_pke_decrypt(m_i, sk, ct);
        if (memcmp(m_r, m_i, PARAMS_KAPPA_BYTES)) {
            ++nr_failed;
            unsigned count = count_bit_errors(m_r, m_i);
            if (count < 8){
                bit_error_counts[count]++;
            } else {
                bit_error_counts[7]++;
            }
        } else {
            bit_error_counts[0]++;
        }
    }
    printf("%s failed %lu out of %ld times (%0.2f%%)\n", CRYPTO_ALGNAME, nr_failed, repeat, 100.0 * nr_failed / repeat);
    printf("Bit error counts: [%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu]\n",bit_error_counts[0],bit_error_counts[1], bit_error_counts[2],bit_error_counts[3],bit_error_counts[4],bit_error_counts[5],bit_error_counts[6],bit_error_counts[7]);
}
