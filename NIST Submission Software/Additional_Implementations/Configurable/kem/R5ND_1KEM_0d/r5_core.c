/*
 * Copyright (c) 2018, Koninklijke Philips N.V.
 */

/**
 * @file
 * Implementation of the core algorithm functions.
 */

#include "r5_core.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "misc.h"
#include "r5_memory.h"
#include "rng.h"
#include "little_endian.h"
#include "drbg.h"
#include "r5_hash.h"
#include "a_fixed.h"
#include "a_random.h"

/*******************************************************************************
 * Private functions
 ******************************************************************************/

/**
 * Create sparse ternary vector(s).
 *
 * @param[out] vectors_idx the generated vector(s) in index form:
 *                         the first h/2 items are indexes with positive values,
 *                         the second h/2 items are indexes with negative values
 * @param[in]  seed        the seed to use for generating the random data
 * @param[in]  nr_vectors  the number of vectors to create
 * @param[in]  params      the algorithm parameters
 * @return __0__ in case of success
 */
static int create_secret_vectors_idx(uint16_t *vectors_idx, const unsigned char *seed, const unsigned nr_vectors, const parameters *params) {
    unsigned i, j;
    uint16_t idx;
    unsigned char *occupied = checked_malloc(params->d);

    /* Initialize drbg */
    drbg_init(seed, params->kappa_bytes);
    DRBG_SAMPLER16_INIT(params->d);

    size_t pos_idx = 0;
    size_t neg_idx = params->h / 2;
    for (j = 0; j < nr_vectors; ++j) {
        memset(occupied, 0, params->d);
        for (i = 0; i < params->h; ++i) {
            do {
                DRBG_SAMPLER16(idx, params->d);
            } while (occupied[idx] != 0);
            occupied[idx] = 1;
            if (i & 1) {
                vectors_idx[neg_idx++] = idx;
            } else {
                vectors_idx[pos_idx++] = idx;
            }
        }
        pos_idx += params->h / 2;
        neg_idx += params->h / 2;
    }

    free(occupied);

    return 0;
}

/**
 * Multiplies a polynomial in the cyclotomic ring times (X - 1), the result can
 * be taken to be in the NTRU ring X^(len+1) - 1.
 *
 * @param[out] ntru_pol  result
 * @param[in]  cyc_pol   polynomial in the cyclotomic ring
 * @param[in]  len       number of coefficients of the cyclotomic polynomial
 * @return __0__ in case of success
 */
static int lift_poly(uint16_t *ntru_pol, const int16_t *cyc_pol, const size_t len) {
    size_t i;

    ntru_pol[0] = (uint16_t) (-cyc_pol[0]);

    for (i = 1; i < len; ++i) {
        ntru_pol[i] = (uint16_t) (cyc_pol[i - 1] - cyc_pol[i]);
    }

    ntru_pol[len] = (uint16_t) cyc_pol[len - 1];

    return 0;
}

int unlift_poly(uint16_t *cyc_pol, const uint16_t *ntru_pol, size_t len) {
    size_t i;

    cyc_pol[0] = (uint16_t) (-ntru_pol[0]);
    for (i = 1; i < len; ++i) {
        cyc_pol[i] = (uint16_t) (cyc_pol[i - 1] - ntru_pol[i]);
    }

    return 0;
}

/**
 * Computes the coefficient with power _l_ of _X = B^T * R_.
 *
 * @param[in] B __B__ (extended and re-arranged)
 * @param[in] R_idx __R__ in index form
 * @param[in] l the power of the coefficient (index)
 * @param[in] params the algorithm parameters in use
 * @return __X[l]__
 */
static inline uint16_t compute_BTR_coefficient_ring(const uint16_t *B, const uint16_t *R_idx, const uint32_t l, const parameters *params) {
    uint32_t k;

    uint16_t B_val;
    uint16_t X_val = 0;

    for (k = 0; k < params->h / 2; ++k) {
        B_val = B[R_idx[k] + l];
        X_val = (uint16_t) (X_val + B_val);
    }

    for (k = params->h / 2; k < params->h; ++k) {
        B_val = B[R_idx[k] + l];
        X_val = (uint16_t) (X_val - B_val);
    }

    return X_val;
}

/**
 * Computes _X = B^T * R_ in case of non-ring parameters.
 *
 * @param[out] X __X__
 * @param[in] B __B__
 * @param[in] R_idx __R__ in index form
 * @param[in] params the algorithm parameters in use
 * @return
 */
static int compute_BTR_non_ring(uint16_t *X, uint16_t *B, const uint16_t *R_idx, const parameters *params) {
    size_t idx;
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t k;
    size_t index;
    uint16_t B_val;
    uint16_t X_val;
    for (idx = 0; idx < params->mu; ++idx) {
        if (j == params->m_bar) {
            j = 0;
            ++i;
        }
        X_val = 0;
        for (k = 0; k < params->h / 2; ++k) {
            index = (size_t) ((R_idx[j * params->h + k]));
            B_val = B[index * params->n_bar + i];
            X_val = (uint16_t) (X_val + B_val);
        }
        for (k = params->h / 2; k < params->h; ++k) {
            index = (size_t) ((R_idx[j * params->h + k]));
            B_val = B[index * params->n_bar + i];
            X_val = (uint16_t) (X_val - B_val);
        }
        X[idx] = X_val;
        ++j;
    }
    return 0;
}

/**
 * Computes _X' = U * S_ in case of non-ring parameters.
 *
 * @param[out] X_prime __X'__
 * @param[in] U_T  __U^T__
 * @param[in] S_idx __S__ in index form
 * @param[in] params the algorithm parameters in use
 * @return
 */
static int compute_US_non_ring(uint16_t *X_prime, uint16_t *U_T, const uint16_t *S_idx, const parameters *params) {
    size_t idx;
    uint32_t i = 0;
    uint32_t j = 0;
    uint32_t k;
    size_t index;
    uint16_t U_T_val;
    uint16_t X_prime_val;
    for (idx = 0; idx < params->mu; ++idx) {
        if (i == params->m_bar) {
            i = 0;
            ++j;
        }
        X_prime_val = 0;
        for (k = 0; k < params->h / 2; ++k) {
            index = (size_t) ((S_idx[j * params->h + k]));
            U_T_val = U_T[index + i * params->d];
            X_prime_val = (uint16_t) (X_prime_val + U_T_val);
        }
        for (k = params->h / 2; k < params->h; ++k) {
            index = (size_t) ((S_idx[j * params->h + k]));
            U_T_val = U_T[index + i * params->d];
            X_prime_val = (uint16_t) (X_prime_val - U_T_val);
        }
        X_prime[idx] = X_prime_val;
        ++i;
    }
    return 0;
}

/**
 * Computes _X = B^T * R_ in case of ring parameters.
 *
 * @param[out] X __X__
 * @param[in] B __B^T__
 * @param[in] R_idx __R__ in index form
 * @param[in] params the algorithm parameters in use
 * @return
 */
static int compute_BTR_ring(uint16_t *X, uint16_t *B, const uint16_t *R_idx, const parameters *params) {
    uint32_t i = 0;
    uint16_t *B_aux;
    uint16_t *X_aux;
    const uint32_t size_B_aux = (uint32_t) (params->d + 1);

    /* First we extend (and lift) B */
    B_aux = checked_malloc((size_t) (2 * size_B_aux) * sizeof (*B_aux));

    if (params->xe == 0 && params->f == 0) {
        /* Move to NTRU ring */
        lift_poly(B_aux, (int16_t*) B, params->d);
    } else {
        memcpy(B_aux, B, params->d * sizeof (*B));
        B_aux[params->d] = 0;
    }

    /* Rearrange elements */
    uint16_t *tmp = checked_malloc(size_B_aux * sizeof (*tmp));
    memcpy(tmp, B_aux, size_B_aux * sizeof (*B_aux));
    for (i = 1; i < size_B_aux; ++i) {
        B_aux[i] = tmp[size_B_aux - i];
    }
    free(tmp);

    /* Duplicate vector to remove need of modulo operation */
    memcpy(B_aux + size_B_aux, B_aux, size_B_aux * sizeof (*B_aux));

    /* Temp variable to store the results */
    X_aux = checked_calloc((size_t) (params->mu + 1), sizeof (*X_aux));

    /* Compute X */
    uint32_t idx;
    X_aux[0] = compute_BTR_coefficient_ring(B_aux, R_idx, 0, params);
    for (idx = 1; idx < params->mu + 1U; ++idx) {
        X_aux[idx] = compute_BTR_coefficient_ring(B_aux, R_idx, params->d + 1U - idx, params);
    }

    if (params->xe == 0 && params->f == 0) {
        /* In case of the ring, convert back to cyclotomic polynomial*/
        unlift_poly(X, X_aux, params->mu);
    } else {
        memcpy(X, X_aux + 1, params->mu * sizeof (*X));
    }

    free(B_aux);
    free(X_aux);

    return 0;
}

/**
 * Generates the row displacements for the A matrix creation variant tau=0.
 * This permutation creates a convolution matrix assuming the first row is the
 * polynomial ordered as a_0, a_(n-1), a_(n-2), ...
 *
 * @param[out] row_disp the row displacements
 * @param[in]  params   the algorithm parameters in use
 * @return __0__ on success
 */
static int permutation_tau_0_ring(uint32_t *row_disp, const parameters *params) {
    uint32_t i;
    uint32_t n_rows = (uint32_t) (params->d + 1);

    for (i = 0; i < n_rows; ++i) {
        row_disp[i] = (uint32_t) (params->d + 1) - i;
    }

    return 0;
}

/**
 * Generates the row displacements for the A matrix creation variant tau=0.
 * Note: This is the identity mapping!
 *
 * @param[out] row_disp the row displacements
 * @param[in]  params   the algorithm parameters in use
 * @return __0__ on success
 */
static int permutation_tau_0_non_ring(uint32_t *row_disp, const parameters *params) {
    uint32_t i;

    for (i = 0; i < params->d; ++i) {
        row_disp[i] = i * params->d;
    }

    return 0;
}

/**
 * The DRBG customization when creating the tau=1 or tau=2 permutations.
 */
static const uint8_t permutation_customization[2] = {0, 1};

/**
 * Generates the row displacements for the A matrix creation variant tau=1.
 *
 * @param[out] row_disp  the row displacements
 * @param[in]  seed      the seed
 * @param[in]  params    the algorithm parameters in use
 * @return __0__ on success
 */
static int permutation_tau_1(uint32_t *row_disp, const unsigned char *seed, const parameters *params) {
    uint32_t i;
    uint16_t rnd;
    DRBG_SAMPLER16_INIT(params->d);

    drbg_init_customization(seed, params->kappa_bytes, permutation_customization, sizeof (permutation_customization));

    for (i = 0; i < params->d; ++i) {
        DRBG_SAMPLER16(rnd, params->d);
        row_disp[i] = 2 * i * params->d + rnd;
    }

    return 0;
}

/**
 * Generates the row displacements for the A matrix creation variant tau=2.
 *
 * @param[out] row_disp  the row displacements
 * @param[in]  seed      the seed
 * @param[in]  params    the algorithm parameters in use
 * @return __0__ on success
 */
static int permutation_tau_2(uint32_t *row_disp, const unsigned char *seed, const parameters *params) {
    uint32_t i;
    uint16_t rnd;
    uint8_t *v = checked_calloc(params->tau2_len, 1);

    drbg_init_customization(seed, params->kappa_bytes, permutation_customization, sizeof (permutation_customization));

    for (i = 0; i < params->k; ++i) {
        do {
            rnd = (uint16_t) (drbg_sampler16_2(params->tau2_len) & (params->tau2_len - 1));
        } while (v[rnd]);
        v[rnd] = 1;
        row_disp[i] = rnd;
    }

    free(v);

    return 0;
}

/**
 * Generates A_master, allocates the necessary space
 *
 * @param[out]  A_master  pointer to a variable holding A_master
 * @param[in]   sigma     seed
 * @param[in]   params    the algorithm parameters in use
 * @return __0__ on success
 */
static int create_A_master(uint16_t **A_master, const unsigned char *sigma, const parameters *params) {
    size_t i;

    switch (params->tau) {
        case 0:
            if (params->k == 1) {
                *A_master = checked_malloc((2 * (size_t) (params->d + 1) * sizeof (**A_master)));
                create_A_random(*A_master, sigma, params);
                uint16_t *aux = checked_malloc((size_t) (params->d + 1) * sizeof (*aux));
                lift_poly(aux, (int16_t *) * A_master, params->d);
                (*A_master)[0] = aux[0];
                for (i = 1; i < (size_t) (params->d + 1); ++i) {
                    (*A_master)[i] = aux[(size_t) (params->d + 1) - i];
                }
                memcpy((*A_master) + (params->d + 1), *A_master, (size_t) (params->d + 1) * sizeof (**A_master));
                free(aux);
            } else {
                *A_master = checked_malloc((size_t) (params->k * params->d) * sizeof (**A_master));
                create_A_random(*A_master, sigma, params);
            }
            break;
        case 1:
            assert(A_fixed != NULL && A_fixed_len == (size_t) (2 * params->d * params->k));
            *A_master = A_fixed;
            break;
        case 2:
            *A_master = checked_malloc((size_t) (params->tau2_len + params->d) * sizeof (**A_master));
            create_A_random(*A_master, sigma, params);
            memcpy((*A_master) + params->tau2_len, *A_master, params->d * sizeof (**A_master));

            break;
    }

    return 0;
}

/*******************************************************************************
 * Public functions
 ******************************************************************************/

int create_A(uint16_t **A_master, uint32_t *A_permutation, const unsigned char *sigma, const parameters *params) {
    /* Create A_master */
    create_A_master(A_master, sigma, params);

    /* Compute the permutation */
    assert(params->tau <= 2);
    switch (params->tau) {
        case 0:
            if (params->k == 1) {
                permutation_tau_0_ring(A_permutation, params);
            } else {
                permutation_tau_0_non_ring(A_permutation, params);
            }
            break;
        case 1:
            permutation_tau_1(A_permutation, sigma, params);
            break;
        case 2:
            permutation_tau_2(A_permutation, sigma, params);
            break;
    }

    return 0;
}

int create_S(uint16_t *S_idx, const unsigned char *sk, const parameters *params) {
    return create_secret_vectors_idx(S_idx, sk, params->n_bar, params);
}

int create_R(uint16_t *R_idx, const unsigned char *rho, const parameters *params) {
    return create_secret_vectors_idx(R_idx, rho, params->m_bar, params);
}

int compute_AS(uint16_t *B, const uint16_t *A_master, const uint32_t *A_permutation, const uint16_t *S_idx, const parameters *params) {
    uint32_t i, j, k;
    uint32_t size_i, size_j;
    size_t idx = 0;
    uint16_t A_val, B_val;
    uint32_t A_idx;
    uint16_t *B_aux;

    if (params->k == 1) {
        /* Ring */
        assert(params->n_bar == 1 && params->m_bar == 1);
        size_i = params->n + 1;
        size_j = 1;
        B_aux = checked_malloc(size_i * sizeof (*B_aux));
    } else {
        /* Non-ring */
        size_i = params->d;
        size_j = params->n_bar;
        B_aux = B;
    }

    for (i = 0; i < size_i; ++i) {
        for (j = 0; j < size_j; ++j) {
            B_val = 0;
            for (k = 0; k < params->h / 2; ++k) {
                /* Positions where S = 1 */
                A_idx = (uint32_t) (S_idx[j * params->h + k] + A_permutation[i]);
                A_val = A_master[A_idx];
                B_val = (uint16_t) (B_val + A_val);
            }
            for (k = params->h / 2; k < params->h; ++k) {
                /* Positions where S = -1 */
                A_idx = (uint32_t) (S_idx[j * params->h + k] + A_permutation[i]);
                A_val = A_master[A_idx];
                B_val = (uint16_t) (B_val - A_val);
            }
            B_aux[idx] = B_val;
            ++idx;
        }
    }

    if (params->k == 1) {
        /* Unlift */
        unlift_poly(B, B_aux, params->n);
        free(B_aux);
    }

    return 0;
}

int compute_RTA(uint16_t *U_T, const uint16_t *A_master, const uint32_t *A_permutation, const uint16_t *R_idx, const parameters *params) {
    if (params->k != 1) {
        /* Non-ring */
        uint32_t i, j, k;
        const size_t len_u = (size_t) (params->d * params->m_bar);

        uint16_t A_val;
        size_t A_row, A_row_idx;

        memset(U_T, 0, len_u * sizeof (*U_T));

        for (j = 0; j < params->m_bar; ++j) {
            for (k = 0; k < params->h / 2; ++k) {
                A_row = (size_t) ((R_idx[j * params->h + k]));
                A_row_idx = A_permutation[A_row];
                for (i = 0; i < params->d; ++i) {
                    A_val = A_master[i + A_row_idx];
                    U_T[j * params->d + i] = (uint16_t) (U_T[j * params->d + i] + A_val);
                }
            }
            for (k = params->h / 2; k < params->h; ++k) {
                A_row = (size_t) ((R_idx[j * params->h + k]));
                A_row_idx = A_permutation[A_row];
                for (i = 0; i < params->d; ++i) {
                    A_val = A_master[i + A_row_idx];
                    U_T[j * params->d + i] = (uint16_t) (U_T[j * params->d + i] - A_val);
                }
            }
        }

        return 0;
    } else {
        /* Ring */
        assert(params->n_bar == 1 && params->m_bar == 1);
        /* With ring, A^T == A and U^T == U so (A^T*R)^T = A*R and since S and R in this case are
         * the same dimensions we can use the same function as when computing B */
        return compute_AS(U_T, A_master, A_permutation, R_idx, params);
    }
}

int compute_BTR(uint16_t *X, uint16_t *B, const uint16_t *R_idx, const parameters *params) {
    if (params->k != 1) {
        /* Non-ring */
        return compute_BTR_non_ring(X, B, R_idx, params);
    } else {
        /* Ring */
        assert(params->n_bar == 1 && params->m_bar == 1);
        return compute_BTR_ring(X, B, R_idx, params);
    }
}

int compute_STU(uint16_t *X_prime, uint16_t *U_T, const uint16_t *S_idx, const parameters *params) {
    if (params->k != 1) {
        /* Non-ring */
        return compute_US_non_ring(X_prime, U_T, S_idx, params);
    } else {
        /* Ring */
        assert(params->n_bar == 1 && params->m_bar == 1);
        /* With ring, X' == X'^T and U^T == U so X'^T = (S^T*U)^T = U^T*S = U*S
         * and since U and B and S and R in this case are the same dimensions we
         * can use the same function as when computing X = B^T*R */
        return compute_BTR_ring(X_prime, U_T, S_idx, params);
    }
}

int decompress_matrix(uint16_t *matrix, const size_t len, const size_t els, const uint16_t a_bits, const uint16_t b_bits) {
    size_t i;

    const uint16_t shift = (uint16_t) (b_bits - a_bits);
    for (i = 0; i < len * els; ++i) {
        matrix[i] = (uint16_t) (matrix[i] << shift);
    }

    return 0;
}

int round_matrix(uint16_t *matrix, const size_t len, const size_t els, const uint16_t a_bits, const uint16_t b_bits, const uint16_t rounding_constant) {
    size_t i;
    const uint16_t shift = (uint16_t) (a_bits - b_bits);
    for (i = 0; i < len * els; ++i) {
        matrix[i] = (uint16_t) ((matrix[i] + rounding_constant) >> shift);
    }

    return 0;
}

