#ifndef EMISSION_CACHE_H
#define EMISSION_CACHE_H

#include <stddef.h>
#include "Logodd.h"
#include "Hmm.h"

/*
 * Precompute emission log-odds for all states × all timesteps.
 * Layout: TIME-MAJOR — cache[t * S + s]
 * This means all states at the same t are contiguous in memory,
 * which is optimal for the Viterbi inner loop access pattern.
 *
 * obs1[t] = pre-encoded 1-mer observation at time t
 * obs2[t] = pre-encoded 2-mer observation at time t  
 * obs3[t] = pre-encoded 3-mer observation at time t
 * (these arrays already exist in Viterbi.c — pass them in directly)
 *
 * Caller owns the returned pointer and must free() it.
 */
LOGODD_T* EmissionCache__build(struct HMM* hmm,
                                size_t T,
                                uint16_t* obs1,
                                uint16_t* obs2,
                                uint16_t* obs3);

#endif /* EMISSION_CACHE_H */
