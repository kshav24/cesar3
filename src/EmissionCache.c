#include <stdlib.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "EmissionCache.h"
#include "EmissionTable.h"
#include "Logodd.h"

LOGODD_T* EmissionCache__build(struct HMM* hmm,
                                size_t T,
                                uint16_t* obs1,
                                uint16_t* obs2,
                                uint16_t* obs3) {
    size_t S = hmm->num_states;
    LOGODD_T* cache = (LOGODD_T*) malloc(S * (T + 1) * sizeof(LOGODD_T));

    /*
     * TIME-MAJOR layout: cache[t * S + s]
     * All states at same t are adjacent — optimal for Viterbi row scan.
     * Each (t, s) cell is independent → perfectly parallel.
     */
    #pragma omp parallel for schedule(dynamic, 8) collapse(2)
    for (size_t t = 0; t <= T; t++) {
        for (size_t s = 0; s < S; s++) {
            struct State* st = &hmm->states[s];
            LOGODD_T e;

            if (st->num_emissions == 0) {
                e = 0;                  /* silent state: no emission */
            } else if (t < st->num_emissions) {
                e = LOGODD_NEGINF;      /* too early to emit */
            } else {
                uint16_t col = 0;
                if      (st->num_emissions == 1) col = obs1[t - 1];
                else if (st->num_emissions == 2) col = obs2[t - 2];
                else                             col = obs3[t - 3];
                e = EmissionTable__get(st->emission_table,
                                       col,
                                       st->reference_uint);
            }
            cache[t * S + s] = e;
        }
    }
    return cache;
}
