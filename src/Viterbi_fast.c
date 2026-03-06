#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Viterbi_fast.h"
#include "EmissionCache.h"
#include "EmissionTable.h"
#include "Logodd.h"

// Helper step function isolated for clean Forward and Traceback Checkpoint regeneration
static inline void Viterbi_fast_step(
        struct HMM* hmm, 
        LOGODD_T* vmatrix, 
        STATE_ID_T* current_row, 
        LOGODD_T* emit_cache, 
        size_t t, 
        size_t S) 
{
    if (t != 0) {
        for (size_t s = 0; s < S; s++) {
            vmatrix[s * 4 + (t % 4)] = LOGODD_NEGINF;
        }
    }

    /* PASS A: Emitting States (Fully parallel) */
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < hmm->num_emitting; i++) {
        struct State* state = &hmm->states[i];
        if (state->num_emissions > t) continue;
        
        LOGODD_T emission_logodd = emit_cache[t * S + i];
        if (emission_logodd == LOGODD_NEGINF) continue;

        LOGODD_T max_logodd = LOGODD_NEGINF;
        STATE_ID_T origin_id = STATE_MAX_ID;

        for (size_t j = 0; j < state->num_incoming; j++) {
            struct Transition tr = state->incoming[j];
            size_t t_prime = t - state->num_emissions;
            LOGODD_T origin_logodd = vmatrix[tr.origin * 4 + (t_prime % 4)];
            if (origin_logodd != LOGODD_NEGINF) {
                LOGODD_T current = Logodd__add(origin_logodd, tr.logodd);
                if (current > max_logodd) {
                    max_logodd = current;
                    origin_id = tr.origin;
                }
            }
        }
        
        if (max_logodd != LOGODD_NEGINF) {
            max_logodd = Logodd__add(max_logodd, emission_logodd);
            if (max_logodd > vmatrix[i * 4 + (t % 4)]) {
                vmatrix[i * 4 + (t % 4)] = max_logodd;
                current_row[i] = origin_id;
            }
        }
    }

    /* PASS B: Silent States (Topologically sorted, strictly serial sweep) */
    for (size_t i = hmm->num_emitting; i < S; i++) {
        struct State* state = &hmm->states[i];
        // Silent states do not emit (emission cache value is typically 0, skip lookup to save ops)
        
        LOGODD_T max_logodd = LOGODD_NEGINF;
        STATE_ID_T origin_id = STATE_MAX_ID;

        for (size_t j = 0; j < state->num_incoming; j++) {
            struct Transition tr = state->incoming[j];
            LOGODD_T origin_logodd = vmatrix[tr.origin * 4 + (t % 4)];
            if (origin_logodd != LOGODD_NEGINF) {
                LOGODD_T current = Logodd__add(origin_logodd, tr.logodd);
                if (current > max_logodd) {
                    max_logodd = current;
                    origin_id = tr.origin;
                }
            }
        }
        if (max_logodd != LOGODD_NEGINF) {
            if (max_logodd > vmatrix[i * 4 + (t % 4)]) {
                vmatrix[i * 4 + (t % 4)] = max_logodd;
                current_row[i] = origin_id;
            }
        }
    }
}


void Viterbi_fast(struct HMM* hmm,
                  size_t T,
                  Literal* observations,
                  uint16_t* obs1, uint16_t* obs2, uint16_t* obs3,
                  size_t pathlength,
                  struct State** path)
{
    size_t S = hmm->num_states;
    size_t required_bytes = (T + 1) * S * sizeof(STATE_ID_T);
    bool use_checkpoints = required_bytes > (250ULL * 1024 * 1024);

    size_t interval;
#ifdef CHECKPOINT_INTERVAL
    interval = CHECKPOINT_INTERVAL;
#else
    interval = 256;
#endif

    LOGODD_T* ckpt_v = NULL;
    STATE_ID_T* ckpt_p0 = NULL;
    STATE_ID_T* dummy_row = NULL;
    STATE_ID_T* full_pmatrix = NULL;

    if (use_checkpoints) {
        size_t num_checkpoints = (T / interval) + 1;
        ckpt_v = (LOGODD_T*) malloc(num_checkpoints * 4 * S * sizeof(LOGODD_T));
        ckpt_p0 = (STATE_ID_T*) malloc(S * sizeof(STATE_ID_T));
        dummy_row = (STATE_ID_T*) malloc(S * sizeof(STATE_ID_T));
        for(size_t s = 0; s < S; s++) ckpt_p0[s] = STATE_MAX_ID;
    } else {
        full_pmatrix = (STATE_ID_T*) malloc((T + 1) * S * sizeof(STATE_ID_T));
        for(size_t _k = 0; _k < (T + 1) * S; _k++) full_pmatrix[_k] = STATE_MAX_ID;
    }

    LOGODD_T* vmatrix = (LOGODD_T*) malloc(S * 4 * sizeof(LOGODD_T));
    for (size_t s = 0; s < S * 4; s++) {
        vmatrix[s] = LOGODD_NEGINF;
    }

    LOGODD_T* emit_cache = EmissionCache__build(hmm, T, obs1, obs2, obs3);

    for (size_t i = 0; i < hmm->num_starts; i++) {
        struct Transition tr = hmm->starts[i];
        vmatrix[tr.origin * 4 + 0] = tr.logodd;
    }

    /* Formward Iteration */
    for (size_t t = 0; t <= T; t++) {
        STATE_ID_T* current_row = use_checkpoints ? dummy_row : &full_pmatrix[t * S];

        if (use_checkpoints) {
            for(size_t s = 0; s < S; s++) current_row[s] = STATE_MAX_ID;
        }
        if (t == 0) {
            for (size_t si = 0; si < hmm->num_starts; si++) {
                current_row[hmm->starts[si].origin] = hmm->starts[si].origin;
            }
        }

        Viterbi_fast_step(hmm, vmatrix, current_row, emit_cache, t, S);

        if (use_checkpoints) {
            if (t == 0) {
                memcpy(ckpt_p0, dummy_row, S * sizeof(STATE_ID_T));
            }
            if (t % interval == 0) {
                size_t c = t / interval;
                for (size_t s = 0; s < S; s++) {
                    for (int r = 0; r < 4; r++) {
                        ckpt_v[c * 4 * S + s * 4 + r] = vmatrix[s * 4 + r];
                    }
                }
            }
        }
    }

    if (use_checkpoints) free(dummy_row);

    LOGODD_T logodd = LOGODD_NEGINF;
    STATE_ID_T best_end = STATE_MAX_ID;
    for (size_t k = 0; k < hmm->num_ends; k++) {
        struct Transition end = hmm->ends[k];
        LOGODD_T current = vmatrix[end.origin * 4 + (T % 4)];
        LOGODD_T sum = Logodd__add(current, end.logodd);
        if (logodd < sum) {
            best_end = end.origin;
            logodd = sum;
        }
    }

    if (best_end == STATE_MAX_ID || logodd == LOGODD_NEGINF) {
        free(vmatrix);
        free(emit_cache);
        if (use_checkpoints) {
            free(ckpt_v);
            free(ckpt_p0);
        } else {
            free(full_pmatrix);
        }
        return;
    }

    size_t trace_buffer_size = pathlength;
    path[trace_buffer_size - 1] = &hmm->states[best_end];

    STATE_ID_T* local_pmatrix = NULL;
    LOGODD_T* vmatrix_rerun = NULL;
    long current_chunk = -1;

    if (use_checkpoints) {
        local_pmatrix = (STATE_ID_T*) malloc(interval * S * sizeof(STATE_ID_T));
        vmatrix_rerun = (LOGODD_T*) malloc(S * 4 * sizeof(LOGODD_T));
    }

    size_t i, t = T;
    for (i = trace_buffer_size - 1; i > 0; i--) {
        struct State* state = path[i];
        STATE_ID_T prev_id = STATE_MAX_ID;

        if (t == 0) {
            prev_id = use_checkpoints ? ckpt_p0[state->id] : full_pmatrix[state->id];
        } else if (!use_checkpoints) {
            prev_id = full_pmatrix[t * S + state->id];
        } else {
            long chunk_needed = (long)((t - 1) / interval);
            if (chunk_needed != current_chunk) {
                current_chunk = chunk_needed;
                
                for (size_t s = 0; s < S; s++) {
                    for (int r = 0; r < 4; r++) {
                        vmatrix_rerun[s * 4 + r] = ckpt_v[current_chunk * 4 * S + s * 4 + r];
                    }
                }
                
                size_t t_start = (size_t)(current_chunk * interval + 1);
                size_t t_end = (size_t)((current_chunk + 1) * interval);
                if (t_end > T) t_end = T;

                for (size_t t_run = t_start; t_run <= t_end; t_run++) {
                    STATE_ID_T* local_row = &local_pmatrix[((t_run - 1) % interval) * S];
                    for (size_t s = 0; s < S; s++) local_row[s] = STATE_MAX_ID;
                    Viterbi_fast_step(hmm, vmatrix_rerun, local_row, emit_cache, t_run, S);
                }
            }
            prev_id = local_pmatrix[((t - 1) % interval) * S + state->id];
        }

        if (prev_id >= hmm->num_states) break;
        path[i - 1] = &hmm->states[prev_id];

        if (t == 0 && state->num_emissions > 0) break;
        t -= state->num_emissions;
    }

    const size_t offset = i;
    for (i = 0; i < trace_buffer_size - offset; i++) {
        if (i < trace_buffer_size - offset) {
            path[i] = path[i + offset];
        } else {
            path[i] = NULL;
        }
    }

    free(emit_cache);
    free(vmatrix);
    if (use_checkpoints) {
        free(vmatrix_rerun);
        free(local_pmatrix);
        free(ckpt_v);
        free(ckpt_p0);
    } else {
        free(full_pmatrix);
    }
}
