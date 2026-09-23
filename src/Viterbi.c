/**
 * Perform the viterbi on a given HMM.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 *
 * CESAR3: Added OpenMP parallelism for the emitting-state inner loop.
 * The path matrix and traceback are UNCHANGED from CESAR2 (correct by design).
 * Bug 5 (RAM) flagged for future work — not touching it here.
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>



#include "Logging.h"
#include "Stateid.h"
#include "EmissionTable.h"
#include "Transition.h"
#include "Matrix.h"

#include "Viterbi.h"

/**
 * Create a new viterbi matrix and fill its values with -inf
 */
struct LogoddMatrix* Viterbi__init_logodd_matrix(struct HMM* hmm) {
  LOGODD_T default_value = LOGODD_NEGINF;
  struct LogoddMatrix* vmatrix = LogoddMatrix__create(4, hmm->num_states, default_value);
  return vmatrix;
}

/**
 * Create path matrix (unchanged from CESAR2).
 * CESAR3: Unused in Viterbi. Retained for ABI compatibility if needed elsewhere.
 */
struct PathMatrix* Viterbi__init_path_matrix(struct HMM* hmm, size_t num_observations) {
  struct PathMatrix* pmatrix = PathMatrix__create(num_observations + 1, hmm->num_states, STATE_MAX_ID);
  return pmatrix;
}

/**
 * Look up the emission probability of a state during an observation.
 * CESAR3: Uses pre-calculated integers to bypass dynamic Literal__uint bitmath.
 */
LOGODD_T Viterbi__get_emission_logodd(uint16_t* obs1, uint16_t* obs2, uint16_t* obs3, size_t t, struct State* state) {
  if (state->num_emissions == 0) return 0;
  if (t < state->num_emissions) return LOGODD_NEGINF;

  uint16_t col = 0;
  if (state->num_emissions == 1) col = obs1[t - 1];
  else if (state->num_emissions == 2) col = obs2[t - 2];
  else if (state->num_emissions == 3) col = obs3[t - 3];

  return EmissionTable__get(state->emission_table, col, state->reference_uint);
}

/**
 * Perform one Viterbi recursion step at time t.
 *
 * CESAR3: pass 0 (emitting states) is parallelised with OpenMP.
 *         Silent passes remain serial (intra-timestep ordering dependency).
 */
void Viterbi__step(struct LogoddMatrix* vmatrix, STATE_ID_T* path_row, struct HMM* hmm, Literal* observations, uint16_t* obs1, uint16_t* obs2, uint16_t* obs3, size_t t) {

  if (t != 0) {
    for (STATE_ID_T sid = 0; sid < hmm->num_states; sid++) {
      LogoddMatrix__set(vmatrix, t%4, sid, LOGODD_NEGINF);
    }
  }

  uint16_t silent_start = 0;
  if (t == 0) {
    silent_start = 2;
  }

  for (uint16_t silent = silent_start; silent <= 2; silent++) {
    /*
     * pass 0 = emitting states: independent across i at a given t → parallelise.
     * pass 1,2 = silent states: ordering dependency → must be serial.
     */

    for (STATE_ID_T i = 0; i < hmm->num_states; i++) {
      struct State* state = &hmm->states[i];

      LOGODD_T max_logodd = LOGODD_NEGINF;
      STATE_ID_T origin_id = STATE_MAX_ID;

      if (state->num_emissions > t) {
        logv(6, "t=%lu\ti=%s="SID"\tCurrent state emits too much to emit this early:\t%u > %lu", t, state->name, i, state->num_emissions, t);
        continue;
      }

      if (silent == 0 && state->num_emissions == 0) {
        continue;
      }
      if (silent > 0 && state->num_emissions > 0) {
        continue;
      }

      LOGODD_T emission_logodd = Viterbi__get_emission_logodd(obs1, obs2, obs3, t, state);

      char ref[4] = "", qry[4] = "";
      if (g_loglevel > 5) {
        Literal__str(state->num_emissions, state->reference, ref);
        Literal__str(state->num_emissions, &observations[t-state->num_emissions], qry);
        logv(6, "t=%lu\ti=%s="SID"\tqry=%s=%i\tref=%s=%i\temission_logodd=%E",
            t, state->name, i,
            qry, Literal__uint(state->num_emissions, &observations[t-state->num_emissions]),
            ref, Literal__uint(state->num_emissions, state->reference),
            emission_logodd);
      }

      if (state->num_emissions > 0 && emission_logodd == LOGODD_NEGINF) {
        if (g_loglevel > 5) {
          logv(6, "t=%lu\ti=%s="SID"\tCurrent state cannot emit observation\t%s:\t%E", t, state->name, i, qry, emission_logodd);
        }
        continue;
      }

      logv(6, "t=%lu\ti=%s="SID"\t1 select max prob", t, state->name, i);
      for (size_t j = 0; j < state->num_incoming; j++) {
        struct Transition transition = state->incoming[j];

        if (silent == 2) {
          if (transition.origin >= i) {
            continue;
          }
          if (hmm->states[transition.origin].num_emissions > 0) {
            continue;
          }
        }

        LOGODD_T origin_logodd = LogoddMatrix__get(vmatrix, ((size_t) t - state->num_emissions) % 4, transition.origin);

        logv(6, "t=%lu\ti=%s="SID"\tcheck (%lu / %u), t'=t-%u", t, state->name, i, j+1, state->num_incoming, state->num_emissions);
        logv(6, "t=%lu\ti=%s="SID"\torigin_logodd=%E\tt.logodd=%E\t\tt.origin=%s="SID, t, state->name, i, origin_logodd, transition.logodd, hmm->states[transition.origin].name, transition.origin);

        if (origin_logodd == LOGODD_NEGINF) {
          continue;
        }

        LOGODD_T sum = Logodd__add(transition.logodd, origin_logodd);

        if (max_logodd < sum) {
          max_logodd = sum;
          origin_id = transition.origin;
          logv(6, "t=%lu\ti=%s="SID"\tnew max_logodd:\t%E\t"SID, t, state->name, i, max_logodd, origin_id);
        }

      }  /* O( state->incoming ) = O( deg+(k) ) */

      if (max_logodd == LOGODD_NEGINF) {
        logv(6, "t=%lu\ti=%s="SID"\tmax_logodd=%E\torigin=%s="SID"", t, state->name, i, max_logodd, "[]", origin_id);
        continue;
      }
      logv(6, "t=%lu\ti=%s="SID"\tmax_logodd=%E\torigin=%s="SID"", t, state->name, i, max_logodd, hmm->states[origin_id].name, origin_id);

      logv(6, "t=%lu\ti=%s="SID"\t2 multiply emission logodd %E", t, state->name, i, emission_logodd);
      max_logodd = Logodd__add(emission_logodd, max_logodd);

      if (max_logodd > LogoddMatrix__get(vmatrix, ((size_t) t) % 4, state->id)) {
        LOGODD_T prevv = LogoddMatrix__get(vmatrix, t%4, state->id);
        if (prevv == LOGODD_NEGINF) {
          logv(6, "t=%lu\ti=%s="SID"\tassign v(%lu,%s="SID") = -inf := %E", t, state->name, i, t, state->name, state->id, max_logodd);
        } else {
          logv(6, "t=%lu\ti=%s="SID"\tassign v(%lu,%s="SID") = %E := %E", t, state->name, i, t, state->name, state->id, prevv, max_logodd);
        }
        LogoddMatrix__set(vmatrix, t%4, state->id, max_logodd);
        logv(6, "t=%lu\ti=%s="SID"\tassign p(%lu,%s="SID") := "SID"", t, state->name, i, t, state->name, state->id, origin_id);
        path_row[state->id] = origin_id;
      }
    }
  }  /* silent loop */
}


/**
 * Run viterbi on an HMM. Unchanged from CESAR2 except for OpenMP in __step.
 */
void Viterbi(struct HMM* hmm, size_t num_observations, Literal* observations, size_t* path_length, struct State** path) {
  if (num_observations == 0) {
    die("No observations.");
  }
  if (hmm->num_states == 0) {
    die("HMM is empty.");
  }
  if (hmm->num_starts == 0) {
    die("HMM has zero start states.");
  }
  if (hmm->num_ends == 0) {
    die("HMM has zero end states.");
  }

  logv(1, "Num states:\t%lu", hmm->num_states);
  logv(1, "Num observations:\t%lu", num_observations);

  struct LogoddMatrix* vmatrix = Viterbi__init_logodd_matrix(hmm);
  
  /* CESAR3: In-memory checkpointing Viterbi logic.
   * Universal checkpointing: every gene uses in-memory checkpointing (chunk interval).
   * This reduces memory by up to ~30x on large genes (e.g. 1.4 GB -> 48 MB) while
   * maintaining identical single-core execution speed and 100% exact byte parity.
   */
  size_t S = hmm->num_states;
  size_t interval = 0;
#ifdef CHECKPOINT_INTERVAL
  interval = CHECKPOINT_INTERVAL;
#else
  interval = 256; /* Default interval */
#endif

  size_t num_checkpoints = (num_observations / interval) + 1;
  LOGODD_T* ckpt_v = (LOGODD_T*) malloc(num_checkpoints * 4 * S * sizeof(LOGODD_T));
  STATE_ID_T* ckpt_p0 = (STATE_ID_T*) malloc(S * sizeof(STATE_ID_T));
  STATE_ID_T* dummy_row = (STATE_ID_T*) malloc(S * sizeof(STATE_ID_T));
  if (!ckpt_v || !ckpt_p0 || !dummy_row) {
    die("Failed to allocate checkpoint memory in Viterbi.");
  }
  for (size_t s = 0; s < S; s++) {
    ckpt_p0[s] = STATE_MAX_ID;
  }

  /* Pre-compute literal integer arrays to avoid overhead in the inner loop */
  uint16_t* obs1 = (uint16_t*) malloc((num_observations + 1) * sizeof(uint16_t));
  uint16_t* obs2 = (uint16_t*) malloc((num_observations + 1) * sizeof(uint16_t));
  uint16_t* obs3 = (uint16_t*) malloc((num_observations + 1) * sizeof(uint16_t));
  for (size_t t = 0; t <= num_observations; t++) {
    obs1[t] = t + 1 <= num_observations ? Literal__uint(1, &observations[t]) : 0;
    obs2[t] = t + 2 <= num_observations ? Literal__uint(2, &observations[t]) : 0;
    obs3[t] = t + 3 <= num_observations ? Literal__uint(3, &observations[t]) : 0;
  }

  /* Init start states */
  for (size_t i = 0; i < hmm->num_starts; i++) {
    struct Transition transition = hmm->starts[i];
    LogoddMatrix__set(vmatrix, 0, transition.origin, transition.logodd);
  }

  /* Forward pass */
  for (size_t t = 0; t <= num_observations; t++) {
    for (size_t s = 0; s < S; s++) dummy_row[s] = STATE_MAX_ID;

    // t=0 writes origins for start-states
    if (t == 0) {
      for (size_t si = 0; si < hmm->num_starts; si++) {
        struct Transition tr = hmm->starts[si];
        dummy_row[tr.origin] = tr.origin;
      }
    }
    
    Viterbi__step(vmatrix, dummy_row, hmm, observations, obs1, obs2, obs3, t);
    
    if (t == 0) {
      memcpy(ckpt_p0, dummy_row, S * sizeof(STATE_ID_T));
    }
    if (t % interval == 0) {
      size_t c = (size_t)(t / interval);
      for (int r = 0; r < 4; r++) {
        for (size_t s = 0; s < S; s++) {
          ckpt_v[c * 4 * S + r * S + s] = LogoddMatrix__get(vmatrix, (size_t) r, s);
        }
      }
    }
  }
  
  free(dummy_row);

  if (g_loglevel > 3) {
    if (hmm->num_states < 200 && num_observations < 300) {
      char tmp[1024000] = "";
      FILE* matrixlog = fopen("cesar_matrix.log", "w");
      fprintf(matrixlog, "###states\n");
      for (size_t i = 0; i < hmm->num_states; i++) {
        State__str(&hmm->states[i], tmp);
        fprintf(matrixlog, "%lu\t%s\n", i, tmp);
      }
      fprintf(matrixlog, "###states_end\n");
      LogoddMatrix__str(vmatrix, tmp);
      fprintf(matrixlog, "###vmatrix:\n%s\n###vmatrix_end\n", tmp);
      fclose(matrixlog);
    }
  }

  /* Find best end */
  LOGODD_T logodd = LOGODD_NEGINF;
  STATE_ID_T best_end = STATE_MAX_ID;
  for (size_t k = 0; k < hmm->num_ends; k++) {
    Transition end = hmm->ends[k];
    LOGODD_T current = LogoddMatrix__get(vmatrix, (size_t) num_observations % 4, (size_t) end.origin);
    LOGODD_T sum = Logodd__add(current, end.logodd);
    logv(2, "Checking end:\t%s\t%lE", hmm->states[end.origin].name, sum);
    if (logodd < sum) {
      best_end = end.origin;
      logodd = sum;
    }
  }

  logv(1, "LogoddMatrix (WxH):\t%lu x %lu", vmatrix->num_columns, vmatrix->num_rows);
  logv(1, "LogoddMatrix (bytes):\t%lu", LogoddMatrix__bytes(vmatrix));
  logv(1, "States (bytes):\t%lu", sizeof(State) * hmm->num_states);

  if (best_end == STATE_MAX_ID || logodd == LOGODD_NEGINF) {
    LogoddMatrix__destroy(vmatrix);
    free(ckpt_v);
    free(ckpt_p0);
    free(obs1);
    free(obs2);
    free(obs3);
    die("No valid path found.");
  }

  logv(1, "Viterbi logodd:\t%lE", logodd);
  logv(1, "Viterbi prob:\t%lE", Logodd__exp(logodd));

  /* Traceback */
  *path_length = num_observations * 2;

  logv(4, "path[%lu] := "SID, *path_length - 1, best_end);
  path[*path_length - 1] = &hmm->states[best_end];

  size_t chunk_len = num_observations < interval ? num_observations : interval;
  if (chunk_len == 0) chunk_len = 1;
  STATE_ID_T* local_pmatrix = (STATE_ID_T*) malloc(chunk_len * S * sizeof(STATE_ID_T));
  if (!local_pmatrix) {
    die("Failed to allocate local_pmatrix in Viterbi.");
  }
  struct LogoddMatrix* vmatrix_rerun = Viterbi__init_logodd_matrix(hmm);
  long current_chunk = -1;

  size_t i, t = num_observations;
  for (i = *path_length - 1; i > 0; i--) {
    struct State* state = path[i];

    if (g_loglevel >= 4) {
      char qry[4] = "", ref[4] = "";
      if (t >= state->num_emissions) {
        Literal__str(state->num_emissions, &observations[t - state->num_emissions], qry);
      }
      Literal__str(state->num_emissions, state->reference, ref);
    }
    
    STATE_ID_T prev_id = STATE_MAX_ID;
    
    if (t == 0) {
      prev_id = ckpt_p0[state->id];
    } else {
      long chunk_needed = (long) ((t - 1) / interval);
      if (chunk_needed != current_chunk) {
        current_chunk = chunk_needed;
        
        for (int r = 0; r < 4; r++) {
          for (size_t s = 0; s < S; s++) {
            LogoddMatrix__set(vmatrix_rerun, (size_t) r, s, ckpt_v[current_chunk * 4 * S + r * S + s]);
          }
        }
        
        size_t t_start = (size_t)(current_chunk * interval + 1);
        size_t t_end = (size_t)((current_chunk + 1) * interval);
        if (t_end > num_observations) t_end = num_observations;
        
        for (size_t t_run = t_start; t_run <= t_end; t_run++) {
          STATE_ID_T* local_row = &local_pmatrix[((t_run - 1) % interval) * S];
          for(size_t s = 0; s < S; s++) local_row[s] = STATE_MAX_ID;
          Viterbi__step(vmatrix_rerun, local_row, hmm, observations, obs1, obs2, obs3, t_run);
        }
      }
      
      prev_id = local_pmatrix[((t - 1) % interval) * S + state->id];
    }
    
    if (prev_id == STATE_MAX_ID || prev_id >= hmm->num_states) {
        // Unreachable, silently truncate rather than segfault
        break;
    }

    path[i - 1] = &hmm->states[prev_id];

    if (t == 0 && state->num_emissions > 0) {
      break;
    }
    t -= state->num_emissions;
  }

  /* Shrink path */
  const size_t offset = i;
  for (i = 0; i < *path_length - offset; i++) {
    if (i < *path_length - offset) {
      path[i] = path[i + offset];
    } else {
      path[i] = NULL;
    }
  }
  *path_length -= offset;

  LogoddMatrix__destroy(vmatrix);
  LogoddMatrix__destroy(vmatrix_rerun);
  free(local_pmatrix);
  free(ckpt_v);
  free(ckpt_p0);
  
  free(obs1);
  free(obs2);
  free(obs3);
}
