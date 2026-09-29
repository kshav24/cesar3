/**
 * Viterbi function declaration.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 *
 * CESAR3: Accelerated Viterbi recursion with checkpointing and state partitioning.
 */
#ifndef VITERBI_H_
#define VITERBI_H_

#include <stdbool.h>
#include "State.h"
#include "Hmm.h"
#include "Matrix.h"

struct LogoddMatrix* Viterbi__init_logodd_matrix(struct HMM* hmm);

struct PathMatrix* Viterbi__init_path_matrix(struct HMM* hmm, size_t num_observations);

LOGODD_T Viterbi__get_emission_logodd(uint16_t* obs1, uint16_t* obs2, uint16_t* obs3, const bool* has_N, Literal* observations, size_t t, struct State* state);

void Viterbi__step(struct LogoddMatrix* vmatrix, STATE_ID_T* path_row, struct HMM* hmm, Literal* observations,
                   uint16_t* obs1, uint16_t* obs2, uint16_t* obs3, const bool* has_N,
                   const STATE_ID_T* emitting_states, size_t num_emitting,
                   const STATE_ID_T* silent_states, size_t num_silent,
                   size_t t);

void Viterbi(struct HMM* hmm, size_t num_observations, Literal* observations, size_t* path_length, struct State** path);

#endif  // VITERBI_H_
