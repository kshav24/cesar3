#ifndef VITERBI_FAST_H
#define VITERBI_FAST_H

#include <stddef.h>
#include <stdint.h>
#include "Hmm.h"
#include "Literal.h"
#include "State.h"

void Viterbi_fast(struct HMM* hmm,
                  size_t T,
                  Literal* observations,
                  uint16_t* obs1, uint16_t* obs2, uint16_t* obs3,
                  size_t pathlength,
                  struct State** path);

#endif /* VITERBI_FAST_H */
