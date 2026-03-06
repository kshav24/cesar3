#ifndef HMM_TOPO_H
#define HMM_TOPO_H

#include "Hmm.h"

/**
 * Topologically sort HMM states so that:
 *   states[0 .. num_emitting-1]  = emitting states (parallelisable)
 *   states[num_emitting .. S-1]  = silent states in dependency order
 * Sets hmm->num_emitting.
 * Remaps all origin IDs in: incoming[], starts[], ends[].
 */
void HMM_topo__sort(struct HMM* hmm);

#endif /* HMM_TOPO_H */
