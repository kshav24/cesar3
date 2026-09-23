#include <stdlib.h>
#include <string.h>
#include "HMM_topo.h"

void HMM_topo__sort(struct HMM* hmm) {
    size_t S = hmm->num_states;
    STATE_ID_T* order   = malloc(S * sizeof(STATE_ID_T));
    STATE_ID_T* remap   = malloc(S * sizeof(STATE_ID_T));
    uint8_t*    placed  = calloc(S, 1);
    size_t      head    = 0;

    /* Pass 1: all emitting states first — fully independent */
    for (STATE_ID_T i = 0; i < (STATE_ID_T)S; i++) {
        if (hmm->states[i].num_emissions > 0) {
            order[head++] = i;
            placed[i] = 1;
        }
    }
    hmm->num_emitting = head;

    /* Pass 2: silent states in topological order (Kahn's algorithm) */
    size_t* indegree = calloc(S, sizeof(size_t));
    for (STATE_ID_T i = 0; i < (STATE_ID_T)S; i++) {
        if (placed[i]) continue;
        for (uint16_t j = 0; j < hmm->states[i].num_incoming; j++) {
            STATE_ID_T src = hmm->states[i].incoming[j].origin;
            if (hmm->states[src].num_emissions == 0)
                indegree[i]++;
        }
    }

    STATE_ID_T* queue = malloc(S * sizeof(STATE_ID_T));
    size_t qhead = 0, qtail = 0;
    for (STATE_ID_T i = 0; i < (STATE_ID_T)S; i++) {
        if (!placed[i] && indegree[i] == 0)
            queue[qtail++] = i;
    }
    while (qhead < qtail) {
        STATE_ID_T u = queue[qhead++];
        order[head++] = u;
        placed[u] = 1;
        for (STATE_ID_T k = 0; k < (STATE_ID_T)S; k++) {
            if (placed[k] || hmm->states[k].num_emissions > 0) continue;
            for (uint16_t j = 0; j < hmm->states[k].num_incoming; j++) {
                if (hmm->states[k].incoming[j].origin == u) {
                    if (--indegree[k] == 0)
                        queue[qtail++] = k;
                }
            }
        }
    }

    /* Fallback: if any silent state was unreachable, place it to complete the bijection */
    if (head < S) {
        for (STATE_ID_T i = 0; i < (STATE_ID_T)S; i++) {
            if (!placed[i]) {
                order[head++] = i;
                placed[i] = 1;
            }
        }
    }

    /* Build remap: old_id -> new position */
    for (size_t i = 0; i < S; i++)
        remap[order[i]] = (STATE_ID_T)i;

    /* Reorder states array in-place */
    struct State* tmp = malloc(S * sizeof(struct State));
    for (size_t i = 0; i < S; i++)
        tmp[i] = hmm->states[order[i]];
    memcpy(hmm->states, tmp, S * sizeof(struct State));
    free(tmp);

    /* Fix state IDs and incoming transition origins */
    for (size_t i = 0; i < S; i++) {
        hmm->states[i].id = (STATE_ID_T)i;
        for (uint16_t j = 0; j < hmm->states[i].num_incoming; j++)
            hmm->states[i].incoming[j].origin =
                remap[hmm->states[i].incoming[j].origin];
    }

    /* Fix starts[] and ends[] — CRITICAL: these also hold origin IDs */
    for (size_t i = 0; i < hmm->num_starts; i++)
        hmm->starts[i].origin = remap[hmm->starts[i].origin];
    for (size_t i = 0; i < hmm->num_ends; i++)
        hmm->ends[i].origin = remap[hmm->ends[i].origin];

    free(order); free(remap); free(placed);
    free(indegree); free(queue);
}
