#ifndef CESAR_WORKER_H
#define CESAR_WORKER_H

#include "Hmm.h"
#include "Fasta.h"
#include "Params.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Distribute sequence alignments across multiple query species
 * using an in-process pthread worker pool with dynamic work-stealing.
 * Results are printed in strictly sequential FASTA query order to stdout.
 */
void cesar_align_all_queries(struct HMM* hmm, struct Fasta* fasta, struct Params* params, int num_threads);

/**
 * Number of query workers that will actually run concurrently.
 * Shared by cesar_align_all_queries() and the memory pre-check in main(), so
 * the two can never disagree about how many Viterbi instances exist at once.
 * @param num_threads requested threads (<= 0 means auto).
 * @param num_queries number of queries in the input.
 * @return effective worker count (>= 1).
 */
int cesar_worker__effective_threads(int num_threads, uint16_t num_queries);

#ifdef __cplusplus
}
#endif

#endif /* CESAR_WORKER_H */
