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

/*
 * Batch file runner for a list of FASTA files.
 */
void cesar_worker_run(const char* jobs_file, int num_threads);

#ifdef __cplusplus
}
#endif

#endif /* CESAR_WORKER_H */
