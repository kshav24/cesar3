#ifndef CESAR_WORKER_H
#define CESAR_WORKER_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Reads a text file of FASTA paths (one per line) and distributes them
 * across a dynamically sized pthread pool utilizing an atomic fetch-add
 * index counter. Prints isolated output segments securely under mutexes.
 */
void cesar_worker_run(const char* jobs_file, int num_threads);

#ifdef __cplusplus
}
#endif

#endif /* CESAR_WORKER_H */
