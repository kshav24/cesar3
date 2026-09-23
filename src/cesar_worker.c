#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "cesar_worker.h"

#define MAX_JOBS 500000
#define MAX_PATH 4096

typedef struct {
    char** jobs;
    size_t num_jobs;
    _Atomic size_t current_job_idx;
    pthread_mutex_t stdout_mtx;
} WorkerPool;

static void* worker_thread(void* arg) {
    WorkerPool* pool = (WorkerPool*)arg;

    while (1) {
        // Atomically claim the next job
        size_t idx = __atomic_fetch_add(&(pool->current_job_idx), 1, __ATOMIC_SEQ_CST);
        
        // If the pool is exhausted, the thread can safely terminate
        if (idx >= pool->num_jobs) {
            break;
        }

        char* job_path = pool->jobs[idx];
        
        // Stub: In Step 6, this boundary invokes Viterbi_fast() on the loaded FASTA.
        // For now, it isolates the output into a strictly sequential temporary chunk buffer.
        char temp_filename[256];
        snprintf(temp_filename, sizeof(temp_filename), "cesar_job_tmp_thread_%zu.out", idx);
        
        FILE* tmp_out = fopen(temp_filename, "w");
        if (tmp_out) {
            fprintf(tmp_out, "[Pthread Allocation] Claimed Index %zu FASTA Target: %s\n", idx, job_path);
            fclose(tmp_out);
        }

        // Lock the global standard output to guarantee isolated continuous job logging
        pthread_mutex_lock(&(pool->stdout_mtx));
        
        FILE* tmp_in = fopen(temp_filename, "r");
        if (tmp_in) {
            char line_buf[MAX_PATH + 512];
            while (fgets(line_buf, sizeof(line_buf), tmp_in)) {
                printf("%s", line_buf);
            }
            fclose(tmp_in);
            remove(temp_filename); // Instantly clean residue
        }
        
        pthread_mutex_unlock(&(pool->stdout_mtx));
    }
    
    return NULL;
}

void cesar_worker_run(const char* jobs_file, int num_threads) {
    FILE* f = fopen(jobs_file, "r");
    if (!f) {
        fprintf(stderr, "Fatal Error: Unable to open the target jobs file segment: %s\n", jobs_file);
        return;
    }

    char** jobs = (char**) malloc(MAX_JOBS * sizeof(char*));
    size_t num_jobs = 0;
    char line[MAX_PATH];

    // Build the dynamic job list
    while (fgets(line, sizeof(line), f) && num_jobs < MAX_JOBS) {
        size_t len = strlen(line);
        // Clean trailing newline residues
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[len-1] = '\0';
            len--;
        }
        if (len > 0) {
            jobs[num_jobs] = strdup(line);
            num_jobs++;
        }
    }
    fclose(f);

    if (num_jobs == 0) {
        free(jobs);
        return;
    }

    WorkerPool pool;
    pool.jobs = jobs;
    pool.num_jobs = num_jobs;
    pool.current_job_idx = 0;
    pthread_mutex_init(&(pool.stdout_mtx), NULL);

    if (num_threads <= 0) num_threads = 1;

    // Spin up the configured threads pool
    pthread_t* threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], NULL, worker_thread, &pool);
    }

    // Await flawless atomic termination
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Teardown
    pthread_mutex_destroy(&(pool.stdout_mtx));
    for (size_t i = 0; i < num_jobs; i++) {
        free(jobs[i]);
    }
    free(jobs);
    free(threads);
}
