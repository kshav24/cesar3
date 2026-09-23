#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include "cesar_worker.h"
#include "Viterbi.h"
#include "Viterbi_fast.h"
#include "Alignment.h"
#include "SafeAlloc.h"
#include "Logging.h"

#define MAX_JOBS 500000
#define MAX_PATH 4096

typedef struct {
    char* ref_aln;
    char* qry_aln;
} QueryResult;

typedef struct {
    struct HMM* hmm;
    struct Fasta* fasta;
    struct Params* params;
    QueryResult* results;
    _Atomic size_t next_q;
    size_t num_queries;
} MultiQueryContext;

static void* multi_query_worker(void* arg) {
    MultiQueryContext* ctx = (MultiQueryContext*) arg;
    while (1) {
        size_t q = __atomic_fetch_add(&ctx->next_q, 1, __ATOMIC_SEQ_CST);
        if (q >= ctx->num_queries) break;

        size_t length = 2 * ctx->fasta->queries[q]->length;
        for (uint16_t i = 0; i < ctx->fasta->num_references; i++) {
            length += ctx->fasta->references[i]->length;
        }

        struct State** path = (struct State**) SAFEMALLOC(sizeof(struct State*) * length);
        size_t path_length = length;

        Viterbi(ctx->hmm, ctx->fasta->queries[q]->length, ctx->fasta->queries[q]->sequence, &path_length, path);

        struct Alignment* alignment = Alignment__create(ctx->fasta, (uint16_t)q, ctx->params, path_length, path);
        if (alignment && alignment->reference && alignment->query) {
            ctx->results[q].ref_aln = strdup(alignment->reference);
            ctx->results[q].qry_aln = strdup(alignment->query);
        } else {
            ctx->results[q].ref_aln = NULL;
            ctx->results[q].qry_aln = NULL;
        }

        if (alignment) Alignment__destroy(alignment);
        free(path);
    }
    return NULL;
}

void cesar_align_all_queries(struct HMM* hmm, struct Fasta* fasta, struct Params* params, int num_threads) {
    if (fasta->num_queries == 0) return;

    if (num_threads <= 0) {
        long nprocs = sysconf(_SC_NPROCESSORS_ONLN);
        num_threads = (int)(0.8 * nprocs);
        if (num_threads < 1) num_threads = 1;
    }
    if (num_threads > fasta->num_queries) {
        num_threads = fasta->num_queries;
    }

    if (num_threads == 1 || fasta->num_queries == 1) {
        for (uint16_t q = 0; q < fasta->num_queries; q++) {
            size_t length = 2 * fasta->queries[q]->length;
            for (uint16_t i = 0; i < fasta->num_references; i++) {
                length += fasta->references[i]->length;
            }

            struct State** path = (struct State**) SAFEMALLOC(sizeof(struct State*) * length);
            size_t path_length = length;

            Viterbi(hmm, fasta->queries[q]->length, fasta->queries[q]->sequence, &path_length, path);

            struct Alignment* alignment = Alignment__create(fasta, q, params, path_length, path);
            if (alignment && alignment->reference && alignment->query) {
                printf(">%s\n%s\n", "referenceExon", alignment->reference);
                printf(">%s\n%s\n", fasta->queries[q]->name, alignment->query);
            }
            if (alignment) Alignment__destroy(alignment);

            free(path);
        }
        return;
    }

    QueryResult* results = (QueryResult*) SAFECALLOC(sizeof(QueryResult), fasta->num_queries);
    MultiQueryContext ctx;
    ctx.hmm = hmm;
    ctx.fasta = fasta;
    ctx.params = params;
    ctx.results = results;
    ctx.next_q = 0;
    ctx.num_queries = fasta->num_queries;

    pthread_t* threads = (pthread_t*) SAFEMALLOC(sizeof(pthread_t) * num_threads);
    for (int t = 0; t < num_threads; t++) {
        pthread_create(&threads[t], NULL, multi_query_worker, &ctx);
    }

    for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
    }
    free(threads);

    for (uint16_t q = 0; q < fasta->num_queries; q++) {
        if (results[q].ref_aln && results[q].qry_aln) {
            printf(">%s\n%s\n", "referenceExon", results[q].ref_aln);
            printf(">%s\n%s\n", fasta->queries[q]->name, results[q].qry_aln);
            free(results[q].ref_aln);
            free(results[q].qry_aln);
        }
    }
    free(results);
}

typedef struct {
    char** jobs;
    size_t num_jobs;
    _Atomic size_t current_job_idx;
    pthread_mutex_t stdout_mtx;
} WorkerPool;

static void* worker_thread(void* arg) {
    WorkerPool* pool = (WorkerPool*)arg;

    while (1) {
        size_t idx = __atomic_fetch_add(&(pool->current_job_idx), 1, __ATOMIC_SEQ_CST);
        if (idx >= pool->num_jobs) {
            break;
        }

        char* job_path = pool->jobs[idx];
        char temp_filename[256];
        snprintf(temp_filename, sizeof(temp_filename), "cesar_job_tmp_thread_%zu.out", idx);

        FILE* tmp_out = fopen(temp_filename, "w");
        if (tmp_out) {
            fprintf(tmp_out, "[Pthread Allocation] Claimed Index %zu FASTA Target: %s\n", idx, job_path);
            fclose(tmp_out);
        }

        pthread_mutex_lock(&(pool->stdout_mtx));
        FILE* tmp_in = fopen(temp_filename, "r");
        if (tmp_in) {
            char line_buf[MAX_PATH + 512];
            while (fgets(line_buf, sizeof(line_buf), tmp_in)) {
                printf("%s", line_buf);
            }
            fclose(tmp_in);
            remove(temp_filename);
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

    while (fgets(line, sizeof(line), f) && num_jobs < MAX_JOBS) {
        size_t len = strlen(line);
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
    pthread_t* threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));

    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], NULL, worker_thread, &pool);
    }
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
    pthread_mutex_destroy(&(pool.stdout_mtx));
    for (size_t i = 0; i < num_jobs; i++) {
        free(jobs[i]);
    }
    free(jobs);
}
