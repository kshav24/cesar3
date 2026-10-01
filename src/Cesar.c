/**
 * CESAR entry point.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <libgen.h>

#include "Fasta.h"
#include "Profile.h"
#include "Logging.h"
#include "SafeAlloc.h"
#include "Viterbi.h"
#include "Alignment.h"
#include "Arguments.h"
#include "EmissionTable.h"
#include "Sequence.h"
#include "Model.h"
#include "cesar_worker.h"


char g_loglevel = LOGLEVEL;

/**
 * The main program.
 */
int main(int argc, char* argv[argc]) {
  struct EmissionTable emission_tables[8];

  struct Params parameters;
  Params__create(&parameters, emission_tables);

  if (!Arguments__read(argc, argv, &parameters)) {
    return 1;
  }

  char BaseDir[PATH_STRING_LENGTH];
  snprintf(BaseDir, PATH_STRING_LENGTH-1, "%s", dirname(argv[0]));
  logv(1, "baseDir: %s\n", BaseDir);
  Params__set_paths(&parameters, BaseDir);

  struct Fasta fasta;
  Fasta__init(&fasta);
  Fasta__read(&fasta, parameters.fasta_file);

  struct Profile** acceptors = SAFEMALLOC(sizeof(Profile*) * fasta.num_references);
  struct Profile** donors = SAFEMALLOC(sizeof(Profile*) * fasta.num_references);

  char prefix[PATH_STRING_LENGTH];
  snprintf(prefix, PATH_STRING_LENGTH-1, "%s/extra/tables/", BaseDir);
  logv(1, "prefix %s\n", prefix);

  /* load profiles for each reference exon */
  for (uint16_t i=0; i < fasta.num_references; i++) {
    struct Sequence* reference = fasta.references[i];

    /* load acceptor or first codon profile */
    if (reference->acceptor[0] != '\0') {      /* filename is given in the input file */
      char name[STATE_NAME_LENGTH];
      snprintf(name, sizeof(name), "ref%iacc", i);
      acceptors[i] = Profile__create(name);

      char fileInsideBinaryLocation[PROFILE_FILENAME_LENGTH];
      snprintf(fileInsideBinaryLocation, sizeof(fileInsideBinaryLocation), "%s%s/%s", prefix, parameters.clade, reference->acceptor);
      char pathInsideBinaryLocation[PROFILE_FILENAME_LENGTH];
      snprintf(pathInsideBinaryLocation, sizeof(pathInsideBinaryLocation), "%s/%s", BaseDir, reference->acceptor);
		logv(1,"file InsideBin %s pathInsideBin %s\n", fileInsideBinaryLocation, pathInsideBinaryLocation);

      /* test whether an absolute or relative path is given. In case of a relative path, test either local directory or the location of the binary */
      if (reference->acceptor[0] == '/') {
        Profile__read(acceptors[i], reference->acceptor);
        logv(1,"read acceptor profile for reference exon %d %s (absolute path)\n", i+1, reference->acceptor);
      }else{
        /* file inside the path where the binary is located in subdir extra/tables/clade/ */
        if (access(fileInsideBinaryLocation, R_OK) != -1) {
          Profile__read(acceptors[i], fileInsideBinaryLocation);
          logv(1,"read acceptor profile for reference exon %d %s (filename given)\n", i+1, fileInsideBinaryLocation);
        /* relative path inside the path where the binary is located */
        } else if (access(pathInsideBinaryLocation, R_OK) != -1) {
          Profile__read(acceptors[i], pathInsideBinaryLocation);
          logv(1,"read acceptor profile for reference exon %d %s (relative path given)\n", i+1, pathInsideBinaryLocation);
        } else {
		    die ("ERROR in reading the acceptor profile of exon %d: %s\n\
			   Four options:\n\
				1) Specify an absolute path\n\
				2) Specify just the filename. Then CESAR expects to find this profile in the dir where the CESAR binary is located in a subdirectory extra/tables/%s\n\
				3) Specify a relative path inside the dir where the CESAR binary is located\n", i+1, reference->acceptor, parameters.clade);
		  }
		}
    } else {
      if (fasta.num_references > 1) {
        warn("Missing acceptor profile for reference %u.", i);
      }
      acceptors[i] = Profile__create("acceptor");

		/* only read the default extra/tables/{clade}/firstCodon_profile.txt if -p {acc_profile} {do_profile} is not given (flag parameters.accSpecified) */
      if (i == 0 && !parameters.acc_do_specified && ((parameters.firstexon && fasta.num_references == 1) || (!parameters.firstexon && fasta.num_references > 1))) {
        Profile__read(acceptors[i], parameters.first_codon_profile);
        logv(1,"read acceptor profile for reference exon %d %s\n", i+1, parameters.first_codon_profile);
      } else {
        Profile__read(acceptors[i], parameters.acc_profile);
        logv(1,"read acceptor profile for reference exon %d %s\n", i+1, parameters.acc_profile);
      }
    }
    logv(1, "Reference %u uses acceptor:\t%s", i, acceptors[i]->filename);

    /* load donor or last codon profile */
    if (reference->donor[0] != '\0') {        /* filename is given in the input file */
      char name[STATE_NAME_LENGTH];
      snprintf(name, sizeof(name), "ref%idon", i);
      donors[i] = Profile__create(name);

      char fileInsideBinaryLocation[PROFILE_FILENAME_LENGTH];
      snprintf(fileInsideBinaryLocation, sizeof(fileInsideBinaryLocation), "%s%s/%s", prefix, parameters.clade, reference->donor);
      char pathInsideBinaryLocation[PROFILE_FILENAME_LENGTH];
      snprintf(pathInsideBinaryLocation, sizeof(pathInsideBinaryLocation), "%s/%s", BaseDir, reference->donor);
		logv(1,"file InsideBin %s pathInsideBin %s\n", fileInsideBinaryLocation, pathInsideBinaryLocation);

      /* test whether an absolute or relative path is given. In case of a relative path, test either local directory or the location of the binary */
      if (reference->donor[0] == '/') {
        Profile__read(donors[i], reference->donor);
        logv(1,"read donor profile for reference exon %d %s (absolute path)\n", i+1, reference->donor);
      }else{
        /* file inside the path where the binary is located in subdir extra/tables/clade/ */
        if (access(fileInsideBinaryLocation, R_OK) != -1) {
          Profile__read(donors[i], fileInsideBinaryLocation);
          logv(1,"read donor profile for reference exon %d %s (filename given)\n", i+1, fileInsideBinaryLocation);
        /* relative path inside the path where the binary is located */
        } else if (access(pathInsideBinaryLocation, R_OK) != -1) {
          Profile__read(donors[i], pathInsideBinaryLocation);
          logv(1,"read donor profile for reference exon %d %s (relative path given)\n", i+1, pathInsideBinaryLocation);
        } else {
		    die ("ERROR in reading the donor profile of exon %d: %s\n\
			   Four options:\n\
				1) Specify an absolute path\n\
				2) Specify just the filename. Then CESAR expects to find this profile in the dir where the CESAR binary is located in a subdirectory extra/tables/%s\n\
				3) Specify a relative path inside the dir where the CESAR binary is located\n", i+1, reference->donor, parameters.clade);
		  }
		}
    } else {
      if (fasta.num_references > 1) {
        warn("Missing donor profile for reference %u.", i);
      }
      donors[i] = Profile__create("donor");
		/* only read the default extra/tables/{clade}/firstCodon_profile.txt if -p {acc_profile} {do_profile} is not given (flag parameters.doSpecified) */
      if (i == fasta.num_references-1 && !parameters.acc_do_specified && ((parameters.lastexon && fasta.num_references == 1) || (!parameters.lastexon && fasta.num_references > 1))) {
        Profile__read(donors[i], parameters.last_codon_profile);
        logv(1,"read donor profile for reference exon %d %s\n", i+1, parameters.last_codon_profile);
      } else {
        Profile__read(donors[i], parameters.do_profile);
        logv(1,"read donor profile for reference exon %d %s\n", i+1, parameters.do_profile);
      }
    }
    logv(1, "Reference %u uses donor:\t%s", i, donors[i]->filename);
  }

  /* estimate memory consumption */
  size_t num_states = 0;
  size_t rlength = 0;
  size_t qlength = 0;
  size_t qlength_max = 0;
  for (uint16_t i=0; i < fasta.num_references; i++) {
    struct Sequence* reference = fasta.references[i];
    num_states += 6 + 6 * reference->num_codons + 1 + 2 + 2 + 22 + 6;  /* 22 and 6 for acc and donor states */

    logv(1, "Reference %u length: %lu", i, fasta.references[i]->length);
    logv(1, "Reference %u split codon lengths: %u %u", i, fasta.references[i]->start_split_length, fasta.references[i]->end_split_length);
    rlength += fasta.references[i]->length;
  }
  for (uint16_t i=0; i < fasta.num_queries; i++) {
    logv(1, "Query %u length: %lu", i, fasta.queries[i]->length);
    qlength += fasta.queries[i]->length;
    if (fasta.queries[i]->length > qlength_max) {
	    qlength_max = fasta.queries[i]->length;
    }
  }

#ifdef CHECKPOINT_INTERVAL
  size_t interval = CHECKPOINT_INTERVAL;
#else
  size_t interval = 256;
#endif
  size_t num_checkpoints = (qlength_max / interval) + 1;

  /* Memory that exists once per process no matter how many queries run at the
   * same time: the HMM itself, the parsed sequences and the small tables. */
  double mem_shared =
     (num_states * sizeof(State))                          +   /* Space for all the state structures, each 304 bytes for a struct State */
	  ((qlength_max + rlength) * 1)                         +   /* Space to store all the sequences */
	  3000000;                                                  /* need some smaller emission tables (32 k) and other small stuff. 
                                                                       Tests have shown that we may underestimate the real memory by up to 2 Mb. Therefore, be conservative. */

  /* Memory that every concurrent worker needs for itself: each worker runs its
   * own Viterbi, with its own vmatrix, checkpoints, traceback chunk and path. */
  double mem_per_worker =
     (num_states * 4 * sizeof(double))                     +   /* Viterbi logodds vmatrix, double has typically 8 bytes */
     (num_checkpoints * 4 * num_states * sizeof(double))   +   /* Viterbi checkpoint logodds slices */
     (interval * num_states * sizeof(uint32_t))            +   /* Viterbi local traceback chunk */
	  ((2 * qlength_max + rlength) * sizeof(struct State*)) +   /* Space for Viterbi path, 8 bytes for a pointer on a 64 bit system */
	  ((qlength_max + rlength) * 2 * sizeof(char));             /* Space for the aligned seqs */

  /* The workers are real threads inside this process and really do hold all of
   * that simultaneously, so the check has to account for them. Use the same
   * helper that cesar_align_all_queries() uses, so the two cannot disagree. */
  int num_workers = cesar_worker__effective_threads(parameters.num_threads, fasta.num_queries);
  double mem = (mem_shared + mem_per_worker * (double)num_workers) / 1000000000.0;  /* in GB */
  logv(1, "Expecting a memory consumption of: %f GB (%d concurrent worker(s), max_memory %f)", mem, num_workers, (double)parameters.max_memory);
  if (mem > (double)parameters.max_memory) {
    die("The memory consumption is limited to %1.4f GB. Your attempt requires %1.4f GB (%d concurrent worker(s) needing %1.4f GB each, plus %1.4f GB shared).\nYou can change the limit via --max-memory, or cap each process with -x when running several species in parallel.",
        (double)parameters.max_memory, mem, num_workers, mem_per_worker/1000000000.0, mem_shared/1000000000.0);
  } else {
    logv(1, "Expecting a memory consumption of: %1.4f GB (%d concurrent worker(s))", mem, num_workers);
  }


  if (g_loglevel >= 7) {
    char* tmp;
    for (uint16_t i=0; i < fasta.num_references; i++) {
      tmp = SAFECALLOC(sizeof(char), SEQUENCENAMELENGTH + fasta.references[i]->length);
      Literal__str(fasta.references[i]->length, fasta.references[i]->sequence, tmp);
      logv(1, ">original %s\n%s", fasta.references[i]->name, tmp);
      free(tmp);
    }

    for (uint16_t i=0; i < fasta.num_queries; i++) {
      tmp = SAFECALLOC(sizeof(char), SEQUENCENAMELENGTH + fasta.queries[i]->length);
      Literal__str(fasta.queries[i]->length, fasta.queries[i]->sequence, tmp);
      logv(1, ">original %s\n%s", fasta.queries[i]->name, tmp);
      free(tmp);
    }
  }

  clock_t time = clock();
  struct HMM* hmm = multi_exon(&parameters, &fasta, acceptors, donors);
  logv(1, "HMM construction (sec):\t%f", (float)(clock() - time) / CLOCKS_PER_SEC);

  if (parameters.dot[0] != '\0') {
    FILE * dotfile = fopen(parameters.dot, "w");
    HMM__dot(hmm, dotfile);
    fclose(dotfile);
  }

  cesar_align_all_queries(hmm, &fasta, &parameters, parameters.num_threads);

  HMM__destroy(hmm);
  for (uint16_t i=0; i<fasta.num_references; i++) {
    if(acceptors[i] != NULL) {
      Profile__destroy(acceptors[i]);
    }
    if(donors[i] != NULL) {
      Profile__destroy(donors[i]);
    }
  }
  free(acceptors);
  free(donors);
  Fasta__destroy(&fasta);
  Params__destroy(&parameters);
  return 0;
}
