/**
 * Logodd definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef LOGODD_H_
#define LOGODD_H_

#include <stddef.h>
#include <stdint.h>
#include <float.h>

//#define LOGODD_T long double
//#define LOGODD_NEGINF -LDBL_MAX
#define LOGODD_T double
#define LOGODD_NEGINF -DBL_MAX

#include <math.h>

#ifdef DEBUG
#include "Logging.h"
#endif

LOGODD_T Logodd__log(LOGODD_T input);
LOGODD_T Logodd__exp(LOGODD_T input);

static inline LOGODD_T Logodd__add(LOGODD_T a, LOGODD_T b) {
#ifdef DEBUG
  if (isnan(a)) {
    die("a is NAN");
  } else if (isnan(b)) {
    die("b is NAN");
  }
#endif
  if (a == LOGODD_NEGINF || b == LOGODD_NEGINF) {
    return LOGODD_NEGINF;
  }
  return a + b;
}

LOGODD_T Logodd__sumexp(LOGODD_T a, LOGODD_T b);
LOGODD_T Logodd__subexp(LOGODD_T a, LOGODD_T b);

#endif  // LOGODD_H_
