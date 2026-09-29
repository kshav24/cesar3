/**
 * Logodd operations.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <math.h>

#include "Logging.h"
#include "Logodd.h"

/**
 * Export a value to log-space. (Natural logarithm)
 * Hereby, the result will be LOGODD_NEGINF if the given value = 0.
 * @param input a value
 * @return the natural logarithm of input
 */
LOGODD_T Logodd__log(LOGODD_T input) {
  if (input == 0) {
    return LOGODD_NEGINF;
  }
  return log(input);
}

/**
 * Import a value from log-space. (Natural logarithm)
 * Hereby, the result will be LOGODD_NEGINF if the given value = LOGODD_NEGINF.
 * @param input a value
 * @return the natural logarithm of input
 */
LOGODD_T Logodd__exp(LOGODD_T input) {
  if (input == LOGODD_NEGINF) {
    return 0;
  }
  return exp(input);
}



/**
 * Compute the sum of two probabilities given in log space.
 * This method tries to avoid loosing precision due to conversion from log space.
 * @param a a log odd.
 * @param b a log odd.
 * @return log(exp(a)+exp(b)).
 */
LOGODD_T Logodd__sumexp(LOGODD_T a, LOGODD_T b) {
  if (isnan(a) || isnan(b)) {
    die("Either a, b or both are NAN.");
  }
	if (a > b) {
		if (b == LOGODD_NEGINF) {
			return a;
		}
		return a + Logodd__log(1 + Logodd__exp(b - a));
	}
	if (a == LOGODD_NEGINF) {
		return b;
	}
	return b + Logodd__log(1 + Logodd__exp(a - b));
}


/**
 * Compute the subtraction of two probabilities given in log space.
 * This method tries to avoid loosing precision due to conversion from log space.
 * @param a a log odd.
 * @param b a log odd.
 * @return log(exp(a)-exp(b)).
 */
LOGODD_T Logodd__subexp(LOGODD_T a, LOGODD_T b) {
  if (isnan(a) || isnan(b)) {
    die("Either a, b or both are NAN.");
  }
	if (a > b) {
		if (b == LOGODD_NEGINF) {
			return a;
		}
		return a + Logodd__log(1 - Logodd__exp(b - a));
	}
	if (a == LOGODD_NEGINF) {
		return b;
	}
	return b + Logodd__log(1 - Logodd__exp(a - b));
}
