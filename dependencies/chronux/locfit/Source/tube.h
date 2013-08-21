/*
 * Copyright 1996-2006 Catherine Loader.
 */
/*
 *   Copyright (c) 1998-2006 Catherine Loader
 *   See README file for details.
 *
 *
 *   Headers for the tube library.
 */

#ifndef I_TUBE_H
#define I_TUBE_H

/*
 * public functions needed by routines calling the tube library.
 */
extern double critval();
extern double tailp(), taild();
extern int tube_constants();
extern int k0_reqd();

/*
 * stuff used internally.
 */

#include "stdlib.h"
#include "mut.h"

#define TUBE_MXDIM 10

/*
 * definitions for integration methods.
 * these match locfit evaluation structures where applicable.
 */

#define ISIMPSON  4    /* grid */
#define ISPHERIC 11    /* circle or sphere */
#define IDERFREE 25    /* derivative free */
#define IMONTE   30    /* monte carlo */

#ifndef PI
#define PI    3.141592653589793238462643

#endif

#define ONE_SIDED 1
#define TWO_SIDED 2

#define UNIF    400
#define GAUSS   401
#define TPROC   402
#endif  /* define I_TUBE_H */
