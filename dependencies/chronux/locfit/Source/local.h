/*
 *  Most of the changes formerly needed here are handled through
 *  the Makefiles and #ifdef's.
 */

#ifndef I_LF_H
#define I_LF_H

/*
 *   DIRSEP: '/' for unix; '\\' for DOS
 */
#ifdef DOS
#define DIRSEP '\\'
#else
#define DIRSEP '/'
#endif

/*
 * Some older math libraries have no lgamma() function, and gamma(arg)
 * actually returns log(gamma(arg)). If so, you need to change
 * LGAMMA macro below.
 */
#define LGAMMA(arg) lgamma(arg)

/******** NOTHING BELOW HERE NEEDS CHANGING **********/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*RVERSION*/

#ifdef SWINVERSION
#define SVERSION
#include "newredef.h"
#endif

#ifdef RVERSION

/* #typedef int Sint is defined in R.h */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#define list_elt(ev,i) VECTOR_PTR(ev)[i]
#define dval2(ev,i,j) NUMERIC_POINTER(list_elt(ev,i))[j]
#define dvec2(ev,i)   NUMERIC_POINTER(list_elt(ev,i))
#define ivec2(ev,i)   INTEGER_POINTER(list_elt(ev,i))
#undef pmatch

#else

#ifdef SVERSION
#include <S.h>
typedef long int Sint;
typedef s_object * SEXP;
#define list_elt(ev,i) LIST_POINTER(ev)[i]
#define dval2(ev,i,j) NUMERIC_POINTER(list_elt(ev,i))[j]
#define dvec2(ev,i)   NUMERIC_POINTER(list_elt(ev,i))
#define ivec2(ev,i)   INTEGER_POINTER(list_elt(ev,i))
#define ALLOW_MODULES
#else
typedef int Sint;
#endif

#endif

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define printf mexPrintf
#endif

#ifdef RVERSION
#undef LGAMMA
#define LGAMMA(arg) Rf_lgammafn(arg)
extern double Rf_lgammafn();
#define SVERSION
#endif

#include "mutil.h"
#include "tube.h"

#include "lfcons.h"

typedef char varname[15];

#ifdef CVERSION
#include "cversion.h"
#endif

#include "lfstruc.h"
#include "design.h"
#include "lffuns.h"

#ifdef CVERSION
#undef printf
#define printf lfprintf
extern int lfprintf(const char *format, ...);
extern int printe(const char *format, ...);
#else
#define printe printf
#endif

#ifdef ERROR
#undef ERROR
#endif

#ifdef WARN
#undef WARN
#endif

#define ERROR(args) {printe("Error: "); printe args; printe("\n"); lf_error=1;}
#define WARN(args)  {printe("Warning: "); printe args; printe("\n"); }

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SGN(x) (((x)>0) ? 1 : -1)
#define SQR(x) ((x)*(x))
#define NOSLN 0.1278433
#define GFACT 2.5
#define EFACT 3.0

#define MAXCOLOR 20
#define MAXWIN 5

#define ISWAP(a,b) { int zz; zz = a; a = b; b = zz; }
extern int lf_error;

#endif /* I_LF_H */
