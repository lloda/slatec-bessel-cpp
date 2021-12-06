/* zmlt.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

int zmlt_(double *ar, double *ai, double *br,
	double *bi, double *cr, double *ci)
{
    double ca, cb;

/* ***BEGIN PROLOGUE  ZMLT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZMLT-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B. */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZMLT */
/* ***FIRST EXECUTABLE STATEMENT  ZMLT */
    ca = *ar * *br - *ai * *bi;
    cb = *ar * *bi + *ai * *br;
    *cr = ca;
    *ci = cb;
    return 0;
} /* zmlt_ */
