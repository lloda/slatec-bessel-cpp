/* zdiv.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

int zdiv_(double const *ar, double const *ai, double const *br,
          double const *bi, double *cr, double *ci)
{
    double ca, cb, cc, cd, bm;

/* ***BEGIN PROLOGUE  ZDIV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZDIV-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     DOUBLE PRECISION COMPLEX DIVIDE C=A/B. */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  ZABS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZDIV */
/* ***FIRST EXECUTABLE STATEMENT  ZDIV */
    bm = 1. / zabs_(br, bi);
    cc = *br * bm;
    cd = *bi * bm;
    ca = (*ar * cc + *ai * cd) * bm;
    cb = (*ai * cc - *ar * cd) * bm;
    *cr = ca;
    *ci = cb;
    return 0;
} /* zdiv_ */
