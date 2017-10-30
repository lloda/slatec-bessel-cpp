/* zlog.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

int zlog_(double *ar, double *ai, double *br,
	double *bi, integer *ierr)
{
    /* Initialized data */

    static double const dpi = 3.141592653589793238462643383;
    static double const dhpi = 1.570796326794896619231321696;

    /* Local variables */
    double zm;
    double dtheta;

/* ***BEGIN PROLOGUE  ZLOG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZLOG-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A) */
/*     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0) */
/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  ZABS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZLOG */
/* ***FIRST EXECUTABLE STATEMENT  ZLOG */
    *ierr = 0;
    if (*ar == 0.) {
	goto L10;
    }
    if (*ai == 0.) {
	goto L20;
    }
    dtheta = atan(*ai / *ar);
    if (dtheta <= 0.) {
	goto L40;
    }
    if (*ar < 0.) {
	dtheta -= dpi;
    }
    goto L50;
L10:
    if (*ai == 0.) {
	goto L60;
    }
    *bi = dhpi;
    *br = log((abs(*ai)));
    if (*ai < 0.) {
	*bi = -(*bi);
    }
    return 0;
L20:
    if (*ar > 0.) {
	goto L30;
    }
    *br = log((abs(*ar)));
    *bi = dpi;
    return 0;
L30:
    *br = log(*ar);
    *bi = 0.;
    return 0;
L40:
    if (*ar < 0.) {
	dtheta += dpi;
    }
L50:
    zm = zabs_(ar, ai);
    *br = log(zm);
    *bi = dtheta;
    return 0;
L60:
    *ierr = 1;
    return 0;
} /* zlog_ */
