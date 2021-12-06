/* zsqrt.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

int zsqrt_(double *ar, double *ai, double *br,
	double *bi)
{
    /* Initialized data */

    static double const drt = .7071067811865475244008443621;
    static double const dpi = 3.141592653589793238462643383;

    /* Local variables */
    double zm;
    double dtheta;

/* ***BEGIN PROLOGUE  ZSQRT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZSQRT-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A) */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  ZABS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZSQRT */
/* ***FIRST EXECUTABLE STATEMENT  ZSQRT */
    zm = zabs_(ar, ai);
    zm = sqrt(zm);
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
    if (*ai > 0.) {
	goto L60;
    }
    if (*ai < 0.) {
	goto L70;
    }
    *br = 0.;
    *bi = 0.;
    return 0;
L20:
    if (*ar > 0.) {
	goto L30;
    }
    *br = 0.;
    *bi = sqrt((abs(*ar)));
    return 0;
L30:
    *br = sqrt(*ar);
    *bi = 0.;
    return 0;
L40:
    if (*ar < 0.) {
	dtheta += dpi;
    }
L50:
    dtheta *= .5;
    *br = zm * cos(dtheta);
    *bi = zm * sin(dtheta);
    return 0;
L60:
    *br = zm * drt;
    *bi = zm * drt;
    return 0;
L70:
    *br = zm * drt;
    *bi = -zm * drt;
    return 0;
} /* zsqrt_ */
