/* dlngam.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__2 = 2;
static integer const c__4 = 4;
static integer const c__3 = 3;
static integer const c__1 = 1;

static double const sq2pil = .91893853320467274178032973640562;
static double const sqpi2l = .225791352644727432363097614947441;
static double const pi = 3.1415926535897932384626433832795;

static double const temp = 1. / log(d1mach_(2));
static double const xmax = temp * d1mach_(2);
static double const dxrel = sqrt(d1mach_(4));

double dlngam_(double const *x)
{
    /* Initialized data */

    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    double y;
    double sinpiy;

/* ***BEGIN PROLOGUE  DLNGAM */
/* ***PURPOSE  Compute the logarithm of the absolute value of the Gamma */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      DOUBLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C) */
/* ***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DLNGAM(X) calculates the double precision logarithm of the */
/* absolute value of the Gamma function for double precision */
/* argument X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9LGMC, DGAMMA, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  DLNGAM */
/* ***FIRST EXECUTABLE STATEMENT  DLNGAM */
    y = abs(*x);
    if (y > 10.) {
        goto L20;
    }

/* LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0 */

    return log(abs(dgamma_(x)));

/* LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0 */

L20:
    if (y > xmax) {
        xermsg_("SLATEC", "DLNGAM", "ABS(X) SO BIG DLNGAM OVERFLOWS", &c__2, &
                c__2, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    }

    if (*x > 0.) {
        return sq2pil + (*x - .5) * log(*x) - *x + d9lgmc_(&y);
    }

    sinpiy = (d__1 = sin(pi * y), abs(d__1));
    if (sinpiy == 0.) {
        xermsg_("SLATEC", "DLNGAM", "X IS A NEGATIVE INTEGER", &c__3, &c__2, (
                    ftnlen)6, (ftnlen)6, (ftnlen)23);
    }

    d__2 = *x - .5;
    if ((d__1 = (*x - f2c::d_int(&d__2)) / *x, abs(d__1)) < dxrel) {
        xermsg_("SLATEC", "DLNGAM", "ANSWER LT HALF PRECISION BECAUSE X TOO \
NEAR NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)60);
    }

    return sqpi2l + (*x - .5) * log(y) - *x - log(sinpiy) - d9lgmc_(&y);

} /* dlngam_ */
