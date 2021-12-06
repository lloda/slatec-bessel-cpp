/* dgamlm.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;
static integer const c__2 = 2;

int dgamlm_(double *xmin, double *xmax)
{
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    integer i__;
    double xln, xold;
    double alnbig, alnsml;

/* ***BEGIN PROLOGUE  DGAMLM */
/* ***PURPOSE  Compute the minimum and maximum bounds for the argument in */
/*            the Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A, R2 */
/* ***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Calculate the minimum and maximum legal bounds for X in gamma(X). */
/* XMIN and XMAX are not the only bounds, but they are the only non- */
/* trivial ones to calculate. */

/*             Output Arguments -- */
/* XMIN   double precision minimum legal value of X in gamma(X).  Any */
/*        smaller value of X might result in underflow. */
/* XMAX   double precision maximum legal value of X in gamma(X).  Any */
/*        larger value of X might cause overflow. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DGAMLM */
/* ***FIRST EXECUTABLE STATEMENT  DGAMLM */
    alnsml = log(d1mach_(1));
    *xmin = -alnsml;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmin;
	xln = log(*xmin);
	*xmin -= *xmin * ((*xmin + .5) * xln - *xmin - .2258 + alnsml) / (*
		xmin * xln + .5);
	if ((d__1 = *xmin - xold, abs(d__1)) < .005) {
	    goto L20;
	}
/* L10: */
    }
    xermsg_("SLATEC", "DGAMLM", "UNABLE TO FIND XMIN", &c__1, &c__2, (ftnlen)
	    6, (ftnlen)6, (ftnlen)19);

L20:
    *xmin = -(*xmin) + .01;

    alnbig = log(d1mach_(2));
    *xmax = alnbig;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmax;
	xln = log(*xmax);
	*xmax -= *xmax * ((*xmax - .5) * xln - *xmax + .9189 - alnbig) / (*
		xmax * xln - .5);
	if ((d__1 = *xmax - xold, abs(d__1)) < .005) {
	    goto L40;
	}
/* L30: */
    }
    xermsg_("SLATEC", "DGAMLM", "UNABLE TO FIND XMAX", &c__2, &c__2, (ftnlen)
	    6, (ftnlen)6, (ftnlen)19);

L40:
    *xmax += -.01;
/* Computing MAX */
    d__1 = *xmin, d__2 = -(*xmax) + 1.;
    *xmin = max(d__1,d__2);

    return 0;
} /* dgamlm_ */
