/* d9lgmc.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__15 = 15;
static integer const c__3 = 3;
static integer const c__1 = 1;
static integer const c__2 = 2;

/* Initialized data */

static double const algmcs[15] = { .1666389480451863247205729650822,
                                       -1.384948176067563840732986059135e-5,
                                       9.810825646924729426157171547487e-9,
                                       -1.809129475572494194263306266719e-11,
                                       6.221098041892605227126015543416e-14,
                                       -3.399615005417721944303330599666e-16,
                                       2.683181998482698748957538846666e-18,
                                       -2.868042435334643284144622399999e-20,
                                       3.962837061046434803679306666666e-22,
                                       -6.831888753985766870111999999999e-24,
                                       1.429227355942498147573333333333e-25,
                                       -3.547598158101070547199999999999e-27,
                                       1.025680058010470912e-28,
                                       -3.401102254316748799999999999999e-30,
                                       1.276642195630062933333333333333e-31 };

static float const r__1 = (float) d1mach_(3);
static integer const nalgm = initds_(algmcs, &c__15, &r__1);
static double const xbig = 1. / sqrt(d1mach_(3));
static double const xmax = exp((min(log(d1mach_(2) / 12.),
                                        -log(d1mach_(1) * 12.))));

double d9lgmc_(double const *x)
{
    /* System generated locals */
    double ret_val, d__1, d__2;

/* ***BEGIN PROLOGUE  D9LGMC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the log Gamma correction factor so that */
/*            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X */
/*            + D9LGMC(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB, */
/*             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the log gamma correction factor for X .GE. 10. so that */
/* LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X) */

/* Series for ALGM       on the interval  0.          to  1.00000E-02 */
/*                                        with weighted error   1.28E-31 */
/*                                         log weighted error  30.89 */
/*                               significant figures required  29.81 */
/*                                    decimal places required  31.48 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  D9LGMC */
/* ***FIRST EXECUTABLE STATEMENT  D9LGMC */
    if (*x < 10.) {
	xermsg_("SLATEC", "D9LGMC", "X MUST BE GE 10", &c__1, &c__2, (ftnlen)
		6, (ftnlen)6, (ftnlen)15);
    }
    if (*x >= xmax) {
	goto L20;
    }

    ret_val = 1. / (*x * 12.);
    if (*x < xbig) {
/* Computing 2nd power */
	d__2 = 10. / *x;
	d__1 = d__2 * d__2 * 2. - 1.;
	ret_val = dcsevl_(&d__1, algmcs, &nalgm) / *x;
    }
    return ret_val;

L20:
    ret_val = 0.;
    xermsg_("SLATEC", "D9LGMC", "X SO BIG D9LGMC UNDERFLOWS", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)26);
    return ret_val;

} /* d9lgmc_ */
