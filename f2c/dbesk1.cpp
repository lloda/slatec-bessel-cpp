/* dbesk1.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__3 = 3;
static integer const c__16 = 16;
static integer const c__2 = 2;
static integer const c__1 = 1;

/* Initialized data */

static double const bk1cs[16] = { .025300227338947770532531120868533,
                                  -.35315596077654487566723831691801,
                                  -.12261118082265714823479067930042,
                                  -.0069757238596398643501812920296083,
                                  -1.7302889575130520630176507368979e-4,
                                  -2.4334061415659682349600735030164e-6,
                                  -2.2133876307347258558315252545126e-8,
                                  -1.4114883926335277610958330212608e-10,
                                  -6.6669016941993290060853751264373e-13,
                                  -2.4274498505193659339263196864853e-15,
                                  -7.023863479386287597178379712e-18,
                                  -1.6543275155100994675491029333333e-20,
                                  -3.2338347459944491991893333333333e-23,
                                  -5.3312750529265274999466666666666e-26,
                                  -7.5130407162157226666666666666666e-29,
                                  -9.1550857176541866666666666666666e-32 };

static float const r__1 = (float) d1mach_(3) * (float).1;
static integer const ntk1 = initds_(bk1cs, &c__16, &r__1);
/* Computing MAX */
static double const xmin = exp(max(+log(d1mach_(1)), -log(d1mach_(2))) + .01);
static double const xsml = sqrt(d1mach_(3) * 4.);
static double const xmaxt = -log(d1mach_(1));
static double const xmax = xmaxt - xmaxt * .5 * log(xmaxt) / (xmaxt + .5);

double dbesk1_(double const *x)
{

    /* System generated locals */
    double d__1;

    /* Local variables */
    double y;

/* ***BEGIN PROLOGUE  DBESK1 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      DOUBLE PRECISION (BESK1-S, DBESK1-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESK1(X) calculates the double precision modified (hyperbolic) */
/* Bessel function of the third kind of order one for double precision */
/* argument X.  The argument must be large enough that the result does */
/* not overflow and small enough that the result does not underflow. */

/* Series for BK1        on the interval  0.          to  4.00000E+00 */
/*                                        with weighted error   9.16E-32 */
/*                                         log weighted error  31.04 */
/*                               significant figures required  30.61 */
/*                                    decimal places required  31.64 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESK1 */
/* ***FIRST EXECUTABLE STATEMENT  DBESK1 */
    if (*x <= 0.) {
        xermsg_("SLATEC", "DBESK1", "X IS ZERO OR NEGATIVE", &c__2, &c__2,
                (ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 2.) {
        if (*x > xmax) {
            xermsg_("SLATEC", "DBESK1", "X SO BIG K1 UNDERFLOWS", &c__1, &c__1,
                    (ftnlen)6, (ftnlen)6, (ftnlen)22);
            return 0.;
        } else {
            return exp(-(*x)) * dbsk1e_(x);
        }
    } else {
        if (*x < xmin) {
            xermsg_("SLATEC", "DBESK1", "X SO SMALL K1 OVERFLOWS", &c__3, &c__2,
                    (ftnlen)6, (ftnlen)6, (ftnlen)23);
        }
        y = 0.;
        if (*x > xsml) {
            y = *x * *x;
        }
        d__1 = y * .5 - 1.;
        return log(*x * .5) * dbesi1_(x) + (dcsevl_(&d__1, bk1cs, &ntk1) + .75) / *x;
    }
} /* dbesk1_ */
