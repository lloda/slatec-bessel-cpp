/* dbesi0.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__3 = 3;
static integer const c__18 = 18;
static integer const c__2 = 2;

/* Initialized data */

static double const bi0cs[18] = { -.07660547252839144951081894976243285,
                                  1.927337953993808269952408750881196,
                                  .2282644586920301338937029292330415,
                                  .01304891466707290428079334210691888,
                                  4.344270900816487451378682681026107e-4,
                                  9.422657686001934663923171744118766e-6,
                                  1.434006289510691079962091878179957e-7,
                                  1.613849069661749069915419719994611e-9,
                                  1.396650044535669699495092708142522e-11,
                                  9.579451725505445344627523171893333e-14,
                                  5.333981859862502131015107744e-16,
                                  2.458716088437470774696785919999999e-18,
                                  9.535680890248770026944341333333333e-21,
                                  3.154382039721427336789333333333333e-23,
                                  9.004564101094637431466666666666666e-26,
                                  2.240647369123670016e-28,
                                  4.903034603242837333333333333333333e-31,
                                  9.508172606122666666666666666666666e-34 };

static float const r__1 = (float) d1mach_(3) * (float).1;
static integer const nti0 = initds_(bi0cs, &c__18, &r__1);
static double const xsml = sqrt(d1mach_(3) * 4.5);
static double const xmax = log(d1mach_(2));

double dbesi0_(double const *x)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double y;

/* ***BEGIN PROLOGUE  DBESI0 */
/* ***PURPOSE  Compute the hyperbolic Bessel function of the first kind */
/*            of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      DOUBLE PRECISION (BESI0-S, DBESI0-D) */
/* ***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESI0(X) calculates the double precision modified (hyperbolic) */
/* Bessel function of the first kind of order zero and double */
/* precision argument X. */

/* Series for BI0        on the interval  0.          to  9.00000E+00 */
/*                                        with weighted error   9.51E-34 */
/*                                         log weighted error  33.02 */
/*                               significant figures required  33.31 */
/*                                    decimal places required  33.65 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DBSI0E, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESI0 */
/* ***FIRST EXECUTABLE STATEMENT  DBESI0 */
    y = abs(*x);
    if (y > 3.) {
        goto L20;
    }

    ret_val = 1.;
    if (y > xsml) {
        d__1 = y * y / 4.5 - 1.;
        ret_val = dcsevl_(&d__1, bi0cs, &nti0) + 2.75;
    }
    return ret_val;

L20:
    if (y > xmax) {
        xermsg_("SLATEC", "DBESI0", "ABS(X) SO BIG I0 OVERFLOWS", &c__2, &
                c__2, (ftnlen)6, (ftnlen)6, (ftnlen)26);
    }

    ret_val = exp(y) * dbsi0e_(x);

    return ret_val;
} /* dbesi0_ */
