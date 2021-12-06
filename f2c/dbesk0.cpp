/* dbesk0.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__3 = 3;
static integer const c__16 = 16;
static integer const c__1 = 1;
static integer const c__2 = 2;

/* Initialized data */

static double const bk0cs[16] = { -.0353273932339027687201140060063153,
                                  .344289899924628486886344927529213,
                                  .0359799365153615016265721303687231,
                                  .00126461541144692592338479508673447,
                                  2.28621210311945178608269830297585e-5,
                                  2.53479107902614945730790013428354e-7,
                                  1.90451637722020885897214059381366e-9,
                                  1.03496952576336245851008317853089e-11,
                                  4.25981614279108257652445327170133e-14,
                                  1.3744654358807508969423832544e-16,
                                  3.57089652850837359099688597333333e-19,
                                  7.63164366011643737667498666666666e-22,
                                  1.36542498844078185908053333333333e-24,
                                  2.07527526690666808319999999999999e-27,
                                  2.7128142180729856e-30,
                                  3.08259388791466666666666666666666e-33 };

static float const r__1 = (float) d1mach_(3) * (float).1;
static integer const ntk0 = initds_(bk0cs, &c__16, &r__1);
static double const	xsml = sqrt(d1mach_(3) * 4.);
static double const xmaxt = -log(d1mach_(1));
static double const	xmax = xmaxt - xmaxt * .5 * log(xmaxt) / (xmaxt + .5);

double dbesk0_(double const *x)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double y;

/* ***BEGIN PROLOGUE  DBESK0 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      DOUBLE PRECISION (BESK0-S, DBESK0-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESK0(X) calculates the double precision modified (hyperbolic) */
/* Bessel function of the third kind of order zero for double */
/* precision argument X.  The argument must be greater than zero */
/* but not so large that the result underflows. */

/* Series for BK0        on the interval  0.          to  4.00000E+00 */
/*                                        with weighted error   3.08E-33 */
/*                                         log weighted error  32.51 */
/*                               significant figures required  32.05 */
/*                                    decimal places required  33.11 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESK0 */
/* ***FIRST EXECUTABLE STATEMENT  DBESK0 */
    if (*x <= 0.) {
        xermsg_("SLATEC", "DBESK0", "X IS ZERO OR NEGATIVE", &c__2, &c__2, (
                    ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 2.) {
        goto L20;
    }

    y = 0.;
    if (*x > xsml) {
        y = *x * *x;
    }
    d__1 = y * .5 - 1.;
    ret_val = -log(*x * .5) * dbesi0_(x) - .25 + dcsevl_(&d__1, bk0cs, &ntk0);
    return ret_val;

L20:
    ret_val = 0.;
    if (*x > xmax) {
        xermsg_("SLATEC", "DBESK0", "X SO BIG K0 UNDERFLOWS", &c__1, &c__1, (
                    ftnlen)6, (ftnlen)6, (ftnlen)22);
    }
    if (*x > xmax) {
        return ret_val;
    }

    ret_val = exp(-(*x)) * dbsk0e_(x);

    return ret_val;
} /* dbesk0_ */
