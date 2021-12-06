/* dbesy1.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__3 = 3;
static integer const c__20 = 20;
static integer const c__2 = 2;
static integer const c__1 = 1;

/* Initialized data */

static double const by1cs[20] = { .0320804710061190862932352018628015,
                                  1.26270789743350044953431725999727,
                                  .00649996189992317500097490637314144,
                                  -.0893616452886050411653144160009712,
                                  .0132508812217570954512375510370043,
                                  -8.97905911964835237753039508298105e-4,
                                  3.64736148795830678242287368165349e-5,
                                  -1.00137438166600055549075523845295e-6,
                                  1.99453965739017397031159372421243e-8,
                                  -3.02306560180338167284799332520743e-10,
                                  3.60987815694781196116252914242474e-12,
                                  -3.48748829728758242414552947409066e-14,
                                  2.78387897155917665813507698517333e-16,
                                  -1.86787096861948768766825352533333e-18,
                                  1.06853153391168259757070336e-20,
                                  -5.27472195668448228943872e-23,
                                  2.27019940315566414370133333333333e-25,
                                  -8.59539035394523108693333333333333e-28,
                                  2.88540437983379456e-30,
                                  -8.64754113893717333333333333333333e-33 };
static double const twodpi = .636619772367581343075535053490057;

static float const r__1 = (float) d1mach_(3) * (float).1;
static integer const nty1 = initds_(by1cs, &c__20, &r__1);
static double const xmin = exp(max(+log(d1mach_(1)), -log(d1mach_(2))) + .01) * 1.571;
static double const xsml = sqrt(d1mach_(3) * 4.);

double dbesy1_(double const *x)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double y;
    double ampl;
    double theta;

/* ***BEGIN PROLOGUE  DBESY1 */
/* ***PURPOSE  Compute the Bessel function of the second kind of order */
/*            one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      DOUBLE PRECISION (BESY1-S, DBESY1-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ONE, SECOND KIND, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESY1(X) calculates the double precision Bessel function of the */
/* second kind of order for double precision argument X. */

/* Series for BY1        on the interval  0.          to  1.60000E+01 */
/*                                        with weighted error   8.65E-33 */
/*                                         log weighted error  32.06 */
/*                               significant figures required  32.17 */
/*                                    decimal places required  32.71 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9B1MP, DBESJ1, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESY1 */
/* ***FIRST EXECUTABLE STATEMENT  DBESY1 */
    if (*x <= 0.) {
	xermsg_("SLATEC", "DBESY1", "X IS ZERO OR NEGATIVE", &c__1, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 4.) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "DBESY1", "X SO SMALL Y1 OVERFLOWS", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)23);
    }
    y = 0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    d__1 = y * .125 - 1.;
    ret_val = twodpi * log(*x * .5) * dbesj1_(x) + (dcsevl_(&d__1, by1cs, &
	    nty1) + .5) / *x;
    return ret_val;

L20:
    d9b1mp_(x, &ampl, &theta);
    ret_val = ampl * sin(theta);
    return ret_val;

} /* dbesy1_ */
