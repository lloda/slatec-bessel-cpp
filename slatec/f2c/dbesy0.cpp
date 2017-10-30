/* dbesy0.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__19 = 19;
static integer const c__1 = 1;
static integer const c__2 = 2;

/* Initialized data */

static double const by0cs[19] = { -.01127783939286557321793980546028,
                                  -.1283452375604203460480884531838,
                                  -.1043788479979424936581762276618,
                                  .02366274918396969540924159264613,
                                  -.002090391647700486239196223950342,
                                  1.039754539390572520999246576381e-4,
                                  -3.369747162423972096718775345037e-6,
                                  7.729384267670667158521367216371e-8,
                                  -1.324976772664259591443476068964e-9,
                                  1.764823261540452792100389363158e-11,
                                  -1.881055071580196200602823012069e-13,
                                  1.641865485366149502792237185749e-15,
                                  -1.19565943860460608574599100672e-17,
                                  7.377296297440185842494112426666e-20,
                                  -3.906843476710437330740906666666e-22,
                                  1.79550366443615794982912e-24,
                                  -7.229627125448010478933333333333e-27,
                                  2.571727931635168597333333333333e-29,
                                  -8.141268814163694933333333333333e-32 };
static double const twodpi = .636619772367581343075535053490057;

static float const r__1 = (float) d1mach_(3) * (float).1;
static integer const nty0 = initds_(by0cs, &c__19, &r__1);
static double const xsml = sqrt(d1mach_(3) * 4.);

double dbesy0_(double const *x)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double y;
    double ampl;
    double theta;

/* ***BEGIN PROLOGUE  DBESY0 */
/* ***PURPOSE  Compute the Bessel function of the second kind of order */
/*            zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      DOUBLE PRECISION (BESY0-S, DBESY0-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESY0(X) calculates the double precision Bessel function of the */
/* second kind of order zero for double precision argument X. */

/* Series for BY0        on the interval  0.          to  1.60000E+01 */
/*                                        with weighted error   8.14E-32 */
/*                                         log weighted error  31.09 */
/*                               significant figures required  30.31 */
/*                                    decimal places required  31.73 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9B0MP, DBESJ0, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESY0 */
/* ***FIRST EXECUTABLE STATEMENT  DBESY0 */
    if (*x <= 0.) {
	xermsg_("SLATEC", "DBESY0", "X IS ZERO OR NEGATIVE", &c__1, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 4.) {
	goto L20;
    }

    y = 0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    d__1 = y * .125 - 1.;
    ret_val = twodpi * log(*x * .5) * dbesj0_(x) + .375 + dcsevl_(&d__1,
	    by0cs, &nty0);
    return ret_val;

L20:
    d9b0mp_(x, &ampl, &theta);
    ret_val = ampl * sin(theta);
    return ret_val;

} /* dbesy0_ */
