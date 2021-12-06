/* dbesj0.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

integer const c__3 = 3;
integer const c__19 = 19;

    /* Initialized data */

double const bj0cs[19] = { .10025416196893913701073127264074,
                               -.66522300776440513177678757831124,
                               .2489837034982813137046046872668,
                               -.033252723170035769653884341503854,
                               .0023114179304694015462904924117729,
                               -9.9112774199508092339048519336549e-5,
                               2.8916708643998808884733903747078e-6,
                               -6.1210858663032635057818407481516e-8,
                               9.8386507938567841324768748636415e-10,
                               -1.2423551597301765145515897006836e-11,
                               1.2654336302559045797915827210363e-13,
                               -1.0619456495287244546914817512959e-15,
                               7.4706210758024567437098915584e-18,
                               -4.4697032274412780547627007999999e-20,
                               2.3024281584337436200523093333333e-22,
                               -1.0319144794166698148522666666666e-24,
                               4.06081782748733227008e-27,
                               -1.4143836005240913919999999999999e-29,
                               4.391090549669888e-32 };

float const r__1 = (float) d1mach_(3) * (float).1;
integer const ntj0 = initds_(bj0cs, &c__19, &r__1);
double const xsml = sqrt(d1mach_(3) * 8.);

double dbesj0_(double const *x)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double y;
    double ampl;
    double theta;

/* ***BEGIN PROLOGUE  DBESJ0 */
/* ***PURPOSE  Compute the Bessel function of the first kind of order */
/*            zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      DOUBLE PRECISION (BESJ0-S, DBESJ0-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESJ0(X) calculates the double precision Bessel function of */
/* the first kind of order zero for double precision argument X. */

/* Series for BJ0        on the interval  0.          to  1.60000E+01 */
/*                                        with weighted error   4.39E-32 */
/*                                         log weighted error  31.36 */
/*                               significant figures required  31.21 */
/*                                    decimal places required  32.00 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9B0MP, DCSEVL, INITDS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DBESJ0 */
/* ***FIRST EXECUTABLE STATEMENT  DBESJ0 */

    y = abs(*x);
    if (y > 4.) {
	goto L20;
    }

    ret_val = 1.;
    if (y > xsml) {
	d__1 = y * .125 * y - 1.;
	ret_val = dcsevl_(&d__1, bj0cs, &ntj0);
    }
    return ret_val;

L20:
    d9b0mp_(&y, &ampl, &theta);
    ret_val = ampl * cos(theta);

    return ret_val;
} /* dbesj0_ */
