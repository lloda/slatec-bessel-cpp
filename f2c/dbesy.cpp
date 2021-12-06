/* dbesy.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__15 = 15;
static integer const c__5 = 5;
static integer const c__1 = 1;
static integer const c__2 = 2;
static integer const c__6 = 6;

int dbesy_(double const *x, double const *fnu, integer const *n, double *y)
{
    /* Initialized data */

    static integer nulim[2] = { 70,100 };

    /* System generated locals */
    integer i__1;

    /* Local variables. Some initialized to avoid -Wmaybe-uninitialized */
    integer i__, j;
    double s, w[2], s1, s2 = 0.;
    integer nb;
    double cn;
    integer nd;
    double fn;
    integer nn;
    double tm = 0., wk[7], w2n, ran;
    integer nud;
    double dnu, azn, trx = 0., xxn, elim;
    integer iflw;
    double xlim, flgjy;

/* ***BEGIN PROLOGUE  DBESY */
/* ***PURPOSE  Implement forward recursion on the three term recursion */
/*            relation for a sequence of non-negative order Bessel */
/*            functions Y/SUB(FNU+I-1)/(X), I=1,...,N for real, positive */
/*            X and non-negative orders FNU. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10A3 */
/* ***TYPE      DOUBLE PRECISION (BESY-S, DBESY-D) */
/* ***KEYWORDS  SPECIAL FUNCTIONS, Y BESSEL FUNCTION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract  **** a double precision routine **** */
/*         DBESY implements forward recursion on the three term */
/*         recursion relation for a sequence of non-negative order Bessel */
/*         functions Y/sub(FNU+I-1)/(X), I=1,N for real X .GT. 0.0D0 and */
/*         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and */
/*         FNU+1 are obtained from DBSYNU which computes by a power */
/*         series for X .LE. 2, the K Bessel function of an imaginary */
/*         argument for 2 .LT. X .LE. 20 and the asymptotic expansion for */
/*         X .GT. 20. */

/*         If FNU .GE. NULIM, the uniform asymptotic expansion is coded */
/*         in DASYJY for orders FNU and FNU+1 to start the recursion. */
/*         NULIM is 70 or 100 depending on whether N=1 or N .GE. 2.  An */
/*         overflow test is made on the leading term of the asymptotic */
/*         expansion before any extensive computation is done. */

/*         The maximum number of significant digits obtainable */
/*         is the smaller of 14 and the number of digits carried in */
/*         double precision arithmetic. */

/*     Description of Arguments */

/*         Input */
/*           X      - X .GT. 0.0D0 */
/*           FNU    - order of the initial Y function, FNU .GE. 0.0D0 */
/*           N      - number of members in the sequence, N .GE. 1 */

/*         Output */
/*           Y      - a vector whose first N components contain values */
/*                    for the sequence Y(I)=Y/sub(FNU+I-1)/(X), I=1,N. */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */

/* ***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate */
/*                 or Large Orders, NPL Mathematical Tables 6, Her */
/*                 Majesty's Stationery Office, London, 1962. */
/*               N. M. Temme, On the numerical evaluation of the modified */
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/*               N. M. Temme, On the numerical evaluation of the ordinary */
/*                 Bessel function of the second kind, Journal of */
/*                 Computational Physics 21, (1976), pp. 343-350. */
/* ***ROUTINES CALLED  D1MACH, DASYJY, DBESY0, DBESY1, DBSYNU, DYAIRY, */
/*                    I1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBESY */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DBESY */
    nn = -i1mach_(15);
    elim = (nn * d1mach_(5) - 3.) * 2.303;
    xlim = d1mach_(1) * 1e3;
    if (*fnu < 0.) {
	goto L140;
    }
    if (*x <= 0.) {
	goto L150;
    }
    if (*x < xlim) {
	goto L170;
    }
    if (*n < 1) {
	goto L160;
    }

/*     ND IS A DUMMY VARIABLE FOR N */

    nd = *n;
    nud = (integer) (*fnu);
    dnu = *fnu - nud;
    nn = min(2,nd);
    fn = *fnu + *n - 1;
    if (fn < 2.) {
	goto L100;
    }

/*     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION) */
/*     FOR THE LAST ORDER, FNU+N-1.GE.NULIM */

    xxn = *x / fn;
    w2n = 1. - xxn * xxn;
    if (w2n <= 0.) {
	goto L10;
    }
    ran = sqrt(w2n);
    azn = log((ran + 1.) / xxn) - ran;
    cn = fn * azn;
    if (cn > elim) {
	goto L170;
    }
L10:
    if (nud < nulim[nn - 1]) {
	goto L20;
    }

/*     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM */

    flgjy = -1.;
    dasyjy_(dyairy_, x, fnu, &flgjy, &nn, &y[1], wk, &iflw);
    if (iflw != 0) {
	goto L170;
    }
    if (nn == 1) {
	return 0;
    }
    trx = 2. / *x;
    tm = (*fnu + *fnu + 2.) / *x;
    goto L80;

L20:
    if (dnu != 0.) {
	goto L30;
    }
    s1 = dbesy0_(x);
    if (nud == 0 && nd == 1) {
	goto L70;
    }
    s2 = dbesy1_(x);
    goto L40;
L30:
    nb = 2;
    if (nud == 0 && nd == 1) {
	nb = 1;
    }
    dbsynu_(x, &dnu, &nb, w);
    s1 = w[0];
    if (nb == 1) {
	goto L70;
    }
    s2 = w[1];
L40:
    trx = 2. / *x;
    tm = (dnu + dnu + 2.) / *x;
/*     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2) */
    if (nd == 1) {
	--nud;
    }
    if (nud > 0) {
	goto L50;
    }
    if (nd > 1) {
	goto L70;
    }
    s1 = s2;
    goto L70;
L50:
    i__1 = nud;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = s2;
	s2 = tm * s2 - s1;
	s1 = s;
	tm += trx;
/* L60: */
    }
    if (nd == 1) {
	s1 = s2;
    }
L70:
    y[1] = s1;
    if (nd == 1) {
	return 0;
    }
    y[2] = s2;
L80:
    if (nd == 2) {
	return 0;
    }
/*     FORWARD RECUR FROM FNU+2 TO FNU+N-1 */
    i__1 = nd;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = tm * y[i__ - 1] - y[i__ - 2];
	tm += trx;
/* L90: */
    }
    return 0;

L100:
/*     OVERFLOW TEST */
    if (fn <= 1.) {
	goto L110;
    }
    if (-fn * (log(*x) - .693) > elim) {
	goto L170;
    }
L110:
    if (dnu == 0.) {
	goto L120;
    }
    dbsynu_(x, fnu, &nd, &y[1]);
    return 0;
L120:
    j = nud;
    if (j == 1) {
	goto L130;
    }
    ++j;
    y[j] = dbesy0_(x);
    if (nd == 1) {
	return 0;
    }
    ++j;
L130:
    y[j] = dbesy1_(x);
    if (nd == 1) {
	return 0;
    }
    trx = 2. / *x;
    tm = trx;
    goto L80;



L140:
    xermsg_("SLATEC", "DBESY", "ORDER, FNU, LESS THAN ZERO", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)26);
    return 0;
L150:
    xermsg_("SLATEC", "DBESY", "X LESS THAN OR EQUAL TO ZERO", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)28);
    return 0;
L160:
    xermsg_("SLATEC", "DBESY", "N LESS THAN ONE", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)5, (ftnlen)15);
    return 0;
L170:
    xermsg_("SLATEC", "DBESY", "OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL",
	    &c__6, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)43);
    return 0;
} /* dbesy_ */
