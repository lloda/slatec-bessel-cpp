/* dbesk.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__15 = 15;
static integer const c__5 = 5;
static integer const c__1 = 1;
static integer const c__2 = 2;
static integer const c__6 = 6;

int dbesk_(double const *x, double const *fnu, integer const *kode, integer const *n, double *y, integer *nz)
{
    /* Initialized data */

    static integer nulim[2] = { 35,70 };

    /* System generated locals */
    integer i__1;

    /* Local variables. Some initialized to avoid -Wmaybe-uninitialized */
    integer i__, j, k;
    double s, t, w[2], s1, s2 = 0.;
    integer nb;
    double cn;
    integer nd;
    double fn;
    integer nn;
    double tm = 0.;
    integer mz;
    double zn, gln, fnn;
    integer nud;
    double dnu, gnu, etx, trx = 0., rtz, elim, xlim, flgik;

/* ***BEGIN PROLOGUE  DBESK */
/* ***PURPOSE  Implement forward recursion on the three term recursion */
/*            relation for a sequence of non-negative order Bessel */
/*            functions K/SUB(FNU+I-1)/(X), or scaled Bessel functions */
/*            EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N for real, positive */
/*            X and non-negative orders FNU. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10B3 */
/* ***TYPE      DOUBLE PRECISION (BESK-S, DBESK-D) */
/* ***KEYWORDS  K BESSEL FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract  **** a double precision routine **** */
/*         DBESK implements forward recursion on the three term */
/*         recursion relation for a sequence of non-negative order Bessel */
/*         functions K/sub(FNU+I-1)/(X), or scaled Bessel functions */
/*         EXP(X)*K/sub(FNU+I-1)/(X), I=1,..,N for real X .GT. 0.0D0 and */
/*         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and */
/*         FNU+1 are obtained from DBSKNU to start the recursion.  If */
/*         FNU .GE. NULIM, the uniform asymptotic expansion is used for */
/*         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or */
/*         70 depending on whether N=1 or N .GE. 2.  Under and overflow */
/*         tests are made on the leading term of the asymptotic expansion */
/*         before any extensive computation is done. */

/*         The maximum number of significant digits obtainable */
/*         is the smaller of 14 and the number of digits carried in */
/*         double precision arithmetic. */

/*     Description of Arguments */

/*         Input      X,FNU are double precision */
/*           X      - X .GT. 0.0D0 */
/*           FNU    - order of the initial K function, FNU .GE. 0.0D0 */
/*           KODE   - a parameter to indicate the scaling option */
/*                    KODE=1 returns Y(I)=       K/sub(FNU+I-1)/(X), */
/*                                        I=1,...,N */
/*                    KODE=2 returns Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), */
/*                                        I=1,...,N */
/*           N      - number of members in the sequence, N .GE. 1 */

/*         Output     Y is double precision */
/*           Y      - a vector whose first N components contain values */
/*                    for the sequence */
/*                    Y(I)=       k/sub(FNU+I-1)/(X), I=1,...,N  or */
/*                    Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N */
/*                    depending on KODE */
/*           NZ     - number of components of Y set to zero due to */
/*                    underflow with KODE=1, */
/*                    NZ=0   , normal return, computation completed */
/*                    NZ .NE. 0, first NZ components of Y set to zero */
/*                             due to underflow, Y(I)=0.0D0, I=1,...,NZ */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */
/*         Underflow with KODE=1 -  a non-fatal error (NZ .NE. 0) */

/* ***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate */
/*                 or Large Orders, NPL Mathematical Tables 6, Her */
/*                 Majesty's Stationery Office, London, 1962. */
/*               N. M. Temme, On the numerical evaluation of the modified */
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/* ***ROUTINES CALLED  D1MACH, DASYIK, DBESK0, DBESK1, DBSK0E, DBSK1E, */
/*                    DBSKNU, I1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBESK */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DBESK */
    nn = -i1mach_(15);
    elim = (nn * d1mach_(5) - 3.) * 2.303;
    xlim = d1mach_(1) * 1e3;
    if (*kode < 1 || *kode > 2) {
	goto L280;
    }
    if (*fnu < 0.) {
	goto L290;
    }
    if (*x <= 0.) {
	goto L300;
    }
    if (*x < xlim) {
	goto L320;
    }
    if (*n < 1) {
	goto L310;
    }
    etx = (double) (*kode - 1);

/*     ND IS A DUMMY VARIABLE FOR N */
/*     GNU IS A DUMMY VARIABLE FOR FNU */
/*     NZ = NUMBER OF UNDERFLOWS ON KODE=1 */

    nd = *n;
    *nz = 0;
    nud = (integer) (*fnu);
    dnu = *fnu - nud;
    gnu = *fnu;
    nn = min(2,nd);
    fn = *fnu + *n - 1;
    fnn = fn;
    if (fn < 2.) {
	goto L150;
    }

/*     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION) */
/*     FOR THE LAST ORDER, FNU+N-1.GE.NULIM */

    zn = *x / fn;
    if (zn == 0.) {
	goto L320;
    }
    rtz = sqrt(zn * zn + 1.);
    gln = log((rtz + 1.) / zn);
    t = rtz * (1. - etx) + etx / (zn + rtz);
    cn = -fn * (t - gln);
    if (cn > elim) {
	goto L320;
    }
    if (nud < nulim[nn - 1]) {
	goto L30;
    }
    if (nn == 1) {
	goto L20;
    }
L10:

/*     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION) */
/*     FOR THE FIRST ORDER, FNU.GE.NULIM */

    fn = gnu;
    zn = *x / fn;
    rtz = sqrt(zn * zn + 1.);
    gln = log((rtz + 1.) / zn);
    t = rtz * (1. - etx) + etx / (zn + rtz);
    cn = -fn * (t - gln);
L20:
    if (cn < -elim) {
	goto L230;
    }

/*     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM */

    flgik = -1.;
    dasyik_(x, &gnu, kode, &flgik, &rtz, &cn, &nn, &y[1]);
    if (nn == 1) {
	goto L240;
    }
    trx = 2. / *x;
    tm = (gnu + gnu + 2.) / *x;
    goto L130;

L30:
    if (*kode == 2) {
	goto L40;
    }

/*     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X) */
/*     FOR ORDER DNU */

    if (*x > elim) {
	goto L230;
    }
L40:
    if (dnu != 0.) {
	goto L80;
    }
    if (*kode == 2) {
	goto L50;
    }
    s1 = dbesk0_(x);
    goto L60;
L50:
    s1 = dbsk0e_(x);
L60:
    if (nud == 0 && nd == 1) {
	goto L120;
    }
    if (*kode == 2) {
	goto L70;
    }
    s2 = dbesk1_(x);
    goto L90;
L70:
    s2 = dbsk1e_(x);
    goto L90;
L80:
    nb = 2;
    if (nud == 0 && nd == 1) {
	nb = 1;
    }
    dbsknu_(x, &dnu, kode, &nb, w, nz);
    s1 = w[0];
    if (nb == 1) {
	goto L120;
    }
    s2 = w[1];
L90:
    trx = 2. / *x;
    tm = (dnu + dnu + 2.) / *x;
/*     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2) */
    if (nd == 1) {
	--nud;
    }
    if (nud > 0) {
	goto L100;
    }
    if (nd > 1) {
	goto L120;
    }
    s1 = s2;
    goto L120;
L100:
    i__1 = nud;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = s2;
	s2 = tm * s2 + s1;
	s1 = s;
	tm += trx;
/* L110: */
    }
    if (nd == 1) {
	s1 = s2;
    }
L120:
    y[1] = s1;
    if (nd == 1) {
	goto L240;
    }
    y[2] = s2;
L130:
    if (nd == 2) {
	goto L240;
    }
/*     FORWARD RECUR FROM FNU+2 TO FNU+N-1 */
    i__1 = nd;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = tm * y[i__ - 1] + y[i__ - 2];
	tm += trx;
/* L140: */
    }
    goto L240;

L150:
/*     UNDERFLOW TEST FOR KODE=1 */
    if (*kode == 2) {
	goto L160;
    }
    if (*x > elim) {
	goto L230;
    }
L160:
/*     OVERFLOW TEST */
    if (fn <= 1.) {
	goto L170;
    }
    if (-fn * (log(*x) - .693) > elim) {
	goto L320;
    }
L170:
    if (dnu == 0.) {
	goto L180;
    }
    dbsknu_(x, fnu, kode, &nd, &y[1], &mz);
    goto L240;
L180:
    j = nud;
    if (j == 1) {
	goto L210;
    }
    ++j;
    if (*kode == 2) {
	goto L190;
    }
    y[j] = dbesk0_(x);
    goto L200;
L190:
    y[j] = dbsk0e_(x);
L200:
    if (nd == 1) {
	goto L240;
    }
    ++j;
L210:
    if (*kode == 2) {
	goto L220;
    }
    y[j] = dbesk1_(x);
    goto L240;
L220:
    y[j] = dbsk1e_(x);
    goto L240;

/*     UPDATE PARAMETERS ON UNDERFLOW */

L230:
    ++nud;
    --nd;
    if (nd == 0) {
	goto L240;
    }
    nn = min(2,nd);
    gnu += 1.;
    if (fnn < 2.) {
	goto L230;
    }
    if (nud < nulim[nn - 1]) {
	goto L230;
    }
    goto L10;
L240:
    *nz = *n - nd;
    if (*nz == 0) {
	return 0;
    }
    if (nd == 0) {
	goto L260;
    }
    i__1 = nd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n - i__ + 1;
	k = nd - i__ + 1;
	y[j] = y[k];
/* L250: */
    }
L260:
    i__1 = *nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.;
/* L270: */
    }
    return 0;



L280:
    xermsg_("SLATEC", "DBESK", "SCALING OPTION, KODE, NOT 1 OR 2", &c__2, &
	    c__1, (ftnlen)6, (ftnlen)5, (ftnlen)32);
    return 0;
L290:
    xermsg_("SLATEC", "DBESK", "ORDER, FNU, LESS THAN ZERO", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)26);
    return 0;
L300:
    xermsg_("SLATEC", "DBESK", "X LESS THAN OR EQUAL TO ZERO", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)28);
    return 0;
L310:
    xermsg_("SLATEC", "DBESK", "N LESS THAN ONE", &c__2, &c__1, (ftnlen)6, (
	    ftnlen)5, (ftnlen)15);
    return 0;
L320:
    xermsg_("SLATEC", "DBESK", "OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL",
	    &c__6, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)43);
    return 0;
} /* dbesk_ */
