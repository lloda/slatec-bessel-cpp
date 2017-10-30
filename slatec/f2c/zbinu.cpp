/* zbinu.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;
static integer const c__2 = 2;

int zbinu_(double *zr, double *zi, double const *fnu,
	integer const *kode, integer const *n, double *cyr, double *cyi, integer *
	nz, double *rl, double *fnul, double *tol, double *
	elim, double *alim)
{
    /* Initialized data */

    static double const zeror = 0.;
    static double const zeroi = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    double az;
    integer nn, nw;
    double cwi[2], cwr[2];
    integer nui, inw;
    double dfnu;
    integer nlast;

/* ***BEGIN PROLOGUE  ZBINU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CBINU-A, ZBINU-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY */
/* ***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZBINU */
    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZBINU */
    *nz = 0;
    az = zabs_(zr, zi);
    nn = *n;
    dfnu = *fnu + (*n - 1);
    if (az <= 2.) {
	goto L10;
    }
    if (az * az * .25 > dfnu + 1.) {
	goto L20;
    }
L10:
/* ----------------------------------------------------------------------- */
/*     POWER SERIES */
/* ----------------------------------------------------------------------- */
    zseri_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol, elim, alim);
    inw = abs(nw);
    *nz += inw;
    nn -= inw;
    if (nn == 0) {
	return 0;
    }
    if (nw >= 0) {
	goto L120;
    }
    dfnu = *fnu + (nn - 1);
L20:
    if (az < *rl) {
	goto L40;
    }
    if (dfnu <= 1.) {
	goto L30;
    }
    if (az + az < dfnu * dfnu) {
	goto L50;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z */
/* ----------------------------------------------------------------------- */
L30:
    zasyi_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, rl, tol, elim, alim)
	    ;
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L40:
    if (dfnu <= 1.) {
	goto L70;
    }
L50:
/* ----------------------------------------------------------------------- */
/*     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM */
/* ----------------------------------------------------------------------- */
    zuoik_(zr, zi, fnu, kode, &c__1, &nn, &cyr[1], &cyi[1], &nw, tol, elim,
	    alim);
    if (nw < 0) {
	goto L130;
    }
    *nz += nw;
    nn -= nw;
    if (nn == 0) {
	return 0;
    }
    dfnu = *fnu + (nn - 1);
    if (dfnu > *fnul) {
	goto L110;
    }
    if (az > *fnul) {
	goto L110;
    }
L60:
    if (az > *rl) {
	goto L80;
    }
L70:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE SERIES */
/* ----------------------------------------------------------------------- */
    zmlri_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L80:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN */
/* ----------------------------------------------------------------------- */
    zuoik_(zr, zi, fnu, kode, &c__2, &c__2, cwr, cwi, &nw, tol, elim, alim);
    if (nw >= 0) {
	goto L100;
    }
    *nz = nn;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cyr[i__] = zeror;
	cyi[i__] = zeroi;
/* L90: */
    }
    return 0;
L100:
    if (nw > 0) {
	goto L130;
    }
    zwrsk_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, cwr, cwi, tol, elim,
	     alim);
    if (nw < 0) {
	goto L130;
    }
    goto L120;
L110:
/* ----------------------------------------------------------------------- */
/*     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD */
/* ----------------------------------------------------------------------- */
    nui = (integer) (*fnul - dfnu + 1);
    nui = max(nui,0);
    zbuni_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, &nui, &nlast, fnul,
	    tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    *nz += nw;
    if (nlast == 0) {
	goto L120;
    }
    nn = nlast;
    goto L60;
L120:
    return 0;
L130:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zbinu_ */
