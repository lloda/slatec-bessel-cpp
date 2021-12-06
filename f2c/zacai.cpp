/* zacai.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;

int zacai_(double *zr, double *zi, double const *fnu,
	integer const *kode, integer *mr, integer const *n, double *yr, double *
	yi, integer *nz, double *rl, double *tol, double *elim,
	double *alim)
{
    /* Initialized data */

    static double const pi = 3.14159265358979324;

    /* Local variables */
    double az;
    integer nn, nw;
    double yy, c1i, c2i, c1r, c2r, arg;
    integer iuf;
    double cyi[2], fmr, sgn;
    integer inu;
    double cyr[2], zni, znr, dfnu;
    double ascle, csgni, csgnr, cspni, cspnr;

/* ***BEGIN PROLOGUE  ZACAI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZAIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CACAI-A, ZACAI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA */

/*         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN) */
/*                 MP=PI*MR*CMPLX(0.0,1.0) */

/*     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT */
/*     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1. */
/*     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND */
/*     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON */
/*     IS CALLED FROM ZAIRY. */

/* ***SEE ALSO  ZAIRY */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZACAI */
/*     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZACAI */
    *nz = 0;
    znr = -(*zr);
    zni = -(*zi);
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
/*     POWER SERIES FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    zseri_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol, elim, alim);
    goto L40;
L20:
    if (az < *rl) {
	goto L30;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    zasyi_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, rl, tol, elim,
	    alim);
    if (nw < 0) {
	goto L80;
    }
    goto L40;
L30:
/* ----------------------------------------------------------------------- */
/*     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION */
/* ----------------------------------------------------------------------- */
    zmlri_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol);
    if (nw < 0) {
	goto L80;
    }
L40:
/* ----------------------------------------------------------------------- */
/*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION */
/* ----------------------------------------------------------------------- */
    zbknu_(&znr, &zni, fnu, kode, &c__1, cyr, cyi, &nw, tol, elim, alim);
    if (nw != 0) {
	goto L80;
    }
    fmr = (double) (*mr);
    sgn = -f2c::d_sign(&pi, &fmr);
    csgnr = 0.;
    csgni = sgn;
    if (*kode == 1) {
	goto L50;
    }
    yy = -zni;
    csgnr = -csgni * sin(yy);
    csgni *= cos(yy);
L50:
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = (integer) (*fnu);
    arg = (*fnu - inu) * sgn;
    cspnr = cos(arg);
    cspni = sin(arg);
    if (inu % 2 == 0) {
	goto L60;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L60:
    c1r = cyr[0];
    c1i = cyi[0];
    c2r = yr[1];
    c2i = yi[1];
    if (*kode == 1) {
	goto L70;
    }
    iuf = 0;
    ascle = d1mach_(1) * 1e3 / *tol;
    zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
    *nz += nw;
L70:
    yr[1] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
    yi[1] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
    return 0;
L80:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zacai_ */
