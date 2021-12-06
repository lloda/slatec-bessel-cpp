/* zrati.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

int zrati_(double *zr, double *zi, double const *fnu,
	integer const *n, double *cyr, double *cyi, double *tol)
{
    /* Initialized data */

    static double const czeror = 0.;
    static double const czeroi = 0.;
    static double const coner = 1.;
    static double const conei = 0.;
    static double const rt2 = 1.41421356237309505;

    /* System generated locals */
    integer i__1;
    double d__1;

    /* Local variables */
    integer i__, k;
    double ak;
    integer id, kk;
    double az, ap1, ap2, p1i, p2i, t1i, p1r, p2r, t1r, arg, rak, rho;
    integer inu;
    double pti, tti, rzi, ptr, ttr, rzr, rap1, flam, dfnu, fdnu;
    integer magz;
    integer idnu;
    double fnup;
    double test, test1, amagz;
    integer itime;
    double cdfnui, cdfnur;

/* ***BEGIN PROLOGUE  ZRATI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CRATI-A, ZRATI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD */
/*     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD */
/*     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B, */
/*     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973, */
/*     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER, */
/*     BY D. J. SOOKNE. */

/* ***SEE ALSO  ZBESH, ZBESI, ZBESK */
/* ***ROUTINES CALLED  ZABS, ZDIV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZRATI */
    /* Parameter adjustments */
    --cyi;
    --cyr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZRATI */
    az = zabs_(zr, zi);
    inu = (integer) (*fnu);
    idnu = inu + *n - 1;
    magz = (integer) az;
    amagz = (double) (magz + 1);
    fdnu = (double) idnu;
    fnup = max(amagz,fdnu);
    id = idnu - magz - 1;
    itime = 1;
    k = 1;
    ptr = 1. / az;
    rzr = ptr * (*zr + *zr) * ptr;
    rzi = -ptr * (*zi + *zi) * ptr;
    t1r = rzr * fnup;
    t1i = rzi * fnup;
    p2r = -t1r;
    p2i = -t1i;
    p1r = coner;
    p1i = conei;
    t1r += rzr;
    t1i += rzi;
    if (id > 0) {
	id = 0;
    }
    ap2 = zabs_(&p2r, &p2i);
    ap1 = zabs_(&p1r, &p1i);
/* ----------------------------------------------------------------------- */
/*     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU */
/*     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT */
/*     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR */
/*     PREMATURELY. */
/* ----------------------------------------------------------------------- */
    arg = (ap2 + ap2) / (ap1 * *tol);
    test1 = sqrt(arg);
    test = test1;
    rap1 = 1. / ap1;
    p1r *= rap1;
    p1i *= rap1;
    p2r *= rap1;
    p2i *= rap1;
    ap2 *= rap1;
L10:
    ++k;
    ap1 = ap2;
    ptr = p2r;
    pti = p2i;
    p2r = p1r - (t1r * ptr - t1i * pti);
    p2i = p1i - (t1r * pti + t1i * ptr);
    p1r = ptr;
    p1i = pti;
    t1r += rzr;
    t1i += rzi;
    ap2 = zabs_(&p2r, &p2i);
    if (ap1 <= test) {
	goto L10;
    }
    if (itime == 2) {
	goto L20;
    }
    ak = zabs_(&t1r, &t1i) * .5;
    flam = ak + sqrt(ak * ak - 1.);
/* Computing MIN */
    d__1 = ap2 / ap1;
    rho = min(d__1,flam);
    test = test1 * sqrt(rho / (rho * rho - 1.));
    itime = 2;
    goto L10;
L20:
    kk = k + 1 - id;
    ak = (double) kk;
    t1r = ak;
    t1i = czeroi;
    dfnu = *fnu + (*n - 1);
    p1r = 1. / ap2;
    p1i = czeroi;
    p2r = czeror;
    p2i = czeroi;
    i__1 = kk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ptr = p1r;
	pti = p1i;
	rap1 = dfnu + t1r;
	ttr = rzr * rap1;
	tti = rzi * rap1;
	p1r = ptr * ttr - pti * tti + p2r;
	p1i = ptr * tti + pti * ttr + p2i;
	p2r = ptr;
	p2i = pti;
	t1r -= coner;
/* L30: */
    }
    if (p1r != czeror || p1i != czeroi) {
	goto L40;
    }
    p1r = *tol;
    p1i = *tol;
L40:
    zdiv_(&p2r, &p2i, &p1r, &p1i, &cyr[*n], &cyi[*n]);
    if (*n == 1) {
	return 0;
    }
    k = *n - 1;
    ak = (double) k;
    t1r = ak;
    t1i = czeroi;
    cdfnur = *fnu * rzr;
    cdfnui = *fnu * rzi;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[k + 1];
	pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[k + 1];
	ak = zabs_(&ptr, &pti);
	if (ak != czeror) {
	    goto L50;
	}
	ptr = *tol;
	pti = *tol;
	ak = *tol * rt2;
L50:
	rak = coner / ak;
	cyr[k] = rak * ptr * rak;
	cyi[k] = -rak * pti * rak;
	t1r -= coner;
	--k;
/* L60: */
    }
    return 0;
} /* zrati_ */
