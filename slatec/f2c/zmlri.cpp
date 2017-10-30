/* zmlri.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;

int zmlri_(double *zr, double *zi, double const *fnu,
	integer const *kode, integer const *n, double *yr, double *yi, integer *
	nz, double *tol)
{
    /* Initialized data */

    static double const zeror = 0.;
    static double const zeroi = 0.;
    static double const coner = 1.;
    static double const conei = 0.;

    /* System generated locals */
    integer i__1, i__2;
    double d__1, d__2, d__3;

    /* Local variables */
    integer i__, k, m;
    double ak, bk, ap, at;
    integer kk, km;
    double az, p1i, p2i, p1r, p2r, ack, cki, fnf, fkk, ckr;
    integer iaz;
    double rho;
    integer inu;
    double pti, raz, sti, rzi, ptr, str, tst, rzr, rho2, flam, fkap, scle,
	     tfnf;
    integer idum;
    integer ifnu;
    double sumi, sumr;
    integer itime;
    double cnormi, cnormr;

/* ***BEGIN PROLOGUE  ZMLRI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CMLRI-A, ZMLRI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE */
/*     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES. */

/* ***SEE ALSO  ZBESI, ZBESK */
/* ***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZMLRI */
/*     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZMLRI */
    scle = d1mach_(1) / *tol;
    *nz = 0;
    az = zabs_(zr, zi);
    iaz = (integer) az;
    ifnu = (integer) (*fnu);
    inu = ifnu + *n - 1;
    at = iaz + 1.;
    raz = 1. / az;
    str = *zr * raz;
    sti = -(*zi) * raz;
    ckr = str * at * raz;
    cki = sti * at * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    p1r = zeror;
    p1i = zeroi;
    p2r = coner;
    p2i = conei;
    ack = (at + 1.) * raz;
    rho = ack + sqrt(ack * ack - 1.);
    rho2 = rho * rho;
    tst = (rho2 + rho2) / ((rho2 - 1.) * (rho - 1.));
    tst /= *tol;
/* ----------------------------------------------------------------------- */
/*     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES */
/* ----------------------------------------------------------------------- */
    ak = at;
    for (i__ = 1; i__ <= 80; ++i__) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r - (ckr * ptr - cki * pti);
	p2i = p1i - (cki * ptr + ckr * pti);
	p1r = ptr;
	p1i = pti;
	ckr += rzr;
	cki += rzi;
	ap = zabs_(&p2r, &p2i);
	if (ap > tst * ak * ak) {
	    goto L20;
	}
	ak += 1.;
/* L10: */
    }
    goto L110;
L20:
    ++i__;
    k = 0;
    if (inu < iaz) {
	goto L40;
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS */
/* ----------------------------------------------------------------------- */
    p1r = zeror;
    p1i = zeroi;
    p2r = coner;
    p2i = conei;
    at = inu + 1.;
    str = *zr * raz;
    sti = -(*zi) * raz;
    ckr = str * at * raz;
    cki = sti * at * raz;
    ack = at * raz;
    tst = sqrt(ack / *tol);
    itime = 1;
    for (k = 1; k <= 80; ++k) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r - (ckr * ptr - cki * pti);
	p2i = p1i - (ckr * pti + cki * ptr);
	p1r = ptr;
	p1i = pti;
	ckr += rzr;
	cki += rzi;
	ap = zabs_(&p2r, &p2i);
	if (ap < tst) {
	    goto L30;
	}
	if (itime == 2) {
	    goto L40;
	}
	ack = zabs_(&ckr, &cki);
	flam = ack + sqrt(ack * ack - 1.);
	fkap = ap / zabs_(&p1r, &p1i);
	rho = min(flam,fkap);
	tst *= sqrt(rho / (rho * rho - 1.));
	itime = 2;
L30:
	;
    }
    goto L110;
L40:
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION */
/* ----------------------------------------------------------------------- */
    ++k;
/* Computing MAX */
    i__1 = i__ + iaz, i__2 = k + inu;
    kk = max(i__1,i__2);
    fkk = (double) kk;
    p1r = zeror;
    p1i = zeroi;
/* ----------------------------------------------------------------------- */
/*     SCALE P2 AND SUM BY SCLE */
/* ----------------------------------------------------------------------- */
    p2r = scle;
    p2i = zeroi;
    fnf = *fnu - ifnu;
    tfnf = fnf + fnf;
    d__1 = fkk + tfnf + 1.;
    d__2 = fkk + 1.;
    d__3 = tfnf + 1.;
    bk = dgamln_(&d__1, &idum) - dgamln_(&d__2, &idum) - dgamln_(&d__3, &idum)
	    ;
    bk = exp(bk);
    sumr = zeror;
    sumi = zeroi;
    km = kk - inu;
    i__1 = km;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
	p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
	p1r = ptr;
	p1i = pti;
	ak = 1. - tfnf / (fkk + tfnf);
	ack = bk * ak;
	sumr += (ack + bk) * p1r;
	sumi += (ack + bk) * p1i;
	bk = ack;
	fkk += -1.;
/* L50: */
    }
    yr[*n] = p2r;
    yi[*n] = p2i;
    if (*n == 1) {
	goto L70;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
	p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
	p1r = ptr;
	p1i = pti;
	ak = 1. - tfnf / (fkk + tfnf);
	ack = bk * ak;
	sumr += (ack + bk) * p1r;
	sumi += (ack + bk) * p1i;
	bk = ack;
	fkk += -1.;
	m = *n - i__ + 1;
	yr[m] = p2r;
	yi[m] = p2i;
/* L60: */
    }
L70:
    if (ifnu <= 0) {
	goto L90;
    }
    i__1 = ifnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
	p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
	p1r = ptr;
	p1i = pti;
	ak = 1. - tfnf / (fkk + tfnf);
	ack = bk * ak;
	sumr += (ack + bk) * p1r;
	sumi += (ack + bk) * p1i;
	bk = ack;
	fkk += -1.;
/* L80: */
    }
L90:
    ptr = *zr;
    pti = *zi;
    if (*kode == 2) {
	ptr = zeror;
    }
    zlog_(&rzr, &rzi, &str, &sti, &idum);
    p1r = -fnf * str + ptr;
    p1i = -fnf * sti + pti;
    d__1 = fnf + 1.;
    ap = dgamln_(&d__1, &idum);
    ptr = p1r - ap;
    pti = p1i;
/* ----------------------------------------------------------------------- */
/*     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW */
/*     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES */
/* ----------------------------------------------------------------------- */
    p2r += sumr;
    p2i += sumi;
    ap = zabs_(&p2r, &p2i);
    p1r = 1. / ap;
    zexp_(&ptr, &pti, &str, &sti);
    ckr = str * p1r;
    cki = sti * p1r;
    ptr = p2r * p1r;
    pti = -p2i * p1r;
    zmlt_(&ckr, &cki, &ptr, &pti, &cnormr, &cnormi);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	str = yr[i__] * cnormr - yi[i__] * cnormi;
	yi[i__] = yr[i__] * cnormi + yi[i__] * cnormr;
	yr[i__] = str;
/* L100: */
    }
    return 0;
L110:
    *nz = -2;
    return 0;
} /* zmlri_ */
