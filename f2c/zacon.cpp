/* zacon.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;
static integer const c__2 = 2;

int zacon_(double *zr, double *zi, double const *fnu,
                            integer const *kode, integer *mr, integer const *n, double *yr, double *
                            yi, integer *nz, double *rl, double *fnul, double *tol,
                            double *elim, double *alim)
{
    /* Initialized data */

    static double const pi = 3.14159265358979324;
    static double const zeror = 0.;
    static double const coner = 1.;

    /* System generated locals */
    integer i__1;

    /* Local variables. Some initialized to avoid -Wmaybe-uninitialized */
    integer i__;
    double fn;
    integer nn, nw;
    double yy, c1i, c2i, c1m, as2, c1r, c2r, s1i, s2i, s1r, s2r, cki, arg,
	     ckr, cpn;
    integer iuf;
    double cyi[2], fmr, csr, azn, sgn;
    integer inu;
    double bry[3], cyr[2], pti, spn, sti, zni, rzi, ptr, str, znr, rzr,
	    sc1i, sc2i = 0., sc1r, sc2r = 0., cscl, cscr;
    double csrr[3], cssr[3], razn;
    integer kflag;
    double ascle, bscle, csgni, csgnr, cspni, cspnr;

/* ***BEGIN PROLOGUE  ZACON */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CACON-A, ZACON-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA */

/*         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN) */
/*                 MP=PI*MR*CMPLX(0.0,1.0) */

/*     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT */
/*     HALF Z PLANE */

/* ***SEE ALSO  ZBESH, ZBESK */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZACON */
/*     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST, */
/*    *S1,S2,Y,Z,ZN */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZACON */
    *nz = 0;
    znr = -(*zr);
    zni = -(*zi);
    nn = *n;
    zbinu_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, rl, fnul, tol,
	    elim, alim);
    if (nw < 0) {
	goto L90;
    }
/* ----------------------------------------------------------------------- */
/*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION */
/* ----------------------------------------------------------------------- */
    nn = min(2,*n);
    zbknu_(&znr, &zni, fnu, kode, &nn, cyr, cyi, &nw, tol, elim, alim);
    if (nw != 0) {
	goto L90;
    }
    s1r = cyr[0];
    s1i = cyi[0];
    fmr = (double) (*mr);
    sgn = -f2c::d_sign(&pi, &fmr);
    csgnr = zeror;
    csgni = sgn;
    if (*kode == 1) {
	goto L10;
    }
    yy = -zni;
    cpn = cos(yy);
    spn = sin(yy);
    zmlt_(&csgnr, &csgni, &cpn, &spn, &csgnr, &csgni);
L10:
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = (integer) (*fnu);
    arg = (*fnu - inu) * sgn;
    cpn = cos(arg);
    spn = sin(arg);
    cspnr = cpn;
    cspni = spn;
    if (inu % 2 == 0) {
	goto L20;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L20:
    iuf = 0;
    c1r = s1r;
    c1i = s1i;
    c2r = yr[1];
    c2i = yi[1];
    ascle = d1mach_(1) * 1e3 / *tol;
    if (*kode == 1) {
	goto L30;
    }
    zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
    *nz += nw;
    sc1r = c1r;
    sc1i = c1i;
L30:
    zmlt_(&cspnr, &cspni, &c1r, &c1i, &str, &sti);
    zmlt_(&csgnr, &csgni, &c2r, &c2i, &ptr, &pti);
    yr[1] = str + ptr;
    yi[1] = sti + pti;
    if (*n == 1) {
	return 0;
    }
    cspnr = -cspnr;
    cspni = -cspni;
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = s2r;
    c1i = s2i;
    c2r = yr[2];
    c2i = yi[2];
    if (*kode == 1) {
	goto L40;
    }
    zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
    *nz += nw;
    sc2r = c1r;
    sc2i = c1i;
L40:
    zmlt_(&cspnr, &cspni, &c1r, &c1i, &str, &sti);
    zmlt_(&csgnr, &csgni, &c2r, &c2i, &ptr, &pti);
    yr[2] = str + ptr;
    yi[2] = sti + pti;
    if (*n == 2) {
	return 0;
    }
    cspnr = -cspnr;
    cspni = -cspni;
    azn = zabs_(&znr, &zni);
    razn = 1. / azn;
    str = znr * razn;
    sti = -zni * razn;
    rzr = (str + str) * razn;
    rzi = (sti + sti) * razn;
    fn = *fnu + 1.;
    ckr = fn * rzr;
    cki = fn * rzi;
/* ----------------------------------------------------------------------- */
/*     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS */
/* ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    cscr = *tol;
    cssr[0] = cscl;
    cssr[1] = coner;
    cssr[2] = cscr;
    csrr[0] = cscr;
    csrr[1] = coner;
    csrr[2] = cscl;
    bry[0] = ascle;
    bry[1] = 1. / ascle;
    bry[2] = d1mach_(2);
    as2 = zabs_(&s2r, &s2i);
    kflag = 2;
    if (as2 > bry[0]) {
	goto L50;
    }
    kflag = 1;
    goto L60;
L50:
    if (as2 < bry[1]) {
	goto L60;
    }
    kflag = 3;
L60:
    bscle = bry[kflag - 1];
    s1r *= cssr[kflag - 1];
    s1i *= cssr[kflag - 1];
    s2r *= cssr[kflag - 1];
    s2i *= cssr[kflag - 1];
    csr = csrr[kflag - 1];
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	str = s2r;
	sti = s2i;
	s2r = ckr * str - cki * sti + s1r;
	s2i = ckr * sti + cki * str + s1i;
	s1r = str;
	s1i = sti;
	c1r = s2r * csr;
	c1i = s2i * csr;
	str = c1r;
	sti = c1i;
	c2r = yr[i__];
	c2i = yi[i__];
	if (*kode == 1) {
	    goto L70;
	}
	if (iuf < 0) {
	    goto L70;
	}
	zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
	*nz += nw;
	sc1r = sc2r;
	sc1i = sc2i;
	sc2r = c1r;
	sc2i = c1i;
	if (iuf != 3) {
	    goto L70;
	}
	iuf = -4;
	s1r = sc1r * cssr[kflag - 1];
	s1i = sc1i * cssr[kflag - 1];
	s2r = sc2r * cssr[kflag - 1];
	s2i = sc2i * cssr[kflag - 1];
	str = sc2r;
	sti = sc2i;
L70:
	ptr = cspnr * c1r - cspni * c1i;
	pti = cspnr * c1i + cspni * c1r;
	yr[i__] = ptr + csgnr * c2r - csgni * c2i;
	yi[i__] = pti + csgnr * c2i + csgni * c2r;
	ckr += rzr;
	cki += rzi;
	cspnr = -cspnr;
	cspni = -cspni;
	if (kflag >= 3) {
	    goto L80;
	}
	ptr = abs(c1r);
	pti = abs(c1i);
	c1m = max(ptr,pti);
	if (c1m <= bscle) {
	    goto L80;
	}
	++kflag;
	bscle = bry[kflag - 1];
	s1r *= csr;
	s1i *= csr;
	s2r = str;
	s2i = sti;
	s1r *= cssr[kflag - 1];
	s1i *= cssr[kflag - 1];
	s2r *= cssr[kflag - 1];
	s2i *= cssr[kflag - 1];
	csr = csrr[kflag - 1];
L80:
	;
    }
    return 0;
L90:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zacon_ */
