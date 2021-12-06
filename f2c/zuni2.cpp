/* zuni2.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;
static integer const c__0 = 0;
static integer const c__2 = 2;

int zuni2_(double *zr, double *zi, double const *fnu,
	integer const *kode, integer const *n, double *yr, double *yi, integer *
	nz, integer *nlast, double *fnul, double *tol, double *
	elim, double *alim)
{
    /* Initialized data */

    static double const zeror = 0.;
    static double const zeroi = 0.;
    static double const coner = 1.;
    static double const cipr[4] = { 1.,0.,-1.,0. };
    static double const cipi[4] = { 0.,1.,0.,-1. };
    static double const hpi = 1.57079632679489662;
    static double const aic = 1.265512123484645396;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, nd;
    double fn;
    integer in, nn, nw;
    double c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, aii, ang, car;
    integer nai;
    double air, zbi, cyi[2], sar;
    integer nuf, inu;
    double bry[3], raz, sti, zbr, zni, cyr[2], rzi, str, znr, rzr, daii,
	    cidi, aarg;
    integer ndai;
    double dair, aphi, argi, cscl, phii, crsc, argr;
    integer idum;
    double phir, csrr[3], cssr[3], rast;
    integer iflag;
    double ascle, asumi, bsumi;
    double asumr, bsumr;
    double zeta1i, zeta2i, zeta1r, zeta2r;

/* ***BEGIN PROLOGUE  ZUNI2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CUNI2-A, ZUNI2-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF */
/*     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I */
/*     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO. */

/*     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC */
/*     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET. */
/*     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER */
/*     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL. */
/*     Y(I)=CZERO FOR I=NLAST+1,N */

/* ***SEE ALSO  ZBESI, ZBESK */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZUCHK, ZUNHJ, ZUOIK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZUNI2 */
/*     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS, */
/*    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZUNI2 */
    *nz = 0;
    nd = *n;
    *nlast = 0;
/* ----------------------------------------------------------------------- */
/*     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG- */
/*     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE, */
/*     EXP(ALIM)=EXP(ELIM)*TOL */
/* ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    crsc = *tol;
    cssr[0] = cscl;
    cssr[1] = coner;
    cssr[2] = crsc;
    csrr[0] = crsc;
    csrr[1] = coner;
    csrr[2] = cscl;
    bry[0] = d1mach_(1) * 1e3 / *tol;
/* ----------------------------------------------------------------------- */
/*     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI */
/* ----------------------------------------------------------------------- */
    znr = *zi;
    zni = -(*zr);
    zbr = *zr;
    zbi = *zi;
    cidi = -coner;
    inu = (integer) (*fnu);
    ang = hpi * (*fnu - inu);
    c2r = cos(ang);
    c2i = sin(ang);
    car = c2r;
    sar = c2i;
    in = inu + *n - 1;
    in = in % 4 + 1;
    str = c2r * cipr[in - 1] - c2i * cipi[in - 1];
    c2i = c2r * cipi[in - 1] + c2i * cipr[in - 1];
    c2r = str;
    if (*zi > 0.) {
	goto L10;
    }
    znr = -znr;
    zbi = -zbi;
    cidi = -cidi;
    c2i = -c2i;
L10:
/* ----------------------------------------------------------------------- */
/*     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER */
/* ----------------------------------------------------------------------- */
    fn = max(*fnu,1.);
    zunhj_(&znr, &zni, &fn, &c__1, tol, &phir, &phii, &argr, &argi, &zeta1r, &
	    zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &bsumi);
    if (*kode == 1) {
	goto L20;
    }
    str = zbr + zeta2r;
    sti = zbi + zeta2i;
    rast = fn / zabs_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti;
    goto L30;
L20:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
L30:
    rs1 = s1r;
    if (abs(rs1) > *elim) {
	goto L150;
    }
L40:
    nn = min(2,nd);
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fn = *fnu + (nd - i__);
	zunhj_(&znr, &zni, &fn, &c__0, tol, &phir, &phii, &argr, &argi, &
		zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &
		bsumi);
	if (*kode == 1) {
	    goto L50;
	}
	str = zbr + zeta2r;
	sti = zbi + zeta2i;
	rast = fn / zabs_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zeta1r + str;
	s1i = -zeta1i + sti + abs(*zi);
	goto L60;
L50:
	s1r = -zeta1r + zeta2r;
	s1i = -zeta1i + zeta2i;
L60:
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (abs(rs1) > *elim) {
	    goto L120;
	}
	if (i__ == 1) {
	    iflag = 2;
	}
	if (abs(rs1) < *alim) {
	    goto L70;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
	aphi = zabs_(&phir, &phii);
	aarg = zabs_(&argr, &argi);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (abs(rs1) > *elim) {
	    goto L120;
	}
	if (i__ == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L70;
	}
	if (i__ == 1) {
	    iflag = 3;
	}
L70:
/* ----------------------------------------------------------------------- */
/*     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR */
/*     EXPONENT EXTREMES */
/* ----------------------------------------------------------------------- */
	zairy_(&argr, &argi, &c__0, &c__2, &air, &aii, &nai, &idum);
	zairy_(&argr, &argi, &c__1, &c__2, &dair, &daii, &ndai, &idum);
	str = dair * bsumr - daii * bsumi;
	sti = dair * bsumi + daii * bsumr;
	str += air * asumr - aii * asumi;
	sti += air * asumi + aii * asumr;
	s2r = phir * str - phii * sti;
	s2i = phir * sti + phii * str;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L80;
	}
	zuchk_(&s2r, &s2i, &nw, bry, tol);
	if (nw != 0) {
	    goto L120;
	}
L80:
	if (*zi <= 0.) {
	    s2i = -s2i;
	}
	str = s2r * c2r - s2i * c2i;
	s2i = s2r * c2i + s2i * c2r;
	s2r = str;
	cyr[i__ - 1] = s2r;
	cyi[i__ - 1] = s2i;
	j = nd - i__ + 1;
	yr[j] = s2r * csrr[iflag - 1];
	yi[j] = s2i * csrr[iflag - 1];
	str = -c2i * cidi;
	c2i = c2r * cidi;
	c2r = str;
/* L90: */
    }
    if (nd <= 2) {
	goto L110;
    }
    raz = 1. / zabs_(zr, zi);
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    bry[1] = 1. / bry[0];
    bry[2] = d1mach_(2);
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    k = nd - 2;
    fn = (double) k;
    i__1 = nd;
    for (i__ = 3; i__ <= i__1; ++i__) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (*fnu + fn) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (*fnu + fn) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[k] = c2r;
	yi[k] = c2i;
	--k;
	fn += -1.;
	if (iflag >= 3) {
	    goto L100;
	}
	str = abs(c2r);
	sti = abs(c2i);
	c2m = max(str,sti);
	if (c2m <= ascle) {
	    goto L100;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	c1r = csrr[iflag - 1];
L100:
	;
    }
L110:
    return 0;
L120:
    if (rs1 > 0.) {
	goto L140;
    }
/* ----------------------------------------------------------------------- */
/*     SET UNDERFLOW AND UPDATE PARAMETERS */
/* ----------------------------------------------------------------------- */
    yr[nd] = zeror;
    yi[nd] = zeroi;
    ++(*nz);
    --nd;
    if (nd == 0) {
	goto L110;
    }
    zuoik_(zr, zi, fnu, kode, &c__1, &nd, &yr[1], &yi[1], &nuf, tol, elim,
	    alim);
    if (nuf < 0) {
	goto L140;
    }
    nd -= nuf;
    *nz += nuf;
    if (nd == 0) {
	goto L110;
    }
    fn = *fnu + (nd - 1);
    if (fn < *fnul) {
	goto L130;
    }
/*      FN = CIDI */
/*      J = NUF + 1 */
/*      K = MOD(J,4) + 1 */
/*      S1R = CIPR(K) */
/*      S1I = CIPI(K) */
/*      IF (FN.LT.0.0D0) S1I = -S1I */
/*      STR = C2R*S1R - C2I*S1I */
/*      C2I = C2R*S1I + C2I*S1R */
/*      C2R = STR */
    in = inu + nd - 1;
    in = in % 4 + 1;
    c2r = car * cipr[in - 1] - sar * cipi[in - 1];
    c2i = car * cipi[in - 1] + sar * cipr[in - 1];
    if (*zi <= 0.) {
	c2i = -c2i;
    }
    goto L40;
L130:
    *nlast = nd;
    return 0;
L140:
    *nz = -1;
    return 0;
L150:
    if (rs1 > 0.) {
	goto L140;
    }
    *nz = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yr[i__] = zeror;
	yi[i__] = zeroi;
/* L160: */
    }
    return 0;
} /* zuni2_ */
