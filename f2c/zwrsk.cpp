/* zwrsk.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__2 = 2;
static integer const c__1 = 1;

int zwrsk_(double *zrr, double *zri, double const *fnu,
                            integer const *kode, integer const *n, double *yr, double *yi, integer * nz,
                            double *cwr, double *cwi, double *tol, double *
                            elim, double *alim)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, nw;
    double c1i, c2i, c1r, c2r, act, acw, cti, ctr, pti, sti, ptr, str, ract;
    double ascle, csclr, cinui, cinur;

/* ***BEGIN PROLOGUE  ZWRSK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CWRSK-A, ZWRSK-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY */
/*     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN */

/* ***SEE ALSO  ZBESI, ZBESK */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZBKNU, ZRATI */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZWRSK */
/*     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR */
/* ***FIRST EXECUTABLE STATEMENT  ZWRSK */
/* ----------------------------------------------------------------------- */
/*     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS */
/*     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE */
/*     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --yi;
    --yr;
    --cwr;
    --cwi;

    /* Function Body */
    *nz = 0;
    zbknu_(zrr, zri, fnu, kode, &c__2, &cwr[1], &cwi[1], &nw, tol, elim, alim)
	    ;
    if (nw != 0) {
	goto L50;
    }
    zrati_(zrr, zri, fnu, n, &yr[1], &yi[1], tol);
/* ----------------------------------------------------------------------- */
/*     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z), */
/*     R(FNU+J-1,Z)=Y(J),  J=1,...,N */
/* ----------------------------------------------------------------------- */
    cinur = 1.;
    cinui = 0.;
    if (*kode == 1) {
	goto L10;
    }
    cinur = cos(*zri);
    cinui = sin(*zri);
L10:
/* ----------------------------------------------------------------------- */
/*     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH */
/*     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE */
/*     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT */
/*     THE RESULT IS ON SCALE. */
/* ----------------------------------------------------------------------- */
    acw = zabs_(&cwr[2], &cwi[2]);
    ascle = d1mach_(1) * 1e3 / *tol;
    csclr = 1.;
    if (acw > ascle) {
	goto L20;
    }
    csclr = 1. / *tol;
    goto L30;
L20:
    ascle = 1. / ascle;
    if (acw < ascle) {
	goto L30;
    }
    csclr = *tol;
L30:
    c1r = cwr[1] * csclr;
    c1i = cwi[1] * csclr;
    c2r = cwr[2] * csclr;
    c2i = cwi[2] * csclr;
    str = yr[1];
    sti = yi[1];
/* ----------------------------------------------------------------------- */
/*     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS */
/*     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT) */
/* ----------------------------------------------------------------------- */
    ptr = str * c1r - sti * c1i;
    pti = str * c1i + sti * c1r;
    ptr += c2r;
    pti += c2i;
    ctr = *zrr * ptr - *zri * pti;
    cti = *zrr * pti + *zri * ptr;
    act = zabs_(&ctr, &cti);
    ract = 1. / act;
    ctr *= ract;
    cti = -cti * ract;
    ptr = cinur * ract;
    pti = cinui * ract;
    cinur = ptr * ctr - pti * cti;
    cinui = ptr * cti + pti * ctr;
    yr[1] = cinur * csclr;
    yi[1] = cinui * csclr;
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ptr = str * cinur - sti * cinui;
	cinui = str * cinui + sti * cinur;
	cinur = ptr;
	str = yr[i__];
	sti = yi[i__];
	yr[i__] = cinur * csclr;
	yi[i__] = cinui * csclr;
/* L40: */
    }
    return 0;
L50:
    *nz = -1;
    if (nw == -2) {
	*nz = -2;
    }
    return 0;
} /* zwrsk_ */
