/* algo-680-erf.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

/* This routine isn't originally from SLATEC, but it makes sense to package it in. */

#include "slatec-internal.hpp"


/*      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM. */
/*      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*      VOL. 16, NO. 1, PP. 47. */
int wofz_(double *xi, double *yi, double *u,
                           double *v, logical *flag__)
{
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    logical a, b;
    double c__, h__;
    integer i__, j, n;
    double x, y, h2 = 0., u1, v1, u2, v2, w1;
    integer nu;
    double rx, ry, sx, sy, tx, ty;
    integer np1, kapn;
    double xabs, yabs, daux, qrho, xaux, xsum, ysum, xabsq, xquad, yquad;

/*  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES */
/*  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z), */
/*  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I */
/*  MEANS SQRT(-1). */
/*  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT */
/*  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT */
/*  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO */
/*  OF THE FUNCTION. */
/*  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION. */


/*  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS : */
/*     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF */
/*                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE */
/*                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION */
/*                FLOATING-POINT ARITHMETIC */
/*     RMAXEXP  = LN(RMAX) - LN(2) */
/*     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION */
/*                GONIOMETRIC FUNCTION (DCOS, DSIN, ...) */
/*  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL */
/*  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS */


/*  PARAMETER LIST */
/*     XI     = REAL      PART OF Z */
/*     YI     = IMAGINARY PART OF Z */
/*     U      = REAL      PART OF W(Z) */
/*     V      = IMAGINARY PART OF W(Z) */
/*     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL */
/*              OCCUR OR NOT; TYPE LOGICAL; */
/*              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING */
/*              MEANING : */
/*              FLAG=.FALSE. : NO ERROR CONDITION */
/*              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE */
/*                             BECOMES INACTIVE */
/*  XI, YI      ARE THE INPUT-PARAMETERS */
/*  U, V, FLAG  ARE THE OUTPUT-PARAMETERS */

/*  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI) */

/*  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE */
/*  PUT TO 0 UPON UNDERFLOW; */

/*  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF */
/*  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE. */





/*      DOUBLE PRECISION D1MACH, FACTOR, RMAX, RMAXREAL, RMAXEXP */

/*      PARAMETER (FACTOR   = 1.12837916709551257389615890312154517D0, */
/*     *           RMAXREAL = 1.340780792994259D+154, */
/*     *           RMAXEXP  = 709.0895657128241D0, */
/*     *           RMAXGONI = 0.6746518850690209D10) */
/*      RMAX = D1MACH(2) */
/*      RMAXREAL = DSQRT(RMAX) */
/*      RMAXEXP = DLOG(RMAX)-DLOG(2D0) */

    *flag__ = FALSE_;

    xabs = abs(*xi);
    yabs = abs(*yi);
    x = xabs / (float)6.3;
    y = yabs / (float)4.4;


/*     THE FOLLOWING IF-STATEMENT PROTECTS */
/*     QRHO = (X**2 + Y**2) AGAINST OVERFLOW */

    if (xabs > 5e153 || yabs > 5e153) {
	goto L100;
    }

/* Computing 2nd power */
    d__1 = x;
/* Computing 2nd power */
    d__2 = y;
    qrho = d__1 * d__1 + d__2 * d__2;

/* Computing 2nd power */
    d__1 = xabs;
    xabsq = d__1 * d__1;
/* Computing 2nd power */
    d__1 = yabs;
    xquad = xabsq - d__1 * d__1;
    yquad = xabs * 2 * yabs;

    a = qrho < .085264;

    if (a) {

/*  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED */
/*  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297) */
/*  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED */
/*  ACCURACY */

	qrho = (1 - y * (float).85) * sqrt(qrho);
	d__1 = qrho * 72 + 6;
	n = f2c::i_dnnt(&d__1);
	j = (n << 1) + 1;
	xsum = (float)1. / j;
	ysum = 0.;
	for (i__ = n; i__ >= 1; --i__) {
	    j += -2;
	    xaux = (xsum * xquad - ysum * yquad) / i__;
	    ysum = (xsum * yquad + ysum * xquad) / i__;
	    xsum = xaux + (float)1. / j;
/* L10: */
	}
	u1 = (xsum * yabs + ysum * xabs) *
		-1.12837916709551257389615890312154517 + (float)1.;
	v1 = (xsum * xabs - ysum * yabs) *
		1.12837916709551257389615890312154517;
	daux = exp(-xquad);
	u2 = daux * cos(yquad);
	v2 = -daux * sin(yquad);

	*u = u1 * u2 - v1 * v2;
	*v = u1 * v2 + v1 * u2;

    } else {

/*  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE */
/*  CONTINUED FRACTION */
/*  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED */
/*  ACCURACY */

/*  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED */
/*  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION */
/*  IS USED TO CALCULATE THE DERIVATIVES OF W(Z) */
/*  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED */
/*  TO OBTAIN THE REQUIRED ACCURACY */
/*  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED */
/*  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY */


	if (qrho > (float)1.) {
	    h__ = 0.;
	    kapn = 0;
	    qrho = sqrt(qrho);
	    nu = (integer) (1442 / (qrho * 26 + 77) + 3);
	} else {
	    qrho = (1 - y) * sqrt(1 - qrho);
	    h__ = qrho * (float)1.88;
	    h2 = h__ * 2;
	    d__1 = qrho * 34 + 7;
	    kapn = f2c::i_dnnt(&d__1);
	    d__1 = qrho * 26 + 16;
	    nu = f2c::i_dnnt(&d__1);
	}

	b = h__ > (float)0.;

    double qlambda = b ? f2c::pow_di(&h2, &kapn) : 0.; /* to avoid -Wmaybe-uninitialized */

	rx = (float)0.;
	ry = (float)0.;
	sx = (float)0.;
	sy = (float)0.;

	for (n = nu; n >= 0; --n) {
	    np1 = n + 1;
	    tx = yabs + h__ + np1 * rx;
	    ty = xabs - np1 * ry;
/* Computing 2nd power */
	    d__1 = tx;
/* Computing 2nd power */
	    d__2 = ty;
	    c__ = (float).5 / (d__1 * d__1 + d__2 * d__2);
	    rx = c__ * tx;
	    ry = c__ * ty;
	    if (b && n <= kapn) {
            tx = qlambda + sx;
            sx = rx * tx - ry * sy;
            sy = ry * tx + rx * sy;
            qlambda /= h2;
	    }
/* L11: */
	}

	if (h__ == (float)0.) {
	    *u = rx * 1.12837916709551257389615890312154517;
	    *v = ry * 1.12837916709551257389615890312154517;
	} else {
	    *u = sx * 1.12837916709551257389615890312154517;
	    *v = sy * 1.12837916709551257389615890312154517;
	}

	if (yabs == (float)0.) {
/* Computing 2nd power */
	    d__1 = xabs;
	    *u = exp(-(d__1 * d__1));
	}

    }



/*  EVALUATION OF W(Z) IN THE OTHER QUADRANTS */


    if (*yi < (float)0.) {

	if (a) {
	    u2 *= 2;
	    v2 *= 2;
	} else {
	    xquad = -xquad;


/*         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2) */
/*         AGAINST OVERFLOW */

	    if (yquad > 3537118876014220. || xquad > 708.503061461606) {
		goto L100;
	    }

	    w1 = exp(xquad) * 2;
	    u2 = w1 * cos(yquad);
	    v2 = -w1 * sin(yquad);
	}

	*u = u2 - *u;
	*v = v2 - *v;
	if (*xi > (float)0.) {
	    *v = -(*v);
	}
    } else {
	if (*xi < (float)0.) {
	    *v = -(*v);
	}
    }

    return 0;

L100:
    *flag__ = TRUE_;
    return 0;

} /* wofz_ */
