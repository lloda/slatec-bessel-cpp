
#pragma once
#include <float.h>
#include <cassert>

/* dlamch and slamch were sourced from Piotr Luszczek at
 * http://icl.cs.utk.edu/lapack-forum/archives/lapack/msg01682.html */

constexpr int FLT_DIGITS = 24;
constexpr int DBL_DIGITS = 53;

constexpr inline
float
slamch_(char const ch)
{
    constexpr float one = 1.0;
    constexpr float zero = 0.0;

    if ('B' == ch || 'b' == ch) {
        return FLT_RADIX;
    } else if ('E' == ch || 'e' == ch) {
        return FLT_EPSILON;
    } else if ('L' == ch || 'l' == ch) {
        return FLT_MAX_EXP;
    } else if ('M' == ch || 'm' == ch) {
        return FLT_MIN_EXP;
    } else if ('N' == ch || 'n' == ch) {
        return FLT_DIGITS;
    } else if ('O' == ch || 'o' == ch) {
        return FLT_MAX;
    } else if ('P' == ch || 'p' == ch) {
        return FLT_EPSILON * FLT_RADIX;
    } else if ('R' == ch || 'r' == ch) {
        return FLT_ROUNDS < 2 ? one : zero;
    } else if ('S' == ch || 's' == ch) {
        /* Use SMALL plus a bit, to avoid the possibility of rounding
           causing overflow
           when computing  1/sfmin. */
        float sfmin = FLT_MIN;
        float small = one / FLT_MAX;
        if (small >= sfmin) sfmin = small * (one + FLT_EPSILON);
        return sfmin;
    } else if ('U' == ch || 'u' == ch) {
        return FLT_MIN;
    }

    return zero;
}

constexpr inline
double
dlamch_(char const ch)
{
    constexpr double one = 1.0;
    constexpr double zero = 0.0;

    if ('B' == ch || 'b' == ch) {
        return FLT_RADIX;
    } else if ('E' == ch || 'e' == ch) {
        return DBL_EPSILON;
    } else if ('L' == ch || 'l' == ch) {
        return DBL_MAX_EXP;
    } else if ('M' == ch || 'm' == ch) {
        return DBL_MIN_EXP;
    } else if ('N' == ch || 'n' == ch) {
        return DBL_DIGITS;
    } else if ('O' == ch || 'o' == ch) {
        return DBL_MAX;
    } else if ('P' == ch || 'p' == ch) {
        return DBL_EPSILON * FLT_RADIX;
    } else if ('R' == ch || 'r' == ch) {
        return FLT_ROUNDS < 2 ? one : zero;
    } else if ('S' == ch || 's' == ch) {
        /* Use SMALL plus a bit, to avoid the possibility of rounding
           causing overflow
           when computing  1/sfmin. */
        double sfmin = DBL_MIN;
        double small = one / DBL_MAX;
        if (small >= sfmin) sfmin = small * (one + DBL_EPSILON);
        return small;
    } else if ('U' == ch || 'u' == ch) {
        return DBL_MIN;
    }

    return zero;
}

constexpr inline
double
d1mach_(int const i__)
{
/* ***BEGIN PROLOGUE  D1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   D1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        D = D1MACH(I) */

/*   where I=1,...,5.  The (output) value of D above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*   D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*   D1MACH( 5) = LOG10(B) */

/*   Assume double precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(14) = T, the number of base-B digits. */
/*   I1MACH(15) = EMIN, the smallest exponent E. */
/*   I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   900911  Added SUN 386i constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/* ***END PROLOGUE  D1MACH */




/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 / */
/*     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 / */
/*     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 / */
/*     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA SMALL(1) / ZC00800000 / */
/*     DATA SMALL(2) / Z000000000 / */
/*     DATA LARGE(1) / ZDFFFFFFFF / */
/*     DATA LARGE(2) / ZFFFFFFFFF / */
/*     DATA RIGHT(1) / ZCC5800000 / */
/*     DATA RIGHT(2) / Z000000000 / */
/*     DATA DIVER(1) / ZCC6800000 / */
/*     DATA DIVER(2) / Z000000000 / */
/*     DATA LOG10(1) / ZD00E730E7 / */
/*     DATA LOG10(2) / ZC77800DC0 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O0000000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O0007777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O7770000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O7777777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA SMALL(1) / Z"3001800000000000" / */
/*     DATA SMALL(2) / Z"3001000000000000" / */
/*     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA LARGE(2) / Z"4FFE000000000000" / */
/*     DATA RIGHT(1) / Z"3FD2800000000000" / */
/*     DATA RIGHT(2) / Z"3FD2000000000000" / */
/*     DATA DIVER(1) / Z"3FD3800000000000" / */
/*     DATA DIVER(2) / Z"3FD3000000000000" / */
/*     DATA LOG10(1) / Z"3FFF9A209A84FBCF" / */
/*     DATA LOG10(2) / Z"3FFFF7988F8959AC" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */
/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */
/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */
/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */
/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn OR -pd8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CC0000000000000' / */
/*     DATA DMACH(4) / Z'3CD0000000000000' / */
/*     DATA DMACH(5) / Z'3FF34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F900000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F910000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE CRAY */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */
/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */
/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */
/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */
/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC DMACH(5) */

/*     DATA SMALL /    20K, 3*0 / */
/*     DATA LARGE / 77777K, 3*177777K / */
/*     DATA RIGHT / 31420K, 3*0 / */
/*     DATA DIVER / 32020K, 3*0 / */
/*     DATA LOG10 / 40423K, 42023K, 50237K, 74776K / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA DMACH(1) / '0000000000000010'X / */
/*     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X / */
/*     DATA DMACH(3) / '0000000000003CC0'X / */
/*     DATA DMACH(4) / '0000000000003CD0'X / */
/*     DATA DMACH(5) / '79FF509F44133FF3'X / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FORMAT */

/*     DATA DMACH(1) / '0010000000000000'X / */
/*     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X / */
/*     DATA DMACH(3) / '3CA0000000000000'X / */
/*     DATA DMACH(4) / '3CB0000000000000'X / */
/*     DATA DMACH(5) / '3FD34413509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/ */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/ */
/*     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/ */
/*     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/ */
/*     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/ */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION) */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '37777577 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000334 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '10237777 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     THREE WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 / */
/*     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B / */
/*     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B / */
/*     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) /  40000B,       0 / */
/*     DATA SMALL(3), SMALL(4) /       0,       1 / */
/*     DATA LARGE(1), LARGE(2) /  77777B, 177777B / */
/*     DATA LARGE(3), LARGE(4) / 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) /  40000B,       0 / */
/*     DATA RIGHT(3), RIGHT(4) /       0,    225B / */
/*     DATA DIVER(1), DIVER(2) /  40000B,       0 / */
/*     DATA DIVER(3), DIVER(4) /       0,    227B / */
/*     DATA LOG10(1), LOG10(2) /  46420B,  46502B / */
/*     DATA LOG10(3), LOG10(4) /  76747B, 176377B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B / */
/*     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B / */
/*     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B / */
/*     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B / */
/*     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTANTS FOR THE IBM PC */
/*     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION */
/*     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087. */

/*     DATA SMALL(1) / 2.23D-308  / */
/*     DATA LARGE(1) / 1.79D+308  / */
/*     DATA RIGHT(1) / 1.11D-16   / */
/*     DATA DIVER(1) / 2.22D-16   / */
/*     DATA LOG10(1) / 0.301029995663981195D0 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    8388608,           0 / */
/*     DATA LARGE(1), LARGE(2) / 2147483647,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /  612368384,           0 / */
/*     DATA DIVER(1), DIVER(2) /  620756992,           0 / */
/*     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 / */

/*     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 / */
/*     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 / */
/*     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 / */
/*     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    128,      0 / */
/*     DATA SMALL(3), SMALL(4) /      0,      0 / */
/*     DATA LARGE(1), LARGE(2) /  32767,     -1 / */
/*     DATA LARGE(3), LARGE(4) /     -1,     -1 / */
/*     DATA RIGHT(1), RIGHT(2) /   9344,      0 / */
/*     DATA RIGHT(3), RIGHT(4) /      0,      0 / */
/*     DATA DIVER(1), DIVER(2) /   9472,      0 / */
/*     DATA DIVER(3), DIVER(4) /      0,      0 / */
/*     DATA LOG10(1), LOG10(2) /  16282,   8346 / */
/*     DATA LOG10(3), LOG10(4) / -31493, -12296 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA SMALL(3), SMALL(4) / O000000, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA LARGE(3), LARGE(4) / O177777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O022200, O000000 / */
/*     DATA RIGHT(3), RIGHT(4) / O000000, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O022400, O000000 / */
/*     DATA DIVER(3), DIVER(4) / O000000, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020232 / */
/*     DATA LOG10(3), LOG10(4) / O102373, O147770 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE SUN 386i */

/*     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' / */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' / */
/*     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF' */
/*     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 / */

/* ***FIRST EXECUTABLE STATEMENT  D1MACH */

    constexpr double dmach[5] = { dlamch_('u'), dlamch_('o'), dlamch_('e'), dlamch_('p'), log10(dlamch_('b')) };
    assert(i__>=1 && i__<=5 && "Bad argument for d1mach");
    return dmach[i__ - 1];
}

constexpr inline
int
i1mach_(int const i__)
{
/* ***BEGIN PROLOGUE  I1MACH */
/* ***PURPOSE  Return integer machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (I1MACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   I1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument and can be referenced as follows: */

/*        K = I1MACH(I) */

/*   where I=1,...,16.  The (output) value of K above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   I/O unit numbers: */
/*     I1MACH( 1) = the standard input unit. */
/*     I1MACH( 2) = the standard output unit. */
/*     I1MACH( 3) = the standard punch unit. */
/*     I1MACH( 4) = the standard error message unit. */

/*   Words: */
/*     I1MACH( 5) = the number of bits per integer storage unit. */
/*     I1MACH( 6) = the number of characters per integer storage unit. */

/*   Integers: */
/*     assume integers are represented in the S-digit, base-A form */

/*                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*                where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*     I1MACH( 7) = A, the base. */
/*     I1MACH( 8) = S, the number of base-A digits. */
/*     I1MACH( 9) = A**S - 1, the largest magnitude. */

/*   Floating-Point Numbers: */
/*     Assume floating-point numbers are represented in the T-digit, */
/*     base-B form */
/*                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*                where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*     I1MACH(10) = B, the base. */

/*   Single-Precision: */
/*     I1MACH(11) = T, the number of base-B digits. */
/*     I1MACH(12) = EMIN, the smallest exponent E. */
/*     I1MACH(13) = EMAX, the largest exponent E. */

/*   Double-Precision: */
/*     I1MACH(14) = T, the number of base-B digits. */
/*     I1MACH(15) = EMIN, the smallest exponent E. */
/*     I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891012  Added VAX G-floating constants.  (WRB) */
/*   891012  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. */
/*           (RWC) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added Convex -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler */
/*           options.  (DWL, RWC and WRB). */
/* ***END PROLOGUE  I1MACH */


/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        129 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1025 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA IMACH( 1) /          7 / */
/*     DATA IMACH( 2) /          2 / */
/*     DATA IMACH( 3) /          2 / */
/*     DATA IMACH( 4) /          2 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         33 / */
/*     DATA IMACH( 9) / Z1FFFFFFFF / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -256 / */
/*     DATA IMACH(13) /        255 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /       -256 / */
/*     DATA IMACH(16) /        255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /        -50 / */
/*     DATA IMACH(16) /         76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /     -32754 / */
/*     DATA IMACH(16) /      32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -4095 / */
/*     DATA IMACH(13) /       4094 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -4095 / */
/*     DATA IMACH(16) /       4094 / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /    6LOUTPUT/ */
/*     DATA IMACH( 5) /         60 / */
/*     DATA IMACH( 6) /         10 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         48 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /       -929 / */
/*     DATA IMACH(13) /       1070 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /       -929 / */
/*     DATA IMACH(16) /       1069 / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / Z'7FFFFFFF' / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16383 / */
/*     DATA IMACH(16) /      16383 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -pd8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 46 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         46 / */
/*     DATA IMACH( 9) / 1777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 64 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 777777777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     DATA IMACH( 1) /         11 / */
/*     DATA IMACH( 2) /         12 / */
/*     DATA IMACH( 3) /          8 / */
/*     DATA IMACH( 4) /         10 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         24 / */
/*     DATA IMACH( 6) /          3 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         23 / */
/*     DATA IMACH( 9) /    8388607 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         38 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /         43 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         63 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         39 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         55 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          7 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1015 / */
/*     DATA IMACH(16) /       1017 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) /  Z7FFFFFFF / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */


/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         54 / */
/*     DATA IMACH(15) /       -101 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         62 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1021 / */
/*     DATA IMACH(13) /       1024 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16381 / */
/*     DATA IMACH(16) /      16384 / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          1 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /      -1024 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA IMACH( 1) /          1 / */
/*     DATA IMACH( 2) /          1 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/* ***FIRST EXECUTABLE STATEMENT  I1MACH */

    constexpr int imach[16] = { 5, 6, 0, 0, 32, 4, 2, 31, 2147483647, 2, 24, -125, 127,
                                53, -1021, 1023 };
    assert(i__>=1 && i__<=16);
    return imach[i__ - 1];
}
