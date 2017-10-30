
// FIXME patch the i/o functions.
// FIXME remove all remaining flag-first-init type patterns in the routines.

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "slatec.hpp"
#include "mach.hpp"

using integer = int;
using logical = int;
using ftnlen = int;
using ftnint = int;
using flag = int;
using address = char *;

constexpr logical TRUE_ = 1;
constexpr logical FALSE_ = 0;

/*internal read, write*/
struct icilist
{
    flag icierr;
    char *iciunit;
    flag iciend;
    char const *icifmt;
    ftnint icirlen;
    ftnint icirnum;
};

/*external read, write*/
struct cilist
{
    flag cierr;
    ftnint ciunit;
    flag ciend;
    char const *cifmt;
    ftnint cirec;
};

using std::log, std::sqrt, std::exp, std::sin, std::cos, std::sinh, std::cosh,
      std::abs, std::atan, std::min, std::max;

namespace f2c {

inline double d_sign(double const *a, double const *b)
{
    double x;
    x = (*a >= 0 ? *a : - *a);
    return (*b >= 0 ? x : -x);
}

inline double d_int(double const *x)
{
    return( (*x>0) ? floor(*x) : -floor(- *x) );
}

inline integer do_fio(integer const *, char const *, ftnlen)
{
    return 0; // FIXME from fmt.c, but replace
}

inline integer e_wsfe() // write sequential formatted external??
{
    return 0;
}

inline integer e_wsfi() // internal I guess
{
    return 0;
}

inline integer s_wsfe(cilist *)
{
    return 0;
}

inline integer s_wsfi(icilist *)
{
    return 0;
}

inline integer i_len(char const *s, ftnlen const n)
{
    return (n);
}

inline integer i_indx(char const *a, char const *b, ftnlen const la, ftnlen const lb)
{
    ftnlen i, n;
    char const *s;
    char const *t;
    char const *bend;

    n = la - lb + 1;
    bend = b + lb;

    for (i = 0 ; i < n ; ++i) {
        s = a + i;
        t = b;
        while(t < bend)
            if(*s++ != *t++)
                goto no;
        return(i+1);
    no: ;
    }
    return(0);
}

/* Unless compiled with -DNO_OVERWRITE, this variant of s_copy allows the
 * target of an assignment to appear on its right-hand side (contrary
 * to the Fortran 77 Standard, but in accordance with Fortran 90),
 * as in  a(2:5) = a(4:7) .
 */
/* assign strings:  a = b */
inline void s_copy(char *a, char const *b, ftnlen const la, ftnlen const lb)
{
    char *aend;
    char const *bend;

    aend = a + la;

    if(la <= lb)
#ifndef NO_OVERWRITE
        if (a <= b || a >= b + la)
#endif
            while(a < aend)
                *a++ = *b++;
#ifndef NO_OVERWRITE
        else
            for(b += la; a < aend; )
                *--aend = *--b;
#endif

    else {
        bend = b + lb;
#ifndef NO_OVERWRITE
        if (a <= b || a >= bend)
#endif
            while(b < bend)
                *a++ = *b++;
#ifndef NO_OVERWRITE
        else {
            a += lb;
            while(b < bend)
                *--a = *--bend;
            a += lb;
        }
#endif
        while(a < aend)
            *a++ = ' ';
    }
}

inline integer s_cmp(char const *a0, char const *b0, ftnlen const la, ftnlen const lb)
{
    unsigned char *a, *aend, *b, *bend;
    a = (unsigned char *)a0;
    b = (unsigned char *)b0;
    aend = a + la;
    bend = b + lb;

    if(la <= lb)
	{
            while(a < aend)
		if(*a != *b)
                    return( *a - *b );
		else
                    { ++a; ++b; }

            while(b < bend)
		if(*b != ' ')
                    return( ' ' - *b );
		else	++b;
	}

    else
	{
            while(b < bend)
		if(*a == *b)
                    { ++a; ++b; }
		else
                    return( *a - *b );
            while(a < aend)
		if(*a != ' ')
                    return(*a - ' ');
		else	++a;
	}
    return(0);
}

inline char *
F77_aloc(integer const Len, char const *whence)
{
     char *rv;
     unsigned int uLen = (unsigned int) Len;	/* for K&R C */

     if (!(rv = (char*)malloc(uLen))) {
         fprintf(stderr, "malloc(%u) failure in %s\n", uLen, whence);
         std::exit(3);
     }
     return rv;
 }

inline void s_cat(char * lp, char *rpp[], ftnint const rnp[], ftnint const *np, ftnlen ll)
{
    ftnlen i, nc;
    char const *rp;
    ftnlen n = *np;
#ifndef NO_OVERWRITE
    ftnlen L, m;
    char *lp0, *lp1;

    lp0 = 0;
    lp1 = lp;
    L = ll;
    i = 0;
    while(i < n) {
        rp = rpp[i];
        m = rnp[i++];
        if (rp >= lp1 || rp + m <= lp) {
            if ((L -= m) <= 0) {
                n = i;
                break;
            }
            lp1 += m;
            continue;
        }
        lp0 = lp;
        lp = lp1 = F77_aloc(L = ll, "s_cat");
        break;
    }
    lp1 = lp;
#endif /* NO_OVERWRITE */
    for(i = 0 ; i < n ; ++i) {
        nc = ll;
        if(rnp[i] < nc)
            nc = rnp[i];
        ll -= nc;
        rp = rpp[i];
        while(--nc >= 0)
            *lp++ = *rp++;
    }
    while(--ll >= 0)
        *lp++ = ' ';
#ifndef NO_OVERWRITE
    if (lp0) {
        std::memcpy(lp0, lp1, L);
        free(lp1);
    }
#endif
}

inline int s_stop(char const *s, ftnlen n)
{
    int i;

    if(n > 0) {
        fprintf(stderr, "STOP ");
        for(i = 0; i<n ; ++i)
            putc(*s++, stderr);
        fprintf(stderr, " statement executed\n");
    }
    exit(0);

    /* We cannot avoid (useless) compiler diagnostics here:		*/
    /* some compilers complain if there is no return statement,	*/
    /* and others complain that this one cannot be reached.		*/

    return 0; /* NOT REACHED */
}

inline double pow_dd(double const *ap, double const *bp)
{
    return(std::pow(*ap, *bp) );
}

inline double pow_di(double const *ap, integer const *bp)
{
    double pow, x;
    integer n;
    unsigned long u;

    pow = 1;
    x = *ap;
    n = *bp;

    if(n != 0) {
        if(n < 0) {
            n = -n;
            x = 1/x;
        }
        for(u = n; ; ) {
            if(u & 01)
                pow *= x;
            if(u >>= 1)
                x *= x;
            else
                break;
        }
    }
    return(pow);
}

inline int i_dnnt(double const *x)
{
    return (integer)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}

} // namespace f2c

extern "C" {

int xercnt_(char *librar, char *subrou, char *messg, integer
            *nerr, integer *level, integer *kontrl, ftnlen librar_len, ftnlen
            subrou_len, ftnlen messg_len);

int xerhlt_(char const *messg, ftnlen const messg_len);

int xermsg_(char const *librar, char const *subrou, char const *messg,
            integer const *nerr, integer const *level,
            ftnlen const librar_len, ftnlen const subrou_len, ftnlen const messg_len);

int xerprn_(char const *prefix, integer const *npref, char const *messg,
            integer const *nwrap, ftnlen prefix_len, ftnlen messg_len);

int xersve_(char const *librar, char const *subrou, char const *messg, integer
            const *kflag, integer const *nerr, integer const *level, integer *icount, ftnlen
            librar_len, ftnlen subrou_len, ftnlen messg_len);

int xgetua_(integer *iunita, integer *n);

integer j4save_(integer const *iwhich, integer const *ivalue, logical const *iset);

int fdump_();

} // extern "C"
