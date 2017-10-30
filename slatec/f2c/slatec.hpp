
// SLATEC exports.
// Decls that don't seem to have use outside SLATEC have been collected in slatec-internal.hpp.

#pragma once

extern "C" {

int d9b0mp_(double const *x, double *ampl, double *theta);

int d9b1mp_(double const *x, double *ampl, double *theta);

double d9lgmc_(double const *x);

int dasyik_(double const *x, double const *fnu, int const *kode, double *flgik,
            double *ra, double *arg, int *in, double *y);

int dasyjy_(int (*funjy)(double *, double *, double *, double *, double *),
            double const *x, double const *fnu, double *flgjy, int *in, double *y,
            double *wk, int *iflw);

int dbesi_(double const *x, double const *alpha, int const *kode, int const *n, double *y, int *nz);

double dbesi0_(double const *x);

double dbesi1_(double const *x);

int dbesj_(double const *x, double const *alpha, int const *n, double *y, int *nz);

double dbesj0_(double const *x);

double dbesj1_(double const *x);

int dbesk_(double const *x, double const *fnu, int const *kode, int const *n, double *y, int *nz);

double dbesk0_(double const *x);

double dbesk1_(double const *x);

int dbesks_(double const *xnu, double const *x, int *nin, double *bk);

int dbesy_(double const *x, double const *fnu, int const *n, double *y);

double dbesy0_(double const *x);

double dbesy1_(double const *x);

double dbsi0e_(double const *x);

double dbsi1e_(double const *x);

double dbsk0e_(double const *x);

double dbsk1e_(double const *x);

int dbsknu_(double const *x, double const *fnu, int const *kode, int const *n, double *y, int *nz);

int dbsynu_(double const *x, double const *fnu, int const *n, double *y);

double dcsevl_(double *x, double const *cs, int const *n);

int dgamlm_(double *xmin, double *xmax);

double dgamln_(double const *z__, int *ierr);

double dgamma_(double const *x);

int djairy_(double *x, double *rx, double *c__, double *ai, double *dai);

double dlngam_(double const *x);

int dyairy_(double *x, double *rx, double *c__, double *bi, double *dbi);

int initds_(double const *os, int const *nos, float const *eta);

double zabs_(double const *zr, double const *zi);

int zacai_(double *zr, double *zi, double const *fnu,
           int const *kode, int *mr, int const *n, double *yr, double *
           yi, int *nz, double *rl, double *tol, double *elim,
           double *alim);

int zacon_(double *zr, double *zi, double const *fnu,
           int const *kode, int *mr, int const *n, double *yr, double *
           yi, int *nz, double *rl, double *fnul, double *tol,
           double *elim, double *alim);

int zairy_(double *zr, double *zi, int const *id,
           int const *kode, double *air, double *aii, int *nz, int
           *ierr);

int zasyi_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *
           nz, double *rl, double *tol, double *elim, double *
           alim);

int zbesh_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *m, int const *n, double *cyr, double *
           cyi, int *nz, int *ierr);

int zbesi_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *cyr, double *cyi, int *
           nz, int *ierr);

int zbesj_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *cyr, double *cyi, int *
           nz, int *ierr);

int zbesk_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *cyr, double *cyi, int *
           nz, int *ierr);

int zbesy_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *cyr, double *cyi, int *
           nz, double *cwrkr, double *cwrki, int *ierr);

int zbinu_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *cyr, double *cyi, int *
           nz, double *rl, double *fnul, double *tol, double *
           elim, double *alim);

int zbknu_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *
           nz, double *tol, double *elim, double *alim);

int zbuni_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *nz, int *nui,
           int *nlast, double *fnul, double *tol, double *elim, double *alim);

int zbunk_(double *zr, double *zi, double const *fnu,
           int const *kode, int *mr, int const *n, double *yr, double *yi, int *nz,
           double *tol, double *elim, double *alim);

int zdiv_(double const *ar, double const *ai, double const *br,
          double const *bi, double *cr, double *ci);

int zexp_(double *ar, double *ai, double *br, double *bi);

int zkscl_(double *zrr, double *zri, double const *fnu,
           int const *n, double *yr, double *yi, int *nz, double *rzr, double *rzi,
           double *ascle, double *tol, double *elim);

int zlog_(double *ar, double *ai, double *br, double *bi, int *ierr);

int zmlri_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *nz, double *tol);

int zmlt_(double *ar, double *ai, double *br,
          double *bi, double *cr, double *ci);

int zrati_(double *zr, double *zi, double const *fnu,
           int const *n, double *cyr, double *cyi, double *tol);

int zs1s2_(double *zrr, double *zri, double *s1r,
           double *s1i, double *s2r, double *s2i, int *nz,
           double *ascle, double *alim, int *iuf);

int zseri_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *nz, double *tol,
           double *elim, double *alim);

int zshch_(double *zr, double *zi, double *cshr, double *cshi, double *cchr, double *cchi);

int zsqrt_(double *ar, double *ai, double *br, double *bi);

int zuchk_(double *yr, double *yi, int *nz, double *ascle, double *tol);

int zunhj_(double *zr, double *zi, double const *fnu,
           int const *ipmtr, double *tol, double *phir, double *phii,
           double *argr, double *argi, double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i,
           double *asumr, double *asumi, double *bsumr, double *bsumi);

int zuni1_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *nz, int *nlast,
           double *fnul, double *tol, double *elim, double *alim);

int zuni2_(double *zr, double *zi, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *nz, int *nlast,
           double *fnul, double *tol, double *elim, double *alim);

int zunik_(double *zrr, double *zri, double const *fnu,
           int const *ikflg, int const *ipmtr, double *tol, int *init,
           double *phir, double *phii, double *zeta1r, double *zeta1i,
           double *zeta2r, double *zeta2i, double *sumr,
           double *sumi, double *cwrkr, double *cwrki);

int zunk1_(double *zr, double *zi, double const *fnu,
           int const *kode, int *mr, int const *n, double *yr, double *yi, int *nz,
           double *tol, double *elim, double *alim);

int zunk2_(double *zr, double *zi, double const *fnu,
           int const *kode, int *mr, int const *n, double *yr, double *yi, int *nz,
           double *tol, double *elim, double *alim);

int zuoik_(double const *zr, double const *zi, double const *fnu,
           int const *kode, int const *ikflg, int const *n, double *yr, double *yi,
           int *nuf, double *tol, double *elim, double *alim);

int zwrsk_(double *zrr, double *zri, double const *fnu,
           int const *kode, int const *n, double *yr, double *yi, int *nz,
           double *cwr, double *cwi, double *tol, double *elim, double *alim);

/* From TOMS 680, not originally in SLATEC */
int wofz_(double *xi, double *yi, double *u, double *v, int *flag__);

} // extern "C"
