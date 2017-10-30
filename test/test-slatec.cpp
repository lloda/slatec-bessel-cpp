
// Part of slatec-bessel-c++
// Check of SLATEC used under OpenMP, to avoid regression.

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <chrono>
#include "slatec/f2c/slatec.hpp"

using std::cout, std::endl, std::flush;

template <class ... A> inline std::string
format(A && ... a)
{
    std::ostringstream o; (o << ... << a); return o.str();
}

template <class X> inline auto sqr(X && x) { return x*x; }

using clck = std::conditional_t<std::chrono::high_resolution_clock::is_steady,
                                std::chrono::high_resolution_clock,
                                std::chrono::steady_clock>;

static double
toseconds(clck::duration const & t)
{
    return std::chrono::duration<float, std::ratio<1, 1>>(t).count();
}

int main()
{
    int const c_1 = 1;
    double const r_1 = 1.;

    int const N = 1000000;

    std::vector<double> x(N), y(N);
    std::mt19937 random_sequence(99);
    std::uniform_real_distribution<double> xydist(-1, 1);
    x[0] = 3.4;
    y[0] = 0.;
    for (unsigned i=1; i<x.size(); ++i) {
        x[i] = xydist(random_sequence);
        y[i] = xydist(random_sequence);
    }

// threaded
    std::vector<double> re(N), im(N);
    {
        auto t0 = clck::now();
#pragma omp parallel for
        for (unsigned i=0; i<x.size(); ++i) {
            int nz, ierr;
            zbesj_(&x[i], &y[i], &r_1, &c_1, &c_1, &re[i], &im[i], &nz, &ierr);
        }
        clck::duration t1 = clck::now()-t0;
        cout << "threaded " << toseconds(t1)*1e3 << "ms. " << endl;
    }

// serial
    std::vector<double> re_ref(N), im_ref(N);
    {
        auto t0 = clck::now();
        for (unsigned i=0; i<x.size(); ++i) {
            int nz, ierr;
            zbesj_(&x[i], &y[i], &r_1, &c_1, &c_1, &re_ref[i], &im_ref[i], &nz, &ierr);
        }
        clck::duration t1 = clck::now()-t0;
        cout << "serial " << toseconds(t1)*1e3 << "ms. " << endl;
    }

// check
    int failures = 0;
    for (unsigned i=0; i<x.size(); ++i) {
        double err = (sqr(re[i]-re_ref[i]) + sqr(im[i]-im_ref[i]));
        if (err!=0) {
            cout << format("error at ", i, " : ", err) << endl;
            ++failures;
        }
    }

// check vs alt impl.
    {
        double err = sqrt(sqr(re[0]-std::cyl_bessel_j(1., x[0])) + sqr(im[0]-0.));
        if (!(err<=1e-15)) {
            cout << format("error omp vs alt : ", err) << endl;
            ++failures;
        }
    }
    {
        double err = sqrt(sqr(re_ref[0]-std::cyl_bessel_j(1., x[0])) + sqr(im_ref[0]-0.));
        if (!(err<=1e-15)) {
            cout << format("error serial vs alt : ", err) << endl;
            ++failures;
        }
    }

// other misc checks
    {
        double const r_07 = 0.7;
        {
            double err = std::abs(dgamma_(&r_07)-std::tgamma(0.7));
            if (!(err<=1e-15)) {
                cout << format("error dgamma_ : ", err) << endl;
                ++failures;
            }
        }
        {
            int ierr = 99;
            double err = std::abs(dgamln_(&r_07, &ierr)-std::lgamma(0.7));
            if (!(err<=1e-15) || ierr!=0) {
                cout << format("error dgamln_ : ", err, " ierr: ", ierr) << endl;
                ++failures;
            }
        }
        double const r_m11 = -11.1;
        {
            double a = dlngam_(&r_m11);
            double b = std::lgamma(-11.1);
            double err = std::abs(a-b);
            cout << format("error dlngam_ a: ", a, " vs b: ", b, ", err: ", err) << endl;
            if (!(err<=2e-14)) {
                ++failures;
            }
        }
    }

    cout << "total " << failures << " failures." << endl;
    return failures;
};
