
This is a translation of part of SLATEC [http://www.netlib.org/slatec/], enough to compute cylindrical Bessel functions of complex arguments.

The main reason for the translation is that the original Fortran 77 routines are not thread-safe. I simply edited the f2c [http://www.netlib.org/f2c/] translation to remove mutable static data and other unsafe patterns. However, the result does *not* depend on f2c anymore. I hacked away the dependence so that I could call the library with const-correct arguments.

I include both the original Fortran source and the translation, and a test that runs against both. The Fortran version fails the threaded test, which is expected. The prototypes in `slatec/f2c/slatec.hpp` can be used with either the Fortran or the C++ version of the library. Other than adding `const` in places, I have not changed the interfaces in any way. All the exports are `extern "C"`.

Bear in mind that:

* The translation requires -std=c++1z on a newish compiler. There's no reason the translation couldn't be, say, C99, or an older version of C++. -std=c++1z is just what I use.
* All the error reporting routines are broken, because I'm too lazy to plunder the right pieces from f2c or to fix the f2c'd I/O. This means that if you call the routines with the wrong arguments, you'll probably get a NaN out, or your program will simply exit. Until this is fixed, it might be useful to link against the Fortran version of the library for development.
* The SLATEC functions `d1mach i1mach dlamch` are *not* exported. The are available (with value arguments) in the header `slatec/f2c/mach.hpp`.

To run the tests, just run `make`.

The whole thing is in the public domain, please feel free to fork or whatever.
