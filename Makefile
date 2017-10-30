
# Part of slatec-bessel-c++
# Makefile for both versions of the library, plus tests

SRCS_FORTRAN := $(addprefix slatec/fortran/, algo-680-erf.f d9lgmc.f dbesi1.f dbesj.f dbsi0e.f dbsknu.f dgamln.f dlngam.f initds.f xercnt.f xersve.f zacon.f zbesi.f zbinu.f zdiv.f zmlri.f zseri.f zunhj.f zunk1.f d1mach.f dasyik.f dbesi.f dbesk0.f dbesy0.f dbsi1e.f dbsynu.f dgamma.f dyairy.f j4save.f xerhlt.f xgetua.f zairy.f zbesj.f zbknu.f zexp.f zmlt.f zshch.f zuni1.f zunk2.f d9b0mp.f dasyjy.f dbesj0.f dbesk1.f dbesy1.f dbsk0e.f dcsevl.f djairy.f fdump.f lsame.f xermsg.f zabs.f zasyi.f zbesk.f zbuni.f zkscl.f zrati.f zsqrt.f zuni2.f zuoik.f d9b1mp.f dbesi0.f dbesj1.f dbesk.f dbesy.f dbsk1e.f dgamlm.f dlamch.f i1mach.f mach.f xerprn.f zacai.f zbesh.f zbesy.f zbunk.f zlog.f zs1s2.f zuchk.f zunik.f zwrsk.f)

SRCS_F2C := $(addprefix slatec/f2c/, algo-680-erf.cpp dasyjy.cpp dbesj1.cpp dbsi1e.cpp dcsevl.cpp dlngam.cpp lsame.cpp xerhlt.cpp zabs.cpp zbesh.cpp zbinu.cpp zexp.cpp zrati.cpp zuchk.cpp zunk1.cpp d9b0mp.cpp dbesi0.cpp dbesj.cpp dbesy0.cpp dbsk0e.cpp dgamlm.cpp dyairy.cpp xermsg.cpp zacai.cpp zbesi.cpp zbknu.cpp zkscl.cpp zs1s2.cpp zunhj.cpp zunk2.cpp d9b1mp.cpp dbesi1.cpp dbesk0.cpp dbesy1.cpp dbsk1e.cpp dgamln.cpp fdump.cpp xerprn.cpp zacon.cpp zbesj.cpp zbuni.cpp zlog.cpp zseri.cpp zuni1.cpp zuoik.cpp d9lgmc.cpp dbesi.cpp dbesk1.cpp dbesy.cpp dbsknu.cpp dgamma.cpp initds.cpp xersve.cpp zairy.cpp zbesk.cpp zbunk.cpp zmlri.cpp zshch.cpp zuni2.cpp zwrsk.cpp dasyik.cpp dbesj0.cpp dbesk.cpp dbsi0e.cpp dbsynu.cpp djairy.cpp j4save.cpp xercnt.cpp xgetua.cpp zasyi.cpp zbesy.cpp zdiv.cpp zmlt.cpp zsqrt.cpp zunik.cpp)

OBJS_FORTRAN := $(SRCS_FORTRAN:.f=.o)
OBJS_F2C := $(SRCS_F2C:.cpp=.o)

# https://stackoverflow.com/a/25817631
print-% : ; @echo $* = $($*)

CXXFLAGS = -std=c++1z -Wno-parentheses -O2 -fopenmp -I .
FORTRANFLAGS = -O2
LDFLAGS_FORTRAN = -lgfortran -fopenmp
LDFLAGS_F2C = -fopenmp

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o : %.f
	$(FORTRAN) $(FORTRANFLAGS) -c $< -o $@

all: test-slatec-f2c test-slatec-fortran

HEADERS = $(addprefix slatec/f2c/, slatec.hpp slatec-internal.hpp mach.hpp)

libslatec-fortran.a: $(OBJS_FORTRAN)
	$(AR) rcs $@ $(OBJS_FORTRAN)

libslatec-f2c.a: $(OBJS_F2C) $(HEADERS)
	$(AR) rcs $@ $(OBJS_F2C)

test-slatec-fortran: libslatec-fortran.a test/test-slatec.o
	$(CXX) $(LDFLAGS_FORTRAN) -o $@ test/test-slatec.o -L. -lslatec-fortran

test-slatec-f2c: libslatec-f2c.a test/test-slatec.o
	$(CXX) $(LDFLAGS_F2C) -o $@ test/test-slatec.o -L. -lslatec-f2c

all: test-slatec-fortran test-slatec-f2c
	printf "\ntest-slatec-fortran\n-------------\n" && ./test-slatec-fortran; \
	printf "\ntest-slatec-f2c\n-------------\n" && ./test-slatec-f2c

.PHONY: all
