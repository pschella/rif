CFLAGS="-O3 -ffast-math"

all: _rif

_rif: _rif.c _rif.pyf
	f2py _rif.pyf _rif.c -c --opt=${CFLAGS} -I/usr/local/include -lfftw3 -lm

clean:
	-rm *.so

