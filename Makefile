all: _rif

_rif: _rif.c _rif.pyf
	f2py _rif.pyf _rif.c -c -I/usr/local/include -lfftw3 -lm

clean:
	-rm *.so

