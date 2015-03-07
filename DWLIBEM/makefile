all: DWNLIBEM

DWNLIBEM : DWNLIBEM.f90
	gfortran -o $@ DWNLIBEM.f90 -Wall -g -cpp -I/usr/local/include -L/usr/local/lib -lfftw3 -I/usr/local/dislin/gf -L/opt/X11/lib -lXt -lXm -L/usr/OpenMotif/lib -L/usr/local/dislin -ldislin -lm -L../preinverted -lLHS -framework Accelerate