# compilador
FF = gfortran
DE = lldb
FFLAGSDEB = -Wall -g -cpp -fbounds-check
FFLAGSOP = -O3 -cpp
DEPS = -I/usr/local/include -L/usr/local/lib -lfftw3 -I/usr/local/dislin/gf -L/opt/X11/lib -lXt -lXm -L/usr/OpenMotif/lib -L/usr/local/dislin -ldislin -lm -framework Accelerate


DWNLIBEM.out : DWNLIBEM.F90 cine.a libWave.a
	$(FF) -o $@ $< libWave.a cine.a $(FFLAGSDEB) $(DEPS)

cine.a: cine.f90 libWave.a
	$(FF) -O3 -c cine.f90 -lWave $(DEPS)
	ar cr cine.a cine.o
	rm *.o

libWave.a: libWave.f90
	$(FF) -O3 -ffast-math -c libWave.f90 $(DEPS)
	ar cr libWave.a libWave.o
	rm *.o

DEB: DWNLIBEM.out
	$(DE) -f $< -s DEB.sh


clean:
	rm *.mod
	rm *.o
	rm libWave.a
	rm DWNLIBEM 

