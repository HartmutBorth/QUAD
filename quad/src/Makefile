F90=gfortran
OPTS=-O2 -ffpe-trap=invalid,zero,overflow  -finit-real=snan
DEBOPTS=-g -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -fcheck=all -finit-real=snan

%.o : %.c
	$(CC) -c $<

%.o : %.f90
	$(F90) -c $(OPTS) $<

quad.x:	quad.o fftmod.o
	$(F90) -o quad.x $(OPTS) quad.o fftmod.o
quad.o:	quad.f90
fftmod.o:	fftmod.f90

clean:
	rm -f quad.x *.o *.mod
