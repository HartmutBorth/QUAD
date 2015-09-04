F90=gfortran
OPTS=-O2 -ffpe-trap=invalid,zero,overflow  -finit-real=snan
DEBOPTS=-g -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -fcheck=all -finit-real=snan

%.o : %.c
	$(CC) -c $<

%.o : %.f90
	$(F90) -c $(OPTS) $<

quad.x:	quad.o fftfor.o fftasm.o
	$(F90) -o quad.x $(OPTS) quad.o fftfor.o fftasm.o
quad.o:	quad.f90
fftfor.o:	fftfor.f90

fftasm.o:	fftasm.s
	as -o fftasm.o fftasm.s

clean:
	rm -f quad.x *.o *.mod

cleanall:
	rm -f quad.x *.o *.mod quad_out quad_out.nc quad_tseri quad_cfl quad_diag

run:
	time quad.x

nc:	
	cdo -f nc copy quad_out quad_out.nc

view:
	ncview quad_out.nc
all:
	make cleanall
	make
	make run
	make nc
	make view
